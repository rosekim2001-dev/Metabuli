
#ifndef GROUP_GENERATOR_H
#define GROUP_GENERATOR_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <cassert>
#include <thread>
#include <atomic>
#include <chrono>
#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <sys/statvfs.h>
#include <dirent.h>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <zstd.h>
#include "IndexCreator.h"
#include "SeqIterator.h"
#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include "QueryIndexer.h"
#include "KmerExtractor.h"
#include "KmerMatcher.h"   // QueryKmerSplit
#include "KSeqWrapper.h"
#include "DeltaIdxReader.h"

#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;

struct Relation {
    uint32_t id1;
    uint32_t id2;
    uint16_t weight;

    Relation(uint32_t a = 0, uint32_t b = 0): id1(a), id2(b), weight(0) {}
    Relation(uint32_t a, uint32_t b, uint16_t w): id1(a), id2(b), weight(w) {}

    static bool compare(const Relation& a, const Relation& b) {
        if (a.id1 != b.id1) return a.id1 < b.id1;
        return a.id2 < b.id2;
    }

    bool operator<(const Relation& other) const {
        if (id1 != other.id1) return id1 < other.id1;
        return id2 < other.id2;
    }

    bool operator==(const Relation& other) const {
        return id1 == other.id1 && id2 == other.id2;
    }
};

// Relation is 12 bytes but only 10 carry payload; the compiler never initialises the trailing
// two. They are written to disk verbatim, so leaving them undefined makes the intermediate files
// differ byte for byte between runs even when every record is identical -- which is exactly what
// happened to relations.bin once the surrounding code changed and the stack slot stopped being
// reused. Build every record that reaches a file through here.
// The void* cast is the documented way to memset a class type without -Wclass-memaccess.
inline Relation makeRelation(uint32_t id1, uint32_t id2, uint16_t weight) {
    Relation rel;
    std::memset(static_cast<void *>(&rel), 0, sizeof(rel));
    rel.id1 = id1;
    rel.id2 = id2;
    rel.weight = weight;
    return rel;
}

// Sequential Relation streams for the subGraph_* files ONLY. relations_* keeps using
// ReadBuffer/WriteBuffer<Relation>: giving the two a different type is the point, so that a
// change to the subgraph representation cannot silently reach the files makeGroups reads.
//
// On-disk layout: one zstd stream per file, nothing else. Records go in as raw Relation bytes,
// so the two padding bytes of every 12-byte struct are zeros the compressor removes for free.
// Streaming (not one-shot) because a flush can carry 200 M records -- one-shot would need a
// 2.4 GB input and a matching output buffer, while the streaming API bounds both.
//
// Every failure is fatal. saveSubGraphToFile used to print to cerr and continue, which let a
// 21-hour run finish with a zero exit code and an empty result after 70,399 writes had failed.
static const int RELATION_ZSTD_LEVEL = 1;             // ~500 MB/s/core; emit runs near 100 MB/s
static const size_t RELATION_WRITE_ELEMS = 65'536;    // 768 KB of records per compressor feed

// How far either side of a common-k-mer hit is discarded along with the hit is --common-kmer-span,
// held in commonKmerSpan below. It was a constant of 1 until 2026-08-24, justified by the
// measurement that a span of 1 removes 17% of query k-mers but only 0.7% of edges: it lowers
// weights rather than deleting pairs, and --min-overlap-ratio scales with the surviving count. That
// measurement counts edges, and edge count is not connectivity -- the edges holding two stretches
// of one genome together are rare, so a small edge loss can still cut them and leave the genome in
// more pieces. On CAMI2 strain-madness, 1 -> 0 took k-mers per read from 73.6 to 81.4 and species
// Recall*c from 0.164 to 0.245 for 0.010 of purity, so the value is worth sweeping and is now a
// parameter.

// Phase 2 never merges two units on a single weak link. One link is what chance sharing across a
// repeat looks like, so the smallest requirement carrying any information is two.
static const uint32_t MERGE_SUPPORT_FLOOR = 2;

static inline void relationStreamFatal(const std::string & path, const char * action) {
    std::cerr << "Error " << action << " file: " << path
              << " (" << strerror(errno) << ")" << std::endl;
    exit(EXIT_FAILURE);
}

static inline void relationStreamZstdFatal(const std::string & path, const char * action,
                                           size_t code) {
    std::cerr << "Error " << action << " file: " << path
              << " (zstd: " << ZSTD_getErrorName(code) << ")" << std::endl;
    exit(EXIT_FAILURE);
}

// Fixed per-stream cost of the codec, on top of the caller's record buffer. Charged to the
// merge memory budget in getMergeBufferElems so the disk representation and the memory
// representation are not silently conflated.
inline size_t relationStreamFixedBytes() {
    return ZSTD_DStreamInSize();
}

struct ZstdCStreamDeleter {
    void operator()(ZSTD_CStream * s) const { ZSTD_freeCStream(s); }
};

struct ZstdDStreamDeleter {
    void operator()(ZSTD_DStream * s) const { ZSTD_freeDStream(s); }
};

class RelationWriter {
public:
    explicit RelationWriter(const std::string & filePath)
        : path(filePath), fp(fopen(filePath.c_str(), "wb")), cstream(ZSTD_createCStream()),
          records(), outBuf(ZSTD_CStreamOutSize()), bytesOnDisk(0) {
        if (fp == nullptr) {
            relationStreamFatal(path, "opening");
        }
        if (!cstream) {
            relationStreamFatal(path, "allocating a compressor for");
        }
        const size_t init = ZSTD_initCStream(cstream.get(), RELATION_ZSTD_LEVEL);
        if (ZSTD_isError(init)) {
            relationStreamZstdFatal(path, "initializing the compressor for", init);
        }
        records.reserve(RELATION_WRITE_ELEMS);
    }

    RelationWriter(const RelationWriter &) = delete;
    RelationWriter & operator=(const RelationWriter &) = delete;

    ~RelationWriter() { finish(); }

    void write(const Relation & rel) {
        records.push_back(rel);
        if (records.size() >= RELATION_WRITE_ELEMS) {
            feedCompressor();
        }
    }

    // Idempotent so the destructor can act as a safety net after an explicit call.
    void finish() {
        if (fp == nullptr) {
            return;
        }
        feedCompressor();

        size_t remaining = 0;
        do {
            ZSTD_outBuffer out = { outBuf.data(), outBuf.size(), 0 };
            remaining = ZSTD_endStream(cstream.get(), &out);
            if (ZSTD_isError(remaining)) {
                relationStreamZstdFatal(path, "finishing", remaining);
            }
            writeOut(out.pos);
        } while (remaining != 0);

        FILE * const closing = fp;
        fp = nullptr;
        if (fclose(closing) != 0) {
            relationStreamFatal(path, "closing");
        }
    }

    size_t compressedBytes() const { return bytesOnDisk; }

private:
    void feedCompressor() {
        if (records.empty()) {
            return;
        }
        ZSTD_inBuffer in = { records.data(), records.size() * sizeof(Relation), 0 };
        while (in.pos < in.size) {
            ZSTD_outBuffer out = { outBuf.data(), outBuf.size(), 0 };
            const size_t code = ZSTD_compressStream(cstream.get(), &out, &in);
            if (ZSTD_isError(code)) {
                relationStreamZstdFatal(path, "compressing into", code);
            }
            writeOut(out.pos);
        }
        records.clear();
    }

    void writeOut(size_t bytes) {
        if (bytes == 0) {
            return;
        }
        if (fwrite(outBuf.data(), 1, bytes, fp) != bytes) {
            relationStreamFatal(path, "writing to");
        }
        bytesOnDisk += bytes;
    }

    std::string path;
    FILE * fp;
    std::unique_ptr<ZSTD_CStream, ZstdCStreamDeleter> cstream;
    std::vector<Relation> records; // reserve() fixes the capacity, so push_back never reallocates
    std::vector<char> outBuf;
    size_t bytesOnDisk;
};

class RelationReader {
public:
    RelationReader(const std::string & filePath, size_t blockElems)
        : path(filePath), fp(fopen(filePath.c_str(), "rb")), dstream(ZSTD_createDStream()),
          inBuf(ZSTD_DStreamInSize()), inFilled(0), inPos(0),
          outBytes((blockElems == 0 ? 1 : blockElems) * sizeof(Relation)),
          outFilled(0), outPos(0), sourceDrained(false) {
        if (fp == nullptr) {
            relationStreamFatal(path, "opening");
        }
        if (!dstream) {
            relationStreamFatal(path, "allocating a decompressor for");
        }
        const size_t init = ZSTD_initDStream(dstream.get());
        if (ZSTD_isError(init)) {
            relationStreamZstdFatal(path, "initializing the decompressor for", init);
        }
    }

    RelationReader(const RelationReader &) = delete;
    RelationReader & operator=(const RelationReader &) = delete;

    ~RelationReader() {
        if (fp != nullptr) {
            fclose(fp);
        }
    }

    // Returns a default Relation at end of stream, matching ReadBuffer<Relation>'s contract so
    // the callers' `if (!(next == Relation()))` termination test still holds. Relation::operator==
    // ignores weight, so the sentinel is really (id1,id2) == (0,0) -- an edge that cannot exist,
    // because writeKmers globalizes query IDs as 1-based (sequenceID += processedReadCnt).
    Relation getNext() {
        if (outFilled - outPos < sizeof(Relation) && !fill()) {
            return Relation();
        }
        Relation rel;
        memcpy(&rel, outBytes.data() + outPos, sizeof(Relation));
        outPos += sizeof(Relation);
        return rel;
    }

private:
    // Decompress until at least one whole record is buffered. A record can straddle two
    // decompressor outputs, so the leftover tail is carried to the front instead of dropped.
    bool fill() {
        const size_t leftover = outFilled - outPos;
        if (leftover > 0 && outPos > 0) {
            memmove(outBytes.data(), outBytes.data() + outPos, leftover);
        }
        outFilled = leftover;
        outPos = 0;

        while (outFilled < sizeof(Relation)) {
            if (inPos >= inFilled) {
                if (sourceDrained) {
                    return false;
                }
                inFilled = fread(inBuf.data(), 1, inBuf.size(), fp);
                inPos = 0;
                if (inFilled == 0) {
                    sourceDrained = true;
                    return false;
                }
            }

            ZSTD_inBuffer in = { inBuf.data(), inFilled, inPos };
            ZSTD_outBuffer out = { outBytes.data(), outBytes.size(), outFilled };
            const size_t code = ZSTD_decompressStream(dstream.get(), &out, &in);
            if (ZSTD_isError(code)) {
                relationStreamZstdFatal(path, "decompressing", code);
            }
            inPos = in.pos;
            outFilled = out.pos;
        }
        return true;
    }

    std::string path;
    FILE * fp;
    std::unique_ptr<ZSTD_DStream, ZstdDStreamDeleter> dstream;
    std::vector<char> inBuf;
    size_t inFilled;
    size_t inPos;
    std::vector<char> outBytes;
    size_t outFilled;
    size_t outPos;
    bool sourceDrained;
};

class DisjointSet {
public:
    std::vector<uint32_t> parent;
    std::vector<uint32_t> rank;
    std::vector<bool> grouped;

    DisjointSet(size_t numQuery) 
        : parent(numQuery + 1), rank(numQuery + 1), grouped(numQuery + 1) {
        for (size_t i = 1; i <= numQuery; ++i) {
            parent[i] = i;
            rank[i] = 0;
            grouped[i] = false;
        }
    }

    uint32_t find(uint32_t element) {
        if (parent[element] != element) {
            parent[element] = find(parent[element]);
        }
        return parent[element];
    }

    void unionSets(uint32_t elem1, uint32_t elem2) {
        uint32_t root1 = find(elem1);
        uint32_t root2 = find(elem2);

        grouped[elem1] = true;
        grouped[elem2] = true;

        if (root1 != root2) {
            if (rank[root1] < rank[root2]) {
                parent[root1] = root2;
            } else if (rank[root1] > rank[root2]) {
                parent[root2] = root1;
            } else {
                if (root1 < root2) {
                    parent[root2] = root1;
                    rank[root1]++;
                }
                else{
                    parent[root1] = root2;
                    rank[root2]++;
                }
            }
        }
    }

    void flatten() {
        for (size_t i = 1; i < parent.size(); ++i) {
            parent[i] = find(i);
        }
    }

    DisjointSet& operator+=(const DisjointSet& rhs) {
        for (size_t i = 1; i < parent.size(); ++i) {
            uint32_t p = rhs.parent[i];
            if (p != i) unionSets(i, p);
        }
        return *this;
    }

    friend DisjointSet operator+(DisjointSet lhs, const DisjointSet& rhs) {
        lhs += rhs;
        return lhs;
    }
};

// Saturating accumulation into a uint16 edge weight: clamps at UINT16_MAX
// instead of wrapping. Relation::weight is uint16_t, so both the per-subgraph
// count and the cross-subgraph sum can exceed 65535; a wrapped weight would
// look small and be dropped by the Phase 1/2 threshold comparisons.
static inline void addSat16(uint16_t & dst, uint32_t add) {
    const uint32_t sum = static_cast<uint32_t>(dst) + add;
    dst = (sum > UINT16_MAX) ? UINT16_MAX : static_cast<uint16_t>(sum);
}

// Memory usable for buffers.
//
// Linux uses /proc/meminfo's MemAvailable, not sysinfo's freeram: freeram excludes the page
// cache even though the cache is reclaimable, so during a multi-hundred-GB write it collapses.
// That collapse used to shrink the flush threshold and split the graph into tens of thousands
// of tiny subgraph files, which then exhausted the process's file descriptors.
// macOS uses hw.memsize, which is TOTAL memory and therefore an over-estimate.
inline size_t getAvailableMemoryBytes() {
#if defined(__linux__)
    // Line-by-line: /proc/meminfo rows do not all carry a unit token, so streaming with >>
    // would read the next row's key as this row's unit.
    std::ifstream meminfo("/proc/meminfo");
    std::string line;
    while (std::getline(meminfo, line)) {
        unsigned long long kb = 0;
        if (sscanf(line.c_str(), "MemAvailable: %llu kB", &kb) == 1) {
            return static_cast<size_t>(kb) * 1024;
        }
    }
    // Pre-3.14 kernels have no MemAvailable; free + reclaimable buffers is the closest stand-in.
    struct sysinfo info;
    sysinfo(&info);
    return (size_t)(info.freeram + info.bufferram) * info.mem_unit;
#elif defined(__APPLE__)
    int64_t freeMemory;
    size_t len = sizeof(freeMemory);
    sysctlbyname("hw.memsize", &freeMemory, &len, nullptr, 0);
    return (size_t)freeMemory; // 근사치
#else
    return 8ULL * 1024 * 1024 * 1024; // fallback 8GB
#endif
}

// Memory budget: what the OS can spare, capped by the user's declared --max-ram. The rest of
// the pipeline already budgets against --max-ram (Buffer::calculateBufferSize); these two
// helpers used to ignore it and trust an instantaneous OS reading alone.
inline size_t getMemoryBudgetBytes(int maxRamGiB) {
    const size_t availableBytes = getAvailableMemoryBytes();
    if (maxRamGiB <= 0) {
        return availableBytes;
    }
    const size_t declared = static_cast<size_t>(maxRamGiB) * 1024 * 1024 * 1024;
    return std::min(availableBytes, declared);
}

// Free space at `path`, or 0 when it cannot be measured. Same platform-guard shape as
// getAvailableMemoryBytes: a value we cannot trust is reported as "unknown", not guessed.
inline size_t getFreeDiskBytes(const std::string & path) {
#if defined(__linux__) || defined(__APPLE__)
    struct statvfs info;
    if (statvfs(path.c_str(), &info) != 0) {
        return 0;
    }
    return static_cast<size_t>(info.f_bavail) * static_cast<size_t>(info.f_frsize);
#else
    (void) path;
    return 0;
#endif
}

// Returned when no cap applies -- either the user asked for none or free space is unknown.
static const size_t TMP_DISK_UNLIMITED = SIZE_MAX;

// Total capacity of the filesystem holding `path`, or 0 when it cannot be measured.
inline size_t getTotalDiskBytes(const std::string & path) {
#if defined(__linux__) || defined(__APPLE__)
    struct statvfs info;
    if (statvfs(path.c_str(), &info) != 0) {
        return 0;
    }
    return static_cast<size_t>(info.f_blocks) * static_cast<size_t>(info.f_frsize);
#else
    (void) path;
    return 0;
#endif
}

// Space the run refuses to consume. Folding bounds the footprint but cannot bound it below the
// merged graph, so a hard floor is still needed -- this is what actually keeps the run from
// filling somebody else's filesystem. 2% of the volume, clamped so it is neither trivial on a
// laptop nor wasteful on a 100 TB array.
inline size_t getDiskReserveBytes(const std::string & path) {
    const size_t MIN_RESERVE = 64ULL * 1024 * 1024;
    const size_t MAX_RESERVE = 1024ULL * 1024 * 1024;
    const size_t total = getTotalDiskBytes(path);
    if (total == 0) {
        return MIN_RESERVE;
    }
    return std::max(MIN_RESERVE, std::min(MAX_RESERVE, total / 50));
}

// Disk counterpart of getMemoryBudgetBytes. Note what this value is NOT: it is not a
// performance target. Anchoring it at a small fraction of free space would make a laptop fold
// constantly while leaving a 100 TB filesystem unprotected, and would make the same data behave
// differently depending on what else happens to sit on the disk. 80% is a safety ceiling: a run
// that fits under it never folds early and takes the same path it takes today.
inline size_t getTmpDiskBudgetBytes(const std::string & outDir, int maxTmpDiskMiB) {
    if (maxTmpDiskMiB > 0) {
        return static_cast<size_t>(maxTmpDiskMiB) * 1024 * 1024;
    }
    const size_t freeBytes = getFreeDiskBytes(outDir);
    if (freeBytes == 0) {
        return TMP_DISK_UNLIMITED;
    }
    return freeBytes / 100 * 80;
}

// The budget below divides by the peak this many rounds can hold, so it has to come first.
//
// Flush rounds that may be harvested but not yet written at the same time.
//
// A round's harvest empties the unit buffers, so the budget is honoured the moment it runs -- but
// the harvested data stays in memory until the write finishes. Capping how many rounds may sit in
// that state is what keeps the peak bounded: resident is at most (this + 1) rounds' worth.
//
// A constant, deliberately not derived from --threads. Deriving it would put the peak back on the
// worker count, which is the dependency the 2026-08-27 and 2026-08-28 work exists to remove: the
// same input would then need a bigger machine on a bigger node. Two lets one round be written
// while the next fills, which is the whole overlap there is to win; a third would buy nothing
// because the harvest is serialised anyway.
inline size_t concurrentFlushRounds() {
    return 2;
}

// Bytes a buffered pair costs at the peak of a round.
//
// 48 is the unordered_map node. The extra 12 is sizeof(Relation): saveSubGraphToFile copies a
// unit's entries into a vector and only then releases that unit's map, so while a worker is
// copying it holds both representations. Across a whole round the vectors replace the maps unit
// by unit rather than adding to them, which bounds the transient at map + vector rather than at
// twice either -- but the bound is 60, not 48.
//
// The [edges] flush writers line reports a larger "vector peak" than 12 B per pair. That figure
// is writers x largest unit x sizeof(Relation), an upper bound that assumes every worker holds a
// copy of the biggest unit; the true total across workers cannot exceed one round's pairs.
inline size_t peakBytesPerBufferedPair() {
    return 48 + sizeof(Relation);
}

// Pairs the emit stage may hold in memory across ALL threads before flushing.
//
// This used to divide the budget by the thread count and hand each thread its own quota. Total
// capacity was the same either way, but each thread flushed independently the moment its own map
// filled, so a quarter-sized map meant four times as many flush events -- and each event writes
// one file per merge unit. CAMI2 strain-madness went from 389 flushes at 16 threads to 2,354 at
// 64, and with 3*partitions units per flush the same 117.29 GB landed in 451,968 files instead of
// 18,672. The merge then folded 23 times instead of once and the run took 8h 08m instead of
// 3h 40m. Shared here so both numbers stop tracking --threads.
//
// The divisor is the whole peak, not one round's share. Until 2026-09-04 it was 48 B with no
// account of concurrentFlushRounds(), and a MAX_THRESHOLD of 200 M silently made up the
// difference: at --max-ram 128 the formula produced 1.72e9 and the clamp cut it to 2e8, so the
// clamp -- not the formula -- was what kept the peak at the 28.8 GB the comment above used to
// quote. Removing the clamp without fixing the divisor would have asked for
// 3 x 1.72e9 x 48 B = 247 GB against a declared 128 GiB.
//
// With the peak in the divisor the clamp is unnecessary and the declared limit means what it
// says: --max-ram 128 now yields about 4.6e8 pairs and a peak near 82 GB, 60% of the declaration.
// That is 2.3x the old budget, so a run needs 2.3x fewer flush rounds and folds fewer times --
// folding is triggered by the round count passing the merge fan-in, not by disk.
//
// This does not change what the run produces. The budget decides only where round boundaries
// fall; which pairs exist and what weight each carries are properties of the k-mers. Weights
// combine with addSat16, and min(65535, sum) over non-negative addends is independent of how the
// sum is split. Verified on the 120x fixture: 110 rounds against 112, and 80 against 92, produced
// byte-identical output.
//
// MIN_THRESHOLD stays. --max-ram 1 now yields about 3.6e6, still above it, but a smaller
// declaration or a machine with little free memory must not drive the budget to zero.
inline size_t getRelationBudget(int maxRamGiB) {
    const size_t availableBytes = getMemoryBudgetBytes(maxRamGiB);

    const double safetyFactor = 0.6;
    const size_t peakBytesPerPair = peakBytesPerBufferedPair() * (concurrentFlushRounds() + 1);

    const size_t budget = (size_t)(availableBytes * safetyFactor) / peakBytesPerPair;

    const size_t MIN_THRESHOLD = 1'000'000;
    return std::max(MIN_THRESHOLD, budget);
}

// What the budget above implies the emit stage will hold at its peak, for the log. Reported so
// the declared --max-ram and the number actually derived from it are visible side by side.
inline size_t projectedEmitPeakBytes(int maxRamGiB) {
    return getRelationBudget(maxRamGiB) * peakBytesPerBufferedPair()
           * (concurrentFlushRounds() + 1);
}

// Which edge set one shared k-mer produces.
//   clique: every pair, C(m,2) edges. An edge weight then counts the k-mers two reads share,
//           which is what --min-overlap-ratio x k-mers-per-read assumes.
//   star:   a hub and its m-1 spokes. Linear in m instead of quadratic, at the cost of that
//           weight meaning: two reads that share a k-mer get no edge from it unless one of them
//           is the hub.
enum EdgeMode {
    EDGE_MODE_CLIQUE = 0,
    EDGE_MODE_STAR = 1
};

// The hub for one k-mer's read list: the read holding the most k-mers, ties going to the smaller
// id. Reading the counts rather than the list order is what keeps the choice independent of how
// the query file happens to be sorted -- ids follow input order, so "first in the list" would not
// be.
//
// ids is sorted ascending, and the comparison is strictly greater, so the first read to reach the
// maximum keeps it: the tie rule needs no separate branch. ids must not be empty.
inline uint32_t pickStarHub(const std::vector<uint32_t> & ids,
                            const std::vector<uint32_t> & kmerCntPerRead) {
    uint32_t hub = ids[0];
    uint32_t best = (hub < kmerCntPerRead.size()) ? kmerCntPerRead[hub] : 0;
    for (size_t i = 1; i < ids.size(); ++i) {
        const uint32_t id = ids[i];
        const uint32_t cnt = (id < kmerCntPerRead.size()) ? kmerCntPerRead[id] : 0;
        if (cnt > best) {
            best = cnt;
            hub = id;
        }
    }
    return hub;
}

// Workers that write one flush round's unit files at the same time.
//
// The unit loop is the round's whole cost and it is embarrassingly parallel: a unit owns its map,
// its sort and its file, and two units never touch the same bytes. Splitting it is what takes a
// round off the critical path.
//
// Bounded by the thread count because that is what the user asked for -- a job given four cores by
// its scheduler gets no more CPU by spawning sixteen writers, only more contention -- and by the
// unit count because there is nothing for a further worker to do. No constant ceiling: the useful
// maximum is set by how lopsided a round is, and that follows --partitions (the cross route is
// sharded that many ways), so a literal here would be right at one partition count and wrong at
// every other.
//
// Deriving this from --threads does not reintroduce the dependency that 0b44e30d and e984f9f1
// removed. Those were about the grouping itself; this only decides who writes which file, and a
// file's name and contents are fixed by its unit. concurrentFlushRounds() is the opposite case and
// stays a constant: it bounds resident memory, which must not follow the core count.
inline size_t subGraphWriters(size_t unitCnt, int threads) {
    const size_t workers = static_cast<size_t>(threads < 1 ? 1 : threads);
    return std::min(unitCnt < 1 ? 1 : unitCnt, workers);
}

// Width of one id range in the relations_* routing scheme. The count is --partitions, not
// --threads: the two were the same number until 2026-08-26, which is what made the result depend
// on the thread count.
//
// Rounded up: read ids run 1..processedReadCnt, so ceil(N/P) keeps id / rangeSize <= P and the
// same-range routes inside partitionCnt+1..2*partitionCnt. floor(N/P) let the last ids reach
// range P+1 or beyond (N=11, P=4: ids 10 and 11 -> route 9 > 2P), past emitUnitCount.
inline size_t getRouteRangeSize(size_t processedReadCnt, int partitionCnt) {
    const size_t parts = static_cast<size_t>(partitionCnt < 1 ? 1 : partitionCnt);
    return (processedReadCnt + parts - 1) / parts;
}

// Which relations_{r}.bin an edge belongs to. Depends only on (id1, id2), which is what lets
// the merge be split by route: a pair can never appear under two different routes. Each route is
// then read by exactly one reader, so a per-route counting map sees every observation of a pair.
// Routes 0..partitionCnt-1 hold the same-residue edges and partitionCnt+1..2*partitionCnt the
// same-range ones; partitionCnt is the cross bucket that spans partitions.
inline size_t routeOf(uint32_t id1, uint32_t id2, int partitionCnt, size_t rangeSize) {
    const uint32_t parts = static_cast<uint32_t>(partitionCnt < 1 ? 1 : partitionCnt);
    if (id1 % parts == id2 % parts) {
        return static_cast<size_t>(id1 % parts);
    }
    // rangeSize is 0 only when there are no reads at all, in which case there are no edges
    // either; guard anyway so the helper can never divide by zero.
    if (rangeSize != 0 && id1 / rangeSize == id2 / rangeSize) {
        // Clamped so an id above processedReadCnt can never index past the last route.
        const size_t range = std::min(static_cast<size_t>(id1 / rangeSize), static_cast<size_t>(parts));
        return range + static_cast<size_t>(parts);
    }
    return static_cast<size_t>(parts);
}

// How many shards a route is split into for merging. Only the cross bucket needs it: the
// routing rule sends ~88% of all edges there (measured), so merging it as one unit caps the
// whole parallel merge at ~1.13x no matter how many threads are available.
inline size_t shardsForRoute(size_t route, int partitionCnt) {
    if (partitionCnt < 1) {
        return 1;
    }
    return (route == static_cast<size_t>(partitionCnt)) ? static_cast<size_t>(partitionCnt) : 1;
}

// Which shard of its route an edge belongs to. Keyed on id1 ALONE -- keying on anything
// derived from both ids would let one pair land in two shards and split its weight.
inline size_t shardOf(uint32_t id1, size_t shardCnt) {
    if (shardCnt <= 1) {
        return 0;
    }
    return static_cast<size_t>(id1 % static_cast<uint32_t>(shardCnt));
}

// Flat, gapless index over the (route, shard) units, keyed on --partitions. Distinct from the
// class's unitIndexOf/unitCount, which derive the same layout from --threads: the two agreed
// until --partitions was split off from --threads, and only this pair is correct once they
// differ. Gapless because every unit here owns a mutex and a map, and the strided form would
// allocate partitionCnt of each per single-shard route.
//
// Layout: routes below the cross bucket hold one shard each and index directly; the cross
// bucket's shards follow; routes above it are offset past them.
inline size_t emitUnitCount(int partitionCnt) {
    const size_t cross = static_cast<size_t>(partitionCnt < 1 ? 1 : partitionCnt);
    const size_t routeCnt = cross * 2 + 1;
    return (routeCnt - 1) + shardsForRoute(cross, partitionCnt);
}

inline size_t emitUnitIndex(size_t route, size_t shard, int partitionCnt) {
    const size_t cross = static_cast<size_t>(partitionCnt < 1 ? 1 : partitionCnt);
    if (route < cross) {
        return route;
    }
    if (route == cross) {
        return cross + shard;
    }
    return cross + shardsForRoute(cross, partitionCnt) + (route - cross - 1);
}

// One emit buffer per (unit, key shard), shared by every thread. Holding the map per unit rather
// than per thread is what takes the flush unit off --threads: a flush empties all of these at
// once, so what lands in a file is the whole budget divided by the unit count, whatever the
// thread count. Non-copyable and non-movable because of the mutex, hence held by unique_ptr in a
// vector.
struct EmitUnitBuffer {
    std::mutex lock;
    std::unordered_map<uint64_t, uint16_t> pairs;
};

// Lock domains a unit's buffer is split into. One buffer per unit put 128 threads on 48 mutexes
// and 74.74% of the takes had to wait for one (CAMI2 strain-madness, job 550906); striping the
// unit's keys over several mutexes multiplies the domains by this factor.
//
// Two domains per worker, so a worker usually finds a free one on its first try. A power of two
// because the shard of a key is picked inside the C(m,2) loop -- the hottest loop in the program
// -- where a 64-bit division would cost more than the hash insert it is there to spread. Capped
// at 16 because a thread's staging batch is split S ways to keep its footprint fixed, and past
// that the sub-batch is too small to amortise the lock it exists to amortise.
//
// --threads 1 gives 1, which is byte-for-byte the pre-shard code path.
inline size_t emitUnitShards(size_t unitCnt, int threads) {
    const size_t workers = static_cast<size_t>(threads < 1 ? 1 : threads);
    const size_t units = (unitCnt < 1) ? 1 : unitCnt;
    const size_t want = (workers * 2 + units - 1) / units;
    size_t shards = 1;
    while (shards < want && shards < 16) {
        shards <<= 1;
    }
    return shards;
}

// Which of a unit's shards a pair belongs to. A function of the key alone, which is the whole
// reason this split cannot change what a file holds: every occurrence of a pair meets in one
// shard, so a unit's shards hold disjoint key sets and their union is the map they replace --
// same keys, same weights, same count. Sharding by thread instead would scatter one pair over
// several maps, turning it into several records and making the budget count it several times.
//
// Multiplied before masking: the low bits of pairKey are id2 and the high bits id1, and within a
// unit id1 is already constrained by routeOf/shardOf, so masking the raw key would leave the
// spread to id2 alone. shardCnt is a power of two, so the mask replaces a division.
inline size_t emitShardOf(uint64_t pairKey, size_t shardCnt) {
    const uint64_t mixed = pairKey * 0x9E3779B97F4A7C15ULL;
    return static_cast<size_t>(mixed >> 32) & (shardCnt - 1);
}

// How many pair occurrences a thread stages for one unit before taking that unit's lock. The
// insertion site is the C(m,2) loop, the hottest in the program -- locking per occurrence would
// cost about as much as the hash insert it guards. Batching amortises it away at the price of
// `emitUnitCount * STAGE * 8` bytes per thread (1.5 MB at 48 units).
static const size_t EMIT_STAGE_BATCH = 4096;

// Soft cap on how many files this process may hold open.
inline size_t getOpenFileLimit() {
    struct rlimit lim;
    if (getrlimit(RLIMIT_NOFILE, &lim) == 0 && lim.rlim_cur != RLIM_INFINITY) {
        return static_cast<size_t>(lim.rlim_cur);
    }
    return 1024; // POSIX-typical default when the limit cannot be read
}

// How many subgraph files one merge pass may open at once. Bounded by the fd limit -- with a
// reserve for the relations_* writers, stdio, and the k-mer readers -- and by a cap that keeps
// each stream's read buffer large enough for sequential I/O. Anything beyond this is handled
// by merging in rounds rather than by opening more files.
// Descriptors one merge unit holds at its peak: `fanIn` readers plus the single file it is
// writing (its output, or a fold intermediate).
// Descriptors this process currently holds, or 0 where it cannot be counted. Predictions of
// fd use are only as good as the model behind them; this is the measurement to check them
// against, and the number that decides whether a run survives on a host with a tighter limit.
inline size_t openFdCount() {
#if defined(__linux__)
    DIR * const dir = opendir("/proc/self/fd");
    if (dir == nullptr) {
        return 0;
    }
    size_t count = 0;
    while (readdir(dir) != nullptr) {
        ++count;
    }
    closedir(dir);
    // ".", ".." and the descriptor opendir itself is holding.
    return (count > 3) ? (count - 3) : 0;
#else
    return 0;   // no portable equivalent; the prediction stands alone there
#endif
}

inline size_t mergeFdsPerUnit(size_t fanIn) {
    return fanIn + 1;
}

// Units maybeFoldEmitted may fold at the same time.
//
// Folding was serial, one unit after another, and the comment justifying it said only that this
// keeps one unit's streams open at a time. That is true but it is not a reason: units are
// independent -- a unit's inputs are its own subGraph_* files and its output name carries its
// route and shard, so no two units touch the same bytes.
//
// The cost of the serial loop is not small. CAMI2 plant-associated folded once and that fold took
// 34,003 s of a 43,463 s emit -- 78% of the stage, 69% of the whole run -- with 48 units averaging
// 708 s each on one thread while 127 cores idled.
//
// Derived from what is actually free, not from a constant. The scan holds threads x splits x 2
// descriptors for the whole emit, so the headroom at fold time is a property of the run: on
// plant-associated the log recorded 38,712 open against a 65,536 limit, leaving room for 52 units
// at 513 descriptors each. On an input with more splits there is less, and on one with a smaller
// limit there may be room for one -- which is the serial loop, unchanged.
inline size_t foldWorkers(size_t unitCnt, size_t fanIn) {
    const size_t units = (unitCnt < 1) ? 1 : unitCnt;
    const size_t limit = static_cast<size_t>(getOpenFileLimit() * 0.75);
    const size_t held = openFdCount();
    const size_t spare = (limit > held) ? (limit - held) : 0;
    const size_t perUnit = mergeFdsPerUnit(fanIn);
    const size_t fits = (perUnit > 0) ? spare / perUnit : 1;
    return std::max<size_t>(1, std::min(units, fits));
}

// Descriptors scanKmerRuns needs, all at once.
//
// Each thread opens both files of every split and merges them in k-mer order, so all of them have
// to be open together -- it is a k-way merge, not a sequential pass. The demand is therefore
// threads x splits x 2, and splits grow with the read count.
//
// Nothing used to check it. CAMI3 human-gut toy at 332,762,102 reads produced 183 splits, which
// needs 128 x 183 x 2 = 46,848 descriptors on top of everything else the process holds; the opens
// past the limit failed one by one and ReadBuffer exited from inside the OpenMP region, so the run
// died with a segmentation fault and a handful of scattered "Error opening file" lines rather than
// a diagnosis (job 551117).
//
// Throttling the threads in that region is NOT a fix: the thread index is a partition index. The
// k-mer files are named kmer_delta_<split>_<thread> and written by a region sized from --threads,
// and both the scan and the cleanup enumerate 0..threads-1. Running the scan with fewer workers
// would leave whole partitions unread. Making it safe means lowering --threads for the entire run,
// which is a decision for the caller, so this reports what is needed and lets the caller act.
inline size_t scanKmerRunsFdDemand(size_t numOfSplits, int threads) {
    const size_t workers = static_cast<size_t>(threads < 1 ? 1 : threads);
    return workers * numOfSplits * 2;
}

// Largest fan-in whose peak descriptor use fits the process limit.
//
// `concurrentMergers` units run at once, each holding mergeFdsPerUnit(fanIn). An earlier
// version reserved a flat slack and forgot the per-unit writer, which put the peak one
// descriptor over an 8192 limit at 16 threads and killed the merge with a bare
// "Error opening file". Only a fraction of the limit is spent so that descriptors held
// elsewhere -- and the writer this arithmetic now counts -- cannot tip it over.
// `concurrentMergers` is how many merges will hold streams open AT THE SAME TIME -- not how
// many threads the run has. The two were passed the same value at the fold site, which made the
// fold's fan-in shrink as --threads grew even though the fold merges one unit at a time. On a
// host with the common `ulimit -n 1024` that took the fan-in from 512 at one thread to 11 at
// sixty-four, and a smaller fan-in means more fold rounds over the same data.
// The old first parameter was already dead (`(void) numThreads`); it is gone so the mistake
// cannot be repeated.
inline size_t getMergeFanIn(size_t concurrentMergers) {
    const size_t FAN_IN_CAP = 512;
    const double utilisation = 0.75;
    const size_t mergers = std::max<size_t>(1, concurrentMergers);
    const size_t budget = static_cast<size_t>(getOpenFileLimit() * utilisation);
    // Solve mergers * (fanIn + 1) <= budget.
    const size_t perUnit = (budget > mergers) ? (budget / mergers - 1) : 0;
    return std::max<size_t>(2, std::min(FAN_IN_CAP, perUnit));
}

// m-distribution histogram: reads-per-k-mer bucketed by floor(log2(m)). m can reach the
// read count (tens of millions), so a flat array is not affordable per thread; log2 buckets
// give cap candidates at factor-2 spacing, which is enough to pick a disk budget.
static const size_t M_HIST_BUCKETS = 64;

static inline size_t mHistBucket(size_t m) {
    size_t bucket = 0;
    while ((m >> bucket) > 1 && bucket + 1 < M_HIST_BUCKETS) {
        ++bucket;
    }
    return bucket; // floor(log2(m)) for m >= 1; 0 for m == 0
}

// Smallest log2-bucket upper bound (3, 7, 15, 31, ...) covering `quantile` of the k-mers with
// m >= 2, i.e. the cap that keeps that share of distinct k-mers and drops the tail above it.
//
// The normaliser is the k-mer count, not the pair count. m does not scale with the size of the
// dataset -- it scales with per-genome coverage, so a fraction of the read count cannot serve
// (which is why --max-kmer-freq-ratio was removed). Weighting by pairs would not work either:
// on CAMI2 strain-madness a single bucket holds 52% of all pair occurrences, so a pair-weighted
// quantile cuts into the body. The k-mer count is what locates the body/tail boundary.
//
// Returns 0 for "no cap": the quantile is outside (0, 1), the histogram is empty, or the
// distribution never reaches the quantile. Pure function of the histogram, so two runs over the
// same data pick the same cap regardless of thread count.
static inline size_t capFromQuantile(const std::vector<uint64_t> & mHistKmers, float quantile) {
    if (!(quantile > 0.0f) || !(quantile < 1.0f)) {
        return 0;
    }
    uint64_t total = 0;
    for (size_t b = 0; b < mHistKmers.size(); ++b) {
        total += mHistKmers[b];
    }
    if (total == 0) {
        return 0;
    }
    const double target = static_cast<double>(total) * static_cast<double>(quantile);
    uint64_t cumulative = 0;
    for (size_t b = 0; b < mHistKmers.size(); ++b) {
        cumulative += mHistKmers[b];
        if (static_cast<double>(cumulative) >= target) {
            // Bucket b holds m in [2^b, 2^(b+1) - 1]; the cap is that range's upper bound.
            // Shifting by 63 or more is undefined, and a cap that large caps nothing anyway.
            if (b + 1 >= 63) {
                return 0;
            }
            return (static_cast<size_t>(1) << (b + 1)) - 1;
        }
    }
    return 0;
}

// Floor for mergeGraph's per-stream read buffer. Reaching it means numOfGraph is so large
// that even the memory budget cannot give each stream a useful buffer; mergeGraph warns,
// and the real fix is multi-round merging (not implemented).
static const size_t MERGE_BUFFER_MIN_ELEMS = 4'096;

// Per-stream ReadBuffer size (in Relation elements) for the k-way subgraph merge.
// The merge opens one buffer per subGraph_* file, so a fixed per-stream size makes total
// memory scale with numOfGraph: the historical 1M elements is 12 MB per stream, i.e. 12 GB
// at 1000 streams. Budget the total against free memory and split it across streams.
inline size_t getMergeBufferElems(size_t numOfGraph, int maxRamGiB, size_t concurrentMergers) {
    const size_t MAX_ELEMS = 1'048'576; // 12 MB per stream; the historical fixed size
    if (numOfGraph == 0) {
        return MAX_ELEMS;
    }

    const double safetyFactor = 0.5;
    const size_t mergers = std::max<size_t>(1, concurrentMergers);
    const size_t budget = (size_t)(getMemoryBudgetBytes(maxRamGiB) * safetyFactor) / mergers;

    // Disk bytes and memory bytes are not the same thing here. Each compressed stream also
    // holds a fixed decompressor input window; charge that first, then split what is left
    // across the record buffers. Conflating the two is how the merge blew past its fd budget
    // once before (see research.md, 2026-08-06 (2)).
    const size_t fixedCost = relationStreamFixedBytes() * numOfGraph;
    const size_t forRecords = (budget > fixedCost) ? (budget - fixedCost) : 0;
    const size_t perStream = forRecords / (numOfGraph * sizeof(Relation));

    return std::max(MERGE_BUFFER_MIN_ELEMS, std::min(perStream, MAX_ELEMS));
}

class GroupGenerator {
protected:
    const LocalParameters & par;
    string commonKmerDB;
    string outDir;
    size_t matchPerKmer;
    int kmerFormat;
    int printLog;
    int commonKmerSpan;   // --common-kmer-span; see the note above the class
    // The one number the routing scheme is built on. Resolved once, in the constructor, from
    // --partitions falling back to --threads, and read everywhere the old code read par.threads
    // for a routing decision. Keeping the two apart is the point: --threads may then change
    // without changing the result, which it used to do through the support pass's cap.
    int partitionCnt;
    // Which edge set a shared k-mer produces. Not a CLI flag: the easy-grouping command sets it,
    // and exposing it as well would make `grouping --edge-mode 1` and `easy-grouping` two names
    // for one thing, with two families of run tags to keep straight.
    int edgeMode;

    // Agents    
    GeneticCode * geneticCode = nullptr;
    QueryIndexer *queryIndexer = nullptr;
    KmerExtractor *kmerExtractor = nullptr;

    size_t numOfSplits = 0;
    size_t numOfGraph = 0;
    // Bytes the subGraph_* files actually occupy. saveSubGraphToFile runs on every emit
    // thread, hence the atomic. Compared against emittedEdgeCnt * sizeof(Relation) it gives
    // the compression ratio the run achieved.
    std::atomic<uint64_t> subGraphBytesOnDisk{0};

    // Compressed size of every subGraph_* file written, for the size distribution the [edges]
    // summary reports. The total alone cannot distinguish "the same data in half as many files"
    // from "half the data" -- and that distinction is the whole question about --threads, since
    // raising it splits one flush budget across more maps and so writes more, smaller files.
    // Appended once per flush, not once per file: saveSubGraphToFile collects a flush's unit
    // sizes locally and takes the lock once. One entry per file at 8 B is a few MB at the
    // largest file counts seen (472,320 files -> 3.8 MB).
    std::mutex subGraphSizeMutex;
    std::vector<uint64_t> subGraphFileBytes;

    // --- Incremental fold state -------------------------------------------------------------
    // The merge used to start only after every flush had been written, so peak disk was the
    // whole emit. Folding a completed prefix while emission continues bounds it instead.
    //
    // A flush index is handed out by counter.fetch_add BEFORE its files exist, so "written" is
    // tracked separately: emitDone[i] is set once all of flush i's unit files are closed, and
    // emitWatermark walks the contiguous prefix of those. Only that prefix may be folded --
    // folding across a hole would list a file that has not been written yet.
    std::mutex foldMutex;                                // guards the five members below
    std::vector<char> emitDone;
    size_t emitWatermark = 0;                            // flushes [0, emitWatermark) are complete
    size_t foldedUpTo = 0;                               // flushes [0, foldedUpTo) folded away
    size_t foldRounds = 0;
    std::vector<std::vector<std::string>> foldedOutputs; // [unit] -> folded file names

    // What folding costs. Nothing reported it before: the [disk] fold line carries file counts,
    // bytes and descriptors but no duration, so the only evidence was the flush write worst round
    // -- and that is the time of a round that happened to overlap a fold, not the fold itself.
    // A fold rewrites every unit's accumulated output, so on a large dataset it is hours; without
    // this the effect of changing the flush budget or parallelising the unit loop cannot be judged.
    std::atomic<uint64_t> foldMillisTotal{0};
    std::atomic<uint64_t> foldMillisMax{0};      // worst single fold
    std::atomic<uint64_t> foldUnitMillisMax{0};  // worst single unit within any fold
    std::atomic<uint64_t> foldUnitsDone{0};      // unit-merges performed, for the per-unit mean

    // Where a flush round's time goes, split four ways, summed over every round in microseconds.
    // The round total is already reported; this says whether it is the sort, the map-to-vector
    // build, zstd or releasing the maps that owns it, which is what decides whether parallelising
    // the unit loop is worth doing and what its ceiling is. Rounds overlap, so these sums can
    // exceed the wall clock. `free` was measured last and was the largest single item on the
    // calling thread before it moved into the workers.
    std::atomic<uint64_t> emitBuildMicros{0};
    std::atomic<uint64_t> emitSortMicros{0};
    std::atomic<uint64_t> emitWriteMicros{0};
    std::atomic<uint64_t> emitFreeMicros{0};

    // How lopsided a round is. Splitting the unit loop across workers cannot go faster than its
    // largest unit, so the largest unit's share of a round is the speed-up ceiling: a round whose
    // biggest unit holds 1/x of the records can at best be x times faster. Summed over rounds so
    // the ratio of the two gives the average share; the worst single round is tracked separately
    // in hundred-thousandths -- coarser scaling truncates the worst round below the average,
    // which reads as a contradiction.
    std::atomic<uint64_t> emitMaxUnitRecords{0};
    std::atomic<uint64_t> emitRoundRecords{0};
    std::atomic<uint64_t> emitWorstShareScaled{0};

    // What parallelising the unit loop adds to the peak: each worker holds one unit's records as a
    // vector while it sorts and writes them, so the bound is writers x the largest unit. Reported
    // next to the round's maps, which hold the same pairs at 48 B each against this 12.
    std::atomic<uint64_t> emitVectorPeakBytes{0};
    std::atomic<size_t> emitWriterCnt{0};

    // Live footprint of subGraph_* and its high-water mark. Signed because a fold subtracts its
    // inputs before adding its output.
    std::atomic<int64_t> subGraphLiveBytes{0};
    std::atomic<uint64_t> subGraphPeakBytes{0};

    // Emit progress, for projecting the total write volume. The [mhist] table cannot serve here:
    // it is only assembled after the parallel region joins. These two are updated as k-mers are
    // consumed, so a projection is available while the emit is still running.
    std::atomic<uint64_t> kmerValuesTotal{0};
    std::atomic<uint64_t> kmerValuesDone{0};
    std::atomic<int> volumeWarned{0};   // the projection warning is printed at most once

    std::vector<uint64_t> kmerBoundaries;
    bool boundariesInitialized = false;
    bool useOnlyTrueRelations = false; // for debug

    // Filter-surviving k-mers, accumulated across splits. Divided by the read count it gives
    // k-mers per read, which is the scale the Phase 1 edge-weight threshold is expressed in:
    // an edge weight is the number of k-mers two reads share, so weight/kmersPerRead is the
    // fraction of a read the two overlap by.
    size_t totalFilteredKmers = 0;

    // The same count broken down per read, indexed by global read id (1-based, so element 0 is
    // unused). Counted over exactly the k-mers that survive the common-k-mer filter, which is the
    // set the graph sees and the set totalFilteredKmers sums -- Sum(kmerCntPerRead) must equal it.
    // Only the star edge mode reads this: it picks the read with the most k-mers as the hub, so
    // the choice is a property of the data rather than of the input file's order. Four bytes per
    // read, held for the whole run.
    std::vector<uint32_t> kmerCntPerRead;

    // Merged edge weights, in the same log2 buckets the m histogram uses. Filled by mergeGraph and
    // read once the core threshold is known, to say where that threshold sits in the distribution
    // it is cutting. The clique and star modes produce very different weights for the same data --
    // star only records an edge when one of the two reads was the hub -- so this is what makes the
    // two comparable at a glance.
    std::vector<uint64_t> edgeWeightHist;

public:
    GroupGenerator(LocalParameters & par);

    ~GroupGenerator();

    void startGroupGeneration(const LocalParameters & par);
    
    // processedReadCnt is the number of reads finished before this split, i.e. the offset that
    // turns a split-local sequence id into the global one writeKmers later stamps into the k-mer
    // files. Needed here because the per-read counts are accumulated globally.
    void filterCommonKmers(Buffer<Kmer>& queryKmerBuffer,
                           Buffer<std::pair<uint32_t, uint32_t>> & matchBuffer,
                           const string & db="",
                           size_t processedReadCnt=0);

    void writeKmers(Buffer<Kmer>& queryKmerBuffer,
                    size_t processedReadCnt);

    std::vector<std::pair<size_t, size_t>> getKmerRanges(const Buffer<Kmer>& kmerBuffer,
                                                         size_t offset);

    void makeSubGraph(size_t processedReadCnt);

    // Walks one thread's kmer_delta_*/kmer_info_* files in k-mer order and hands each distinct
    // k-mer's sorted, deduplicated read ids to `onKmer`. Template rather than std::function: the
    // callback fires once per distinct k-mer (452 M times on CAMI2 strain-madness) and the scan
    // steps through 4.9e9 k-mer instances, so an indirect call per k-mer is not free.
    //
    // `valuesScanned` is incremented, never assigned. makeSubGraph publishes it into
    // kmerValuesDone at every flush and then zeroes it, and the disk projection depends on that
    // cadence -- the callback runs inside this loop, so it can reset the counter itself.
    //
    // kmerValuesTotal is accumulated here, once per reader. A caller that runs the scan twice
    // has to reset it in between, or the projection's denominator doubles.
    template <typename OnKmer>
    void scanKmerRuns(int threadIdx, size_t processedReadCnt, uint64_t & valuesScanned,
                      OnKmer && onKmer) {
        std::vector<std::unique_ptr<DeltaIdxReader>> readers;
        readers.reserve(numOfSplits);
        std::vector<Kmer> currentKmers;
        currentKmers.reserve(numOfSplits);
        for (size_t i = 0; i < numOfSplits; ++i) {
            const std::string diffFile = outDir + "/kmer_delta_" + std::to_string(i)
                                       + "_" + std::to_string(threadIdx);
            const std::string infoFile = outDir + "/kmer_info_" + std::to_string(i)
                                       + "_" + std::to_string(threadIdx);
            readers.push_back(std::unique_ptr<DeltaIdxReader>(
                new DeltaIdxReader(diffFile, infoFile, 1024 * 1024, 1024 * 1024)));
            currentKmers.push_back(readers.back()->next());
            // Denominator of the emit progress used to project the write volume.
            kmerValuesTotal.fetch_add(readers.back()->getTotalValueNum(), std::memory_order_relaxed);
        }

        std::vector<uint32_t> queryIds;
        queryIds.reserve(1024);

        while (true) {
            // Find the smallest k-mer
            uint64_t minKmer = UINT64_MAX;
            for (size_t file = 0; file < numOfSplits; ++file) {
                if (!readers[file]->isCompleted()) {
                    minKmer = std::min(minKmer, currentKmers[file].value);
                }
            }
            if (minKmer == UINT64_MAX) break;

            queryIds.clear();
            for (size_t file = 0; file < numOfSplits; ++file) {
                while (currentKmers[file].value == minKmer) {
                    const uint32_t seqId = currentKmers[file].tInfo.taxId; // query ID is in taxId
                    if (seqId != UINT32_MAX && seqId <= processedReadCnt) {
                        queryIds.emplace_back(seqId);
                    }
                    ++valuesScanned;
                    currentKmers[file] = readers[file]->next();
                    if (readers[file]->isCompleted()) {
                        currentKmers[file].value = UINT64_MAX;
                        break;
                    }
                }
            }

            std::sort(queryIds.begin(), queryIds.end());
            queryIds.erase(std::unique(queryIds.begin(), queryIds.end()), queryIds.end());
            onKmer(minKmer, queryIds);
        }
    }

    // Counting pre-pass: the same scan with no edge emission, so the reads-per-k-mer distribution
    // is known before makeSubGraph has to choose a cap. Only k-mers with m >= 2 are counted,
    // which is the domain the [mhist] table reports.
    void scanMHistogram(size_t processedReadCnt,
                        std::vector<uint64_t> & mHistKmers,
                        std::vector<uint64_t> & mHistPairs);


    // Writes one file per unit for one flush index: subGraph_{counter}_{route}_{shard}.
    // Takes the units already partitioned, because the emit buffers are now kept per unit --
    // the routing that used to happen here now happens at insertion time. `unitKey[u]` is the
    // (route, shard) of unit u, so the caller's flat index survives into the file name.
    // Every unit gets a file even when empty: the merge lists its inputs by name without
    // checking existence, and a missing file would be read as a short stream.
    // `units[u]` is the unit's key shards, whose key sets are disjoint (see emitShardOf), so the
    // file is built by concatenating them and sorting once -- no re-summing of weights.
    // Takes `units` by mutable reference and empties each unit as soon as its records are copied
    // out: a round holds up to RELATION_BUDGET individually allocated nodes, and freeing them on
    // the calling thread afterwards cost 3.9x the parallel write it followed while still holding
    // the round's in-flight slot. Emptying here spreads that over the same workers and drops the
    // unit's resident bytes from 60 per pair to the vector's 12 before the sort begins.
    void saveSubGraphToFile(
            std::vector<std::vector<std::unordered_map<uint64_t, uint16_t>>> & units,
            const std::vector<std::pair<size_t, size_t>> & unitKey,
            const size_t counter_now);

    // No processedReadCnt: routing now happens in saveSubGraphToFile, which is where the
    // read count is needed. This function only folds and merges what is already partitioned.
    void mergeGraph();

    // Merge one batch of sorted subgraph files into a single sorted file, summing the weights
    // of duplicate (id1,id2) pairs. Returns the number of records written. Used to fold many
    // subgraphs down to a mergeable fan-in; it does no partitioning or histogramming, which
    // only the final pass performs.
    size_t mergeSubGraphBatch(const std::vector<std::string>& inputs,
                              const std::string& output,
                              size_t bufElems);

    // Reduce `files` to at most maxFanIn entries by merging in rounds, deleting each batch's
    // inputs as soon as it is folded so peak disk grows by one batch rather than a full copy.
    // `route`/`shard` only name the intermediates: units fold concurrently, so a shared name
    // would let two of them write the same file and silently mix their edges.
    void reduceSubGraphFanIn(std::vector<std::string>& files, size_t maxFanIn,
                             size_t route, size_t shard);

    // Merge every subgraph belonging to one (route, shard) unit into `outPath`. Purely
    // sequential -- mergeGraph runs these concurrently, one unit per thread.
    // `files` comes from unitInputFiles(): incremental folding replaces a prefix of the flushes
    // with folded files, so the unit's inputs can no longer be derived from numOfGraph alone.
    void mergeUnit(size_t route, size_t shard, const std::string& outPath,
                   const std::vector<std::string>& files,
                   size_t bufElems, size_t maxFanIn,
                   size_t& mergedOut, size_t& ceilingOut,
                   std::vector<uint64_t>& weightHistOut);

    // --- Incremental fold ---------------------------------------------------------------
    // Unit numbering lives in emitUnitCount/emitUnitIndex, above. There used to be a second
    // copy here that derived the same layout from --threads; the two agreed only while
    // --partitions followed --threads, and once they could differ the fold indexed
    // foldedOutputs with one and sized it with the other. With --partitions above --threads
    // that wrote past the end of the vector, and unitInputFiles -- which is bounds-checked --
    // silently read another unit's folded output instead. Removed rather than fixed in place,
    // so there is one definition to be wrong about.
    std::string subGraphName(size_t flushIdx, size_t route, size_t shard) const {
        return outDir + "/subGraph_" + std::to_string(flushIdx) + "_" + std::to_string(route)
             + "_" + std::to_string(shard);
    }

    // Live-footprint bookkeeping. `delta` is positive for a file written, negative for one
    // removed; the high-water mark is what the [disk] summary reports.
    void noteSubGraphBytes(int64_t delta);

    // Called once every unit file of `flushIdx` is closed. Advances the contiguous watermark.
    void markEmitComplete(size_t flushIdx);

    // Two guards, checked at every flush.
    //  - projection: extrapolate the write volume from emit progress and warn once if it will
    //    not fit. Only a warning: folding routinely keeps the footprint far below the total,
    //    so a projected overflow does not mean the run must fail.
    //  - headroom: if the filesystem's actual free space falls to the reserve, stop. This one
    //    depends on no estimate, which is what makes it the real guarantee.
    void checkDiskHeadroom();

    // Fold the completed-but-unfolded prefix when the live footprint reaches `tmpDiskBudget`,
    // or when that prefix alone would exceed the merge's fan-in. Returns without doing anything
    // if another thread is already folding -- the caller just keeps emitting.
    void maybeFoldEmitted(size_t tmpDiskBudget, size_t maxFanIn);

    // Every file the merge must read for one unit: folded outputs first, then the flushes that
    // were never folded.
    std::vector<std::string> unitInputFiles(size_t route, size_t shard) const;

    // How many streams one unit contributes to the merge. Every unit folds on the same rounds
    // and keeps the same unfolded tail, so one unit's count answers for all of them. This --
    // not numOfGraph -- is what the fan-in and buffer budgets have to be sized against once
    // folding has replaced part of the flush range.
    size_t unitStreamCount() const {
        const size_t folded = foldedOutputs.empty() ? 0 : foldedOutputs[0].size();
        return folded + (numOfGraph - foldedUpTo);
    }

    void mergeGraph_one(size_t processedReadCnt);
    
    void makeGroups(int groupKmerThr,
                    size_t processedReadCnt,
                    unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
                    vector<uint32_t> &queryGroupInfo);
    // Phase 2: merge Phase-1 units that are joined by several independent weak edges.
    // A single weak edge (few shared k-mers) is indistinguishable from chance sharing across
    // repeats or conserved regions; several between the same pair of units are not. Counting
    // support at the unit level -- rather than testing node-level triangles -- is what makes
    // this affordable: two singletons can only ever have one edge between them, so every pair
    // worth counting involves at least one multi-read core group.
    //
    // `supportRatio` > 0 requires the support to scale with the smaller unit's read count, and
    // switches the count from weak edges to the distinct reads on the smaller side that carry
    // one. Both changes matter: chance links grow with |u| * |v|, so a fixed count is met
    // automatically on high-coverage data, and one repeat-bearing read linked to many reads of
    // the other unit produces many edges but only one distinct read. 0 keeps the plain edge
    // count with the floor alone, which is the pre-ratio behaviour.
    // `maxUnitReads` > 0 restricts the merge to pairs whose two units are both that small, and the
    // component size is re-checked as merges land, so nothing leaves this phase holding
    // 2 * maxUnitReads reads or more. Purity is lost per read: one wrong join between large units
    // costs more than many wrong joins between small ones. 0 leaves the merge unbounded.
    void mergeBySupport(int coreThr,
                        int weakThr,
                        size_t processedReadCnt,
                        std::unordered_map<uint32_t, std::unordered_set<uint32_t>> &groupInfo,
                        std::vector<uint32_t> &queryGroupInfo);


    void makeGroupsPhase3(int groupKmerThr,
                          size_t processedReadCnt,
                          const std::vector<bool>& isSingleton,
                          unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
                          vector<uint32_t>& queryGroupInfo);

    void saveGroupsToFile(const unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
                          const vector<uint32_t>& queryGroupInfo);
    
};


// ---------------------------------------------------------------------------
// Grouping driver, shared by the grouping commands and `em --read-groupings 1`
// (implemented in workflow/groupGeneration.cpp).
// ---------------------------------------------------------------------------

// Defaults of the `grouping` command. Call before parseParameters so that flags win.
void setGroupGenerationDefaults(LocalParameters & par);

// Derives the k-mer format from --syncmer and validates the grouping thresholds. When the common
// k-mer DB records its k-mer settings (kmer_params), they are either adopted into par
// (adoptDbKmerSettings, for callers without --syncmer/--smer-len flags such as em) or checked
// against par, failing on a mismatch. Returns false with a human-readable message in `error`.
bool prepareGroupGeneration(LocalParameters & par, const std::string & commonKmerDb, std::string & error,
                            bool adoptDbKmerSettings = false);

// Groups the reads of `queryFiles` (one file, or two for paired-end; par.seqMode must agree)
// against `commonKmerDb`, writing groups/groupMap into `outDir`. GroupGenerator and
// QueryIndexer read their inputs from par.filenames in the grouping command's positional
// layout, so that layout is presented to them and the caller's filenames are restored.
bool runReadGrouping(LocalParameters & par,
                     const std::vector<std::string> & queryFiles,
                     const std::string & commonKmerDb,
                     const std::string & outDir);

#endif // GROUP_GENERATOR_H