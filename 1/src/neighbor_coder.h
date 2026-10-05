#pragma once

#include "range_coder.h"

#include <cstdint>
#include <vector>

namespace NGraphCompressor {
    // Number of degree buckets and of distinct "already coded edges" counters
    // that define the category (cell) of a candidate neighbor.
    inline constexpr int DEGREE_BUCKET_COUNT = 64;
    inline constexpr int CODED_LEVEL_COUNT = 12;
    inline constexpr int CELL_COUNT = DEGREE_BUCKET_COUNT * CODED_LEVEL_COUNT;

    int GetDegreeBucket(uint32_t degree);
    int GetCell(uint32_t degree, uint32_t codedEdgeCount);

    // Codes the graph vertex by vertex. Vertices are identified by positions
    // 0..n-1. For each vertex the coder transmits the sorted set of its
    // neighbors with larger positions; edges to smaller positions are already
    // known.
    //
    // A candidate is drawn with probability proportional to the weight of its
    // cell (degree bucket, number of its edges coded so far). The sorted set
    // of k neighbors is coded as successive order statistics of k independent
    // draws.
    class TNeighborCoder {
    private:
        // Fenwick tree node: total weight and number of live candidates.
        struct TNode {
            uint64_t Weight = 0;
            int64_t Live = 0;
        };

    public:
        TNeighborCoder(const std::vector<uint32_t>& degrees, const std::vector<uint64_t>& cellWeights);

        // Starts the vertex (vertices must go in increasing order) and returns
        // the number of its neighbors that have not been coded yet.
        uint32_t BeginVertex(int vertex);
        void EncodeNeighbor(int neighbor, TRangeEncoder* encoder);
        // Throws TGraphCompressorException if the stream is inconsistent.
        int DecodeNeighbor(TRangeDecoder* decoder);
        void EndVertex();

    private:
        uint64_t GetCurrentWeight(int vertex) const;
        void SetWeight(int vertex, uint64_t weight);
        TNode GetPrefix(int end) const;
        // Prepares the distribution of the next neighbor after `Previous_`.
        void PrepareStep();
        // Cumulative frequency of all candidates whose weight sum (relative to
        // the start of the current range) is `weight` and count is `live`.
        uint64_t GetBoundary(int64_t weight, int64_t live) const;
        void Accept(int neighbor);

    private:
        int VertexCount_;
        int TopBit_;
        std::vector<uint32_t> Degrees_;
        std::vector<uint32_t> CodedEdgeCounts_;
        std::vector<uint64_t> Weights_;
        std::vector<uint64_t> CellWeights_;
        std::vector<TNode> Tree_;
        uint64_t TotalWeight_;
        int64_t TotalLive_;

        // State of the vertex being coded.
        int Previous_;
        uint32_t Remaining_;
        TNode Base_;
        int64_t RangeWeight_;
        uint64_t FreeTotal_;
        std::vector<int> Chosen_;
    };
} // namespace NGraphCompressor
