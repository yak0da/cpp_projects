#include "neighbor_coder.h"

#include "exception.h"

#include <algorithm>
#include <cmath>

namespace NGraphCompressor {
    namespace {
        constexpr uint32_t EXACT_DEGREE_BUCKET_COUNT = 40;
        constexpr uint64_t STRUCTURE_TOTAL = MAX_CODER_TOTAL;
        // Every live candidate owns this many units of the coder total, so its
        // interval stays non-empty despite rounding of the continuous part.
        constexpr uint64_t UNITS_PER_CANDIDATE = 2;

        [[noreturn]] void ThrowCorrupted() {
            throw TGraphCompressorException("corrupted input");
        }
    } // anonymous namespace

    int GetDegreeBucket(uint32_t degree) {
        if (degree < EXACT_DEGREE_BUCKET_COUNT) {
            return static_cast<int>(degree);
        }
        int bucket = EXACT_DEGREE_BUCKET_COUNT;
        uint64_t limit = 48;
        while (degree >= limit && bucket < DEGREE_BUCKET_COUNT - 1) {
            ++bucket;
            limit = limit * 6 / 5;
        }
        return bucket;
    }

    int GetCell(uint32_t degree, uint32_t codedEdgeCount) {
        const uint32_t level = std::min<uint32_t>(codedEdgeCount, CODED_LEVEL_COUNT - 1);
        return GetDegreeBucket(degree) * CODED_LEVEL_COUNT + static_cast<int>(level);
    }

    TNeighborCoder::TNeighborCoder(const std::vector<uint32_t>& degrees, const std::vector<uint64_t>& cellWeights)
        : VertexCount_(static_cast<int>(degrees.size()))
        , TopBit_(1)
        , Degrees_(degrees)
        , CodedEdgeCounts_(degrees.size(), 0)
        , Weights_(degrees.size(), 0)
        , CellWeights_(cellWeights)
        , Tree_(degrees.size() + 1)
        , TotalWeight_(0)
        , TotalLive_(0)
        , Previous_(-1)
        , Remaining_(0)
        , Base_()
        , RangeWeight_(0)
        , FreeTotal_(0)
        , Chosen_()
    {
        while (TopBit_ * 2 <= VertexCount_) {
            TopBit_ *= 2;
        }
        for (int vertex = 0; vertex < VertexCount_; ++vertex) {
            Weights_[vertex] = GetCurrentWeight(vertex);
            const int64_t live = Weights_[vertex] > 0 ? 1 : 0;
            TotalWeight_ += Weights_[vertex];
            TotalLive_ += live;
            const int index = vertex + 1;
            Tree_[index].Weight += Weights_[vertex];
            Tree_[index].Live += live;
            const int parent = index + (index & -index);
            if (parent <= VertexCount_) {
                Tree_[parent].Weight += Tree_[index].Weight;
                Tree_[parent].Live += Tree_[index].Live;
            }
        }
    }

    uint32_t TNeighborCoder::BeginVertex(int vertex) {
        SetWeight(vertex, 0);
        Previous_ = vertex;
        Remaining_ = Degrees_[vertex] - CodedEdgeCounts_[vertex];
        return Remaining_;
    }

    void TNeighborCoder::EncodeNeighbor(int neighbor, TRangeEncoder* encoder) {
        PrepareStep();
        const TNode before = GetPrefix(neighbor);
        const int64_t weight = static_cast<int64_t>(before.Weight - Base_.Weight);
        const int64_t live = before.Live - Base_.Live;
        const uint64_t low = GetBoundary(weight, live);
        const uint64_t high = GetBoundary(weight + static_cast<int64_t>(Weights_[neighbor]), live + 1);
        encoder->Encode(low, high - low, STRUCTURE_TOTAL);
        Accept(neighbor);
    }

    int TNeighborCoder::DecodeNeighbor(TRangeDecoder* decoder) {
        if (Remaining_ == 0) {
            ThrowCorrupted();
        }
        PrepareStep();
        if (RangeWeight_ <= 0) {
            ThrowCorrupted();
        }
        const uint64_t target = decoder->GetTarget(STRUCTURE_TOTAL);
        // Finds the first position whose upper boundary exceeds `target`.
        int position = 0;
        TNode accumulated;
        for (int step = TopBit_; step > 0; step >>= 1) {
            const int next = position + step;
            if (next > VertexCount_) {
                continue;
            }
            TNode candidate = accumulated;
            candidate.Weight += Tree_[next].Weight;
            candidate.Live += Tree_[next].Live;
            const int64_t weight = static_cast<int64_t>(candidate.Weight - Base_.Weight);
            if (GetBoundary(weight, candidate.Live - Base_.Live) <= target) {
                position = next;
                accumulated = candidate;
            }
        }
        const int neighbor = position;
        if (neighbor >= VertexCount_ || neighbor <= Previous_ || Weights_[neighbor] == 0) {
            ThrowCorrupted();
        }
        const int64_t weight = static_cast<int64_t>(accumulated.Weight - Base_.Weight);
        const int64_t live = accumulated.Live - Base_.Live;
        const uint64_t low = GetBoundary(weight, live);
        const uint64_t high = GetBoundary(weight + static_cast<int64_t>(Weights_[neighbor]), live + 1);
        if (target < low || target >= high) {
            ThrowCorrupted();
        }
        decoder->Consume(low, high - low);
        Accept(neighbor);
        return neighbor;
    }

    void TNeighborCoder::EndVertex() {
        for (int neighbor : Chosen_) {
            ++CodedEdgeCounts_[neighbor];
            SetWeight(neighbor, GetCurrentWeight(neighbor));
        }
        Chosen_.clear();
    }

    uint64_t TNeighborCoder::GetCurrentWeight(int vertex) const {
        if (CodedEdgeCounts_[vertex] >= Degrees_[vertex]) {
            return 0;
        }
        return CellWeights_[GetCell(Degrees_[vertex], CodedEdgeCounts_[vertex])];
    }

    void TNeighborCoder::SetWeight(int vertex, uint64_t weight) {
        const uint64_t deltaWeight = weight - Weights_[vertex];
        const int64_t deltaLive = (weight > 0 ? 1 : 0) - (Weights_[vertex] > 0 ? 1 : 0);
        if (deltaWeight == 0 && deltaLive == 0) {
            return;
        }
        Weights_[vertex] = weight;
        TotalWeight_ += deltaWeight;
        TotalLive_ += deltaLive;
        for (int i = vertex + 1; i <= VertexCount_; i += i & -i) {
            Tree_[i].Weight += deltaWeight;
            Tree_[i].Live += deltaLive;
        }
    }

    TNeighborCoder::TNode TNeighborCoder::GetPrefix(int end) const {
        TNode sum;
        for (int i = end; i > 0; i -= i & -i) {
            sum.Weight += Tree_[i].Weight;
            sum.Live += Tree_[i].Live;
        }
        return sum;
    }

    void TNeighborCoder::PrepareStep() {
        Base_ = GetPrefix(Previous_ + 1);
        RangeWeight_ = static_cast<int64_t>(TotalWeight_ - Base_.Weight);
        const int64_t rangeLive = TotalLive_ - Base_.Live;
        FreeTotal_ = STRUCTURE_TOTAL - UNITS_PER_CANDIDATE * static_cast<uint64_t>(rangeLive);
    }

    uint64_t TNeighborCoder::GetBoundary(int64_t weight, int64_t live) const {
        if (live <= 0) {
            return 0;
        }
        const uint64_t units = UNITS_PER_CANDIDATE * static_cast<uint64_t>(live);
        if (weight >= RangeWeight_) {
            return FreeTotal_ + units;
        }
        // Probability that the minimum of `Remaining_` draws lies in the prefix.
        const double tail = static_cast<double>(RangeWeight_ - weight) / static_cast<double>(RangeWeight_);
        const double cdf = 1.0 - std::pow(tail, static_cast<double>(Remaining_));
        const double scaled = std::max(0.0, std::floor(static_cast<double>(FreeTotal_) * cdf));
        const uint64_t continuous = std::min<uint64_t>(static_cast<uint64_t>(scaled), FreeTotal_);
        return continuous + units;
    }

    void TNeighborCoder::Accept(int neighbor) {
        Previous_ = neighbor;
        --Remaining_;
        Chosen_.push_back(neighbor);
    }
} // namespace NGraphCompressor
