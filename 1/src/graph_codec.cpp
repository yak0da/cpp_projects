#include "graph_codec.h"

#include "adaptive_model.h"
#include "exception.h"
#include "neighbor_coder.h"
#include "range_coder.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace NGraphCompressor {
    namespace {
        // Binary layout:
        //   4 bytes   number of edges E (little endian);
        //   E bytes   edge weights in canonical order;
        //   the rest  range coded stream with vertex ids, degrees, loops and
        //             the adjacency structure.
        constexpr size_t HEADER_SIZE = 4;
        constexpr int GAP_DIRECT_SYMBOL_COUNT = 256;
        constexpr int DEGREE_DIRECT_SYMBOL_COUNT = 512;
        constexpr uint64_t ID_SPACE_SIZE = uint64_t{1} << 32;
        // Cell weights are stored as 2^(level / WEIGHT_STEPS) with level in
        // [0, MAX_WEIGHT_LEVEL].
        constexpr int WEIGHT_STEPS = 16;
        constexpr int MAX_WEIGHT_LEVEL = 383;
        constexpr int WEIGHT_LEVEL_BITS = 9;
        constexpr int FIT_ITERATION_COUNT = 3;

        struct TGraphStructure {
            // Vertex ids in increasing order; the "rank" of a vertex is its index here.
            std::vector<uint32_t> Ids;
            // Degree of each rank, loops excluded.
            std::vector<uint32_t> Degrees;
            // Ranks of the vertices with loops, increasing.
            std::vector<uint32_t> Loops;
        };

        // Adjacency restricted to neighbors with larger positions (CSR layout).
        struct TLaterNeighbors {
            std::vector<uint64_t> Offsets;
            std::vector<int> Targets;
            std::vector<uint8_t> Weights;
        };

        [[noreturn]] void ThrowCorrupted() {
            throw TGraphCompressorException("corrupted input");
        }

        // Vertices are coded in order of decreasing degree; ties go by rank.
        std::vector<int> GetCodingOrder(const std::vector<uint32_t>& degrees) {
            std::vector<int> order(degrees.size());
            std::iota(order.begin(), order.end(), 0);
            std::stable_sort(order.begin(), order.end(), [&degrees](int left, int right) {
                return degrees[left] > degrees[right];
            });
            return order;
        }

        int GetGapShift(size_t idCount) {
            if (idCount == 0) {
                return 0;
            }
            int shift = 0;
            while ((ID_SPACE_SIZE >> (shift + 1)) >= idCount) {
                ++shift;
            }
            // Now 2^shift approximates the mean gap; the quotient keeps ~4 bits.
            return std::max(0, shift - 4);
        }

        void EncodeIds(const std::vector<uint32_t>& ids, TRangeEncoder* encoder) {
            const int shift = GetGapShift(ids.size());
            TEscapedIntModel quotientModel(GAP_DIRECT_SYMBOL_COUNT);
            int64_t previous = -1;
            for (uint32_t id : ids) {
                const uint64_t gap = static_cast<uint64_t>(id - previous - 1);
                quotientModel.Encode(gap >> shift, encoder);
                encoder->EncodeBits(gap, shift);
                previous = id;
            }
        }

        std::vector<uint32_t> DecodeIds(size_t idCount, TRangeDecoder* decoder) {
            const int shift = GetGapShift(idCount);
            TEscapedIntModel quotientModel(GAP_DIRECT_SYMBOL_COUNT);
            std::vector<uint32_t> ids(idCount);
            int64_t previous = -1;
            for (size_t i = 0; i < idCount; ++i) {
                const uint64_t quotient = quotientModel.Decode(decoder);
                const uint64_t gap = (quotient << shift) | decoder->DecodeBits(shift);
                const int64_t id = previous + 1 + static_cast<int64_t>(gap);
                if (id >= static_cast<int64_t>(ID_SPACE_SIZE)) {
                    ThrowCorrupted();
                }
                ids[i] = static_cast<uint32_t>(id);
                previous = id;
            }
            return ids;
        }

        // Cells that can ever be used by vertices with the given degrees.
        std::vector<int> GetCanonicalCells(const std::vector<uint32_t>& degrees) {
            std::vector<uint32_t> maxDegrees(DEGREE_BUCKET_COUNT, 0);
            for (uint32_t degree : degrees) {
                const int bucket = GetDegreeBucket(degree);
                maxDegrees[bucket] = std::max(maxDegrees[bucket], degree);
            }
            std::vector<int> cells;
            for (int bucket = 0; bucket < DEGREE_BUCKET_COUNT; ++bucket) {
                const uint32_t levelCount = std::min<uint32_t>(maxDegrees[bucket], CODED_LEVEL_COUNT);
                for (uint32_t level = 0; level < levelCount; ++level) {
                    cells.push_back(bucket * CODED_LEVEL_COUNT + static_cast<int>(level));
                }
            }
            return cells;
        }

        uint64_t GetWeightFromLevel(int level) {
            const double weight = std::exp2(static_cast<double>(level) / WEIGHT_STEPS);
            return std::max<uint64_t>(1, static_cast<uint64_t>(std::llround(weight)));
        }

        // Fits cell weights by maximizing the likelihood of the neighbor choices
        // (multiplicative EM updates). Returns quantized weight levels; `usedCells`
        // marks the cells that are ever chosen.
        std::vector<int> FitCellWeights(
            const std::vector<uint32_t>& degrees,
            const TLaterNeighbors& later,
            std::vector<bool>* usedCells)
        {
            const int vertexCount = static_cast<int>(degrees.size());
            std::vector<double> cellWeights(CELL_COUNT, 1.0);
            std::vector<uint32_t> bucketDegrees(DEGREE_BUCKET_COUNT, 1);
            for (uint32_t degree : degrees) {
                bucketDegrees[GetDegreeBucket(degree)] = degree;
            }
            for (int cell = 0; cell < CELL_COUNT; ++cell) {
                const double degree = bucketDegrees[cell / CODED_LEVEL_COUNT];
                cellWeights[cell] = std::max(1.0, degree - cell % CODED_LEVEL_COUNT);
            }

            std::vector<double> chosenCounts(CELL_COUNT, 0.0);
            for (int iteration = 0; iteration < FIT_ITERATION_COUNT; ++iteration) {
                std::vector<uint32_t> codedEdgeCounts(vertexCount, 0);
                std::vector<double> liveCounts(CELL_COUNT, 0.0);
                // Time integral of liveCounts[cell] * k / S, maintained lazily.
                std::vector<double> integrals(CELL_COUNT, 0.0);
                std::vector<double> stamps(CELL_COUNT, 0.0);
                double clock = 0.0;
                double liveWeight = 0.0;
                std::fill(chosenCounts.begin(), chosenCounts.end(), 0.0);
                auto changeLiveCount = [&](int cell, double delta) {
                    integrals[cell] += liveCounts[cell] * (clock - stamps[cell]);
                    stamps[cell] = clock;
                    liveCounts[cell] += delta;
                    liveWeight += delta * cellWeights[cell];
                };
                for (int vertex = 0; vertex < vertexCount; ++vertex) {
                    if (degrees[vertex] > 0) {
                        changeLiveCount(GetCell(degrees[vertex], 0), 1.0);
                    }
                }
                for (int vertex = 0; vertex < vertexCount; ++vertex) {
                    if (codedEdgeCounts[vertex] < degrees[vertex]) {
                        changeLiveCount(GetCell(degrees[vertex], codedEdgeCounts[vertex]), -1.0);
                    }
                    const uint64_t begin = later.Offsets[vertex];
                    const uint64_t end = later.Offsets[vertex + 1];
                    if (begin == end) {
                        continue;
                    }
                    if (liveWeight > 0) {
                        clock += static_cast<double>(end - begin) / liveWeight;
                    }
                    for (uint64_t i = begin; i < end; ++i) {
                        const int neighbor = later.Targets[i];
                        const int cell = GetCell(degrees[neighbor], codedEdgeCounts[neighbor]);
                        chosenCounts[cell] += 1.0;
                        changeLiveCount(cell, -1.0);
                        ++codedEdgeCounts[neighbor];
                        if (codedEdgeCounts[neighbor] < degrees[neighbor]) {
                            changeLiveCount(GetCell(degrees[neighbor], codedEdgeCounts[neighbor]), 1.0);
                        }
                    }
                }
                for (int cell = 0; cell < CELL_COUNT; ++cell) {
                    const double integral = integrals[cell] + liveCounts[cell] * (clock - stamps[cell]);
                    const double expected = cellWeights[cell] * integral;
                    if (chosenCounts[cell] > 0 && expected > 0) {
                        cellWeights[cell] *= chosenCounts[cell] / expected;
                    }
                }
            }

            usedCells->assign(CELL_COUNT, false);
            double maxWeight = 0.0;
            for (int cell = 0; cell < CELL_COUNT; ++cell) {
                if (chosenCounts[cell] > 0) {
                    (*usedCells)[cell] = true;
                    maxWeight = std::max(maxWeight, cellWeights[cell]);
                }
            }
            std::vector<int> levels(CELL_COUNT, 0);
            for (int cell = 0; cell < CELL_COUNT; ++cell) {
                if (!(*usedCells)[cell]) {
                    continue;
                }
                const double level = std::round(WEIGHT_STEPS * std::log2(cellWeights[cell] / maxWeight)) +
                    MAX_WEIGHT_LEVEL;
                levels[cell] = static_cast<int>(std::clamp(level, 0.0, static_cast<double>(MAX_WEIGHT_LEVEL)));
            }
            return levels;
        }

        void PutUint32(uint32_t value, std::vector<uint8_t>* data) {
            for (int i = 0; i < 4; ++i) {
                data->push_back((value >> (8 * i)) & 0xFF);
            }
        }

        uint32_t GetUint32(const uint8_t* data) {
            uint32_t value = 0;
            for (int i = 0; i < 4; ++i) {
                value |= static_cast<uint32_t>(data[i]) << (8 * i);
            }
            return value;
        }
    } // anonymous namespace

    std::vector<uint8_t> SerializeGraph(const std::vector<TEdge>& edges) {
        TGraphStructure graph;
        graph.Ids.reserve(edges.size() * 2);
        for (const TEdge& edge : edges) {
            graph.Ids.push_back(edge.From);
            graph.Ids.push_back(edge.To);
        }
        std::sort(graph.Ids.begin(), graph.Ids.end());
        graph.Ids.erase(std::unique(graph.Ids.begin(), graph.Ids.end()), graph.Ids.end());
        const int vertexCount = static_cast<int>(graph.Ids.size());
        auto getRank = [&graph](uint32_t id) {
            return static_cast<int>(std::lower_bound(graph.Ids.begin(), graph.Ids.end(), id) - graph.Ids.begin());
        };

        // Ranks of edge endpoints; loops are kept apart.
        std::vector<std::pair<int, int>> rankedEdges(edges.size());
        std::vector<std::pair<uint32_t, uint8_t>> loops;
        graph.Degrees.assign(vertexCount, 0);
        for (size_t i = 0; i < edges.size(); ++i) {
            const int fromRank = getRank(edges[i].From);
            const int toRank = getRank(edges[i].To);
            rankedEdges[i] = {fromRank, toRank};
            if (fromRank == toRank) {
                loops.push_back({static_cast<uint32_t>(fromRank), edges[i].Weight});
            } else {
                ++graph.Degrees[fromRank];
                ++graph.Degrees[toRank];
            }
        }
        std::sort(loops.begin(), loops.end());
        for (size_t i = 1; i < loops.size(); ++i) {
            if (loops[i].first == loops[i - 1].first) {
                throw TGraphCompressorException("multiple loops at one vertex are not supported");
            }
        }

        const std::vector<int> order = GetCodingOrder(graph.Degrees);
        std::vector<int> positions(vertexCount);
        std::vector<uint32_t> positionDegrees(vertexCount);
        for (int i = 0; i < vertexCount; ++i) {
            positions[order[i]] = i;
            positionDegrees[i] = graph.Degrees[order[i]];
        }

        TLaterNeighbors later;
        later.Offsets.assign(vertexCount + 1, 0);
        for (const auto& [fromRank, toRank] : rankedEdges) {
            if (fromRank != toRank) {
                ++later.Offsets[std::min(positions[fromRank], positions[toRank]) + 1];
            }
        }
        for (int i = 0; i < vertexCount; ++i) {
            later.Offsets[i + 1] += later.Offsets[i];
        }
        const uint64_t simpleEdgeCount = later.Offsets[vertexCount];
        {
            std::vector<std::pair<int, uint8_t>> slots(simpleEdgeCount);
            std::vector<uint64_t> fillPositions(later.Offsets.begin(), later.Offsets.end() - 1);
            for (size_t i = 0; i < edges.size(); ++i) {
                const auto [fromRank, toRank] = rankedEdges[i];
                if (fromRank == toRank) {
                    continue;
                }
                const int low = std::min(positions[fromRank], positions[toRank]);
                const int high = std::max(positions[fromRank], positions[toRank]);
                slots[fillPositions[low]++] = {high, edges[i].Weight};
            }
            rankedEdges.clear();
            rankedEdges.shrink_to_fit();
            for (int vertex = 0; vertex < vertexCount; ++vertex) {
                const auto begin = slots.begin() + later.Offsets[vertex];
                const auto end = slots.begin() + later.Offsets[vertex + 1];
                std::sort(begin, end);
                auto isSameTarget = [](const auto& left, const auto& right) {
                    return left.first == right.first;
                };
                if (std::adjacent_find(begin, end, isSameTarget) != end) {
                    throw TGraphCompressorException("multiple edges are not supported");
                }
            }
            later.Targets.resize(simpleEdgeCount);
            later.Weights.resize(simpleEdgeCount);
            for (uint64_t i = 0; i < simpleEdgeCount; ++i) {
                later.Targets[i] = slots[i].first;
                later.Weights[i] = slots[i].second;
            }
        }

        std::vector<uint8_t> data;
        data.reserve(HEADER_SIZE + edges.size() * 4);
        PutUint32(static_cast<uint32_t>(edges.size()), &data);
        data.insert(data.end(), later.Weights.begin(), later.Weights.end());
        for (const auto& loop : loops) {
            data.push_back(loop.second);
        }

        TRangeEncoder encoder;
        encoder.EncodeUniversal(static_cast<uint64_t>(vertexCount));
        EncodeIds(graph.Ids, &encoder);
        TEscapedIntModel degreeModel(DEGREE_DIRECT_SYMBOL_COUNT);
        for (uint32_t degree : graph.Degrees) {
            degreeModel.Encode(degree, &encoder);
        }
        encoder.EncodeUniversal(loops.size());
        uint32_t previousLoop = 0;
        for (const auto& loop : loops) {
            encoder.EncodeUniversal(loop.first - previousLoop);
            previousLoop = loop.first;
        }

        std::vector<bool> usedCells;
        const std::vector<int> levels = FitCellWeights(positionDegrees, later, &usedCells);
        std::vector<uint64_t> cellWeights(CELL_COUNT, 1);
        for (int cell : GetCanonicalCells(positionDegrees)) {
            encoder.EncodeBits(usedCells[cell] ? 1 : 0, 1);
            if (usedCells[cell]) {
                encoder.EncodeBits(static_cast<uint64_t>(levels[cell]), WEIGHT_LEVEL_BITS);
                cellWeights[cell] = GetWeightFromLevel(levels[cell]);
            }
        }

        TNeighborCoder coder(positionDegrees, cellWeights);
        for (int vertex = 0; vertex < vertexCount; ++vertex) {
            const uint32_t neighborCount = coder.BeginVertex(vertex);
            const uint64_t begin = later.Offsets[vertex];
            if (neighborCount != later.Offsets[vertex + 1] - begin) {
                throw TGraphCompressorException("internal error: neighbor count mismatch");
            }
            for (uint32_t i = 0; i < neighborCount; ++i) {
                coder.EncodeNeighbor(later.Targets[begin + i], &encoder);
            }
            coder.EndVertex();
        }

        const std::vector<uint8_t> stream = encoder.Finish();
        data.insert(data.end(), stream.begin(), stream.end());
        return data;
    }

    std::vector<TEdge> DeserializeGraph(const std::vector<uint8_t>& data) {
        if (data.size() < HEADER_SIZE) {
            ThrowCorrupted();
        }
        const uint64_t edgeCount = GetUint32(data.data());
        if (data.size() < HEADER_SIZE + edgeCount) {
            ThrowCorrupted();
        }
        const uint8_t* weights = data.data() + HEADER_SIZE;
        TRangeDecoder decoder(weights + edgeCount, data.size() - HEADER_SIZE - edgeCount);

        TGraphStructure graph;
        const uint64_t vertexCount = decoder.DecodeUniversal();
        if (vertexCount > 2 * edgeCount) {
            ThrowCorrupted();
        }
        graph.Ids = DecodeIds(vertexCount, &decoder);
        TEscapedIntModel degreeModel(DEGREE_DIRECT_SYMBOL_COUNT);
        graph.Degrees.resize(vertexCount);
        uint64_t degreeSum = 0;
        for (uint64_t i = 0; i < vertexCount; ++i) {
            const uint64_t degree = degreeModel.Decode(&decoder);
            if (degree >= vertexCount) {
                ThrowCorrupted();
            }
            graph.Degrees[i] = static_cast<uint32_t>(degree);
            degreeSum += degree;
        }
        const uint64_t loopCount = decoder.DecodeUniversal();
        if (loopCount > edgeCount || loopCount > vertexCount || degreeSum != 2 * (edgeCount - loopCount)) {
            ThrowCorrupted();
        }
        const uint64_t simpleEdgeCount = edgeCount - loopCount;
        graph.Loops.resize(loopCount);
        uint64_t previousLoop = 0;
        for (uint64_t i = 0; i < loopCount; ++i) {
            previousLoop += decoder.DecodeUniversal();
            if (previousLoop >= vertexCount) {
                ThrowCorrupted();
            }
            graph.Loops[i] = static_cast<uint32_t>(previousLoop);
        }

        const std::vector<int> order = GetCodingOrder(graph.Degrees);
        std::vector<uint32_t> positionDegrees(vertexCount);
        for (uint64_t i = 0; i < vertexCount; ++i) {
            positionDegrees[i] = graph.Degrees[order[i]];
        }
        std::vector<uint64_t> cellWeights(CELL_COUNT, 1);
        for (int cell : GetCanonicalCells(positionDegrees)) {
            if (decoder.DecodeBits(1) != 0) {
                cellWeights[cell] = GetWeightFromLevel(static_cast<int>(decoder.DecodeBits(WEIGHT_LEVEL_BITS)));
            }
        }

        std::vector<TEdge> edges;
        edges.reserve(edgeCount);
        TNeighborCoder coder(positionDegrees, cellWeights);
        for (uint64_t vertex = 0; vertex < vertexCount; ++vertex) {
            const uint32_t neighborCount = coder.BeginVertex(static_cast<int>(vertex));
            for (uint32_t i = 0; i < neighborCount; ++i) {
                const int neighbor = coder.DecodeNeighbor(&decoder);
                if (edges.size() >= simpleEdgeCount) {
                    ThrowCorrupted();
                }
                edges.push_back({graph.Ids[order[vertex]], graph.Ids[order[neighbor]], weights[edges.size()]});
            }
            coder.EndVertex();
        }
        if (edges.size() != simpleEdgeCount) {
            ThrowCorrupted();
        }
        for (uint32_t loop : graph.Loops) {
            edges.push_back({graph.Ids[loop], graph.Ids[loop], weights[edges.size()]});
        }
        return edges;
    }
} // namespace NGraphCompressor
