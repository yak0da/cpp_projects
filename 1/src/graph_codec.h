#pragma once

#include "tsv_io.h"

#include <cstdint>
#include <vector>

namespace NGraphCompressor {
    // Compresses an undirected weighted graph without multiple edges.
    // Throws TGraphCompressorException if the input violates this assumption.
    std::vector<uint8_t> SerializeGraph(const std::vector<TEdge>& edges);

    // Restores the edges written by SerializeGraph(). Their order and the order
    // of endpoints within an edge may differ from the original ones.
    // Throws TGraphCompressorException if the data is corrupted.
    std::vector<TEdge> DeserializeGraph(const std::vector<uint8_t>& data);
} // namespace NGraphCompressor
