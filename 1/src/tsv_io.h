#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace NGraphCompressor {
    struct TEdge {
        uint32_t From = 0;
        uint32_t To = 0;
        uint8_t Weight = 0;
    };

    // All functions throw TGraphCompressorException on failure.
    std::vector<TEdge> ReadEdgesTsv(const std::string& path);
    void WriteEdgesTsv(const std::string& path, const std::vector<TEdge>& edges);

    std::vector<uint8_t> ReadBinaryFile(const std::string& path);
    void WriteBinaryFile(const std::string& path, const std::vector<uint8_t>& data);
} // namespace NGraphCompressor
