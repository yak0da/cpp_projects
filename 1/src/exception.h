#pragma once

#include <stdexcept>

namespace NGraphCompressor {
    // Thrown on malformed input files and I/O failures.
    class TGraphCompressorException: public std::runtime_error {
    public:
        using std::runtime_error::runtime_error;
    };
} // namespace NGraphCompressor
