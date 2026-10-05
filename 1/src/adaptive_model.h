#pragma once

#include "range_coder.h"

#include <cstdint>
#include <vector>

namespace NGraphCompressor {
    // Adaptive frequency model over the alphabet [0, symbolCount).
    class TAdaptiveModel {
    public:
        explicit TAdaptiveModel(int symbolCount);

        void Encode(int symbol, TRangeEncoder* encoder);
        int Decode(TRangeDecoder* decoder);

    private:
        uint64_t Prefix(int end) const;
        void Add(int symbol, uint64_t delta);

    private:
        int SymbolCount_;
        int TopBit_;
        std::vector<uint64_t> Tree_;
        std::vector<uint64_t> Frequencies_;
        uint64_t Total_;
    };

    // Codes non-negative integers: small values through an adaptive model and
    // larger ones through an escape symbol followed by a universal code.
    class TEscapedIntModel {
    public:
        explicit TEscapedIntModel(int directSymbolCount);

        void Encode(uint64_t value, TRangeEncoder* encoder);
        uint64_t Decode(TRangeDecoder* decoder);

    private:
        int DirectSymbolCount_;
        TAdaptiveModel Model_;
    };
} // namespace NGraphCompressor
