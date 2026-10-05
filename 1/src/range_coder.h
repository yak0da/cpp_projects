#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace NGraphCompressor {
    // Largest `total` accepted by the coders. The range never drops below 2^56,
    // so every interval keeps at least 16 bits of precision.
    inline constexpr uint64_t MAX_CODER_TOTAL = uint64_t{1} << 40;

    // Range encoder with a 64-bit range and byte-wise carry propagation.
    class TRangeEncoder {
    public:
        TRangeEncoder();

        // Encodes the interval [start, start + size) out of [0, total).
        void Encode(uint64_t start, uint64_t size, uint64_t total);
        // Encodes `bits` (at most 40) low bits of `value` with uniform probability.
        void EncodeBits(uint64_t value, int bits);
        // Encodes an arbitrary 64-bit value (Elias gamma style).
        void EncodeUniversal(uint64_t value);
        // Flushes the state and returns the produced bytes.
        std::vector<uint8_t> Finish();

    private:
        void ShiftLow();

    private:
        uint64_t Low_;
        bool Carry_;
        uint64_t Range_;
        uint8_t Cache_;
        uint64_t CacheSize_;
        std::vector<uint8_t> Output_;
    };

    class TRangeDecoder {
    public:
        TRangeDecoder(const uint8_t* data, size_t size);

        // Returns a value in [0, total) that identifies the encoded interval.
        // Must be followed by Consume() with the same `total`.
        uint64_t GetTarget(uint64_t total);
        // Removes the interval [start, start + size) found after GetTarget().
        void Consume(uint64_t start, uint64_t size);
        uint64_t DecodeBits(int bits);
        uint64_t DecodeUniversal();

    private:
        uint8_t NextByte();

    private:
        const uint8_t* Data_;
        size_t Size_;
        size_t Position_;
        uint64_t Code_;
        uint64_t Range_;
        uint64_t Step_;
    };
} // namespace NGraphCompressor
