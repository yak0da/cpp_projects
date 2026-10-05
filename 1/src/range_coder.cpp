#include "range_coder.h"

#include <bit>
#include <utility>

namespace NGraphCompressor {
    namespace {
        constexpr uint64_t TOP_VALUE = uint64_t{1} << 56;
        constexpr int STATE_BYTES = 9;
    } // anonymous namespace

    TRangeEncoder::TRangeEncoder()
        : Low_(0)
        , Carry_(false)
        , Range_(~uint64_t{0})
        , Cache_(0)
        , CacheSize_(1)
        , Output_()
    {
    }

    void TRangeEncoder::Encode(uint64_t start, uint64_t size, uint64_t total) {
        const uint64_t step = Range_ / total;
        const uint64_t addition = step * start;
        Low_ += addition;
        if (Low_ < addition) {
            Carry_ = true;
        }
        Range_ = step * size;
        while (Range_ < TOP_VALUE) {
            Range_ <<= 8;
            ShiftLow();
        }
    }

    void TRangeEncoder::EncodeBits(uint64_t value, int bits) {
        if (bits == 0) {
            return;
        }
        Encode(value & ((uint64_t{1} << bits) - 1), 1, uint64_t{1} << bits);
    }

    void TRangeEncoder::EncodeUniversal(uint64_t value) {
        // Values up to 2^64 - 2 are supported.
        const uint64_t shifted = value + 1;
        const int width = static_cast<int>(std::bit_width(shifted));
        EncodeBits(width - 1, 6);
        if (width - 1 > 32) {
            EncodeBits(shifted >> 32, width - 1 - 32);
            EncodeBits(shifted, 32);
        } else {
            EncodeBits(shifted, width - 1);
        }
    }

    std::vector<uint8_t> TRangeEncoder::Finish() {
        for (int i = 0; i < STATE_BYTES; ++i) {
            ShiftLow();
        }
        return std::move(Output_);
    }

    void TRangeEncoder::ShiftLow() {
        if (Low_ < 0xFF00000000000000ULL || Carry_) {
            const uint8_t carry = Carry_ ? 1 : 0;
            uint8_t byte = Cache_;
            do {
                Output_.push_back(static_cast<uint8_t>(byte + carry));
                byte = 0xFF;
            } while (--CacheSize_ != 0);
            Cache_ = static_cast<uint8_t>(Low_ >> 56);
        }
        ++CacheSize_;
        Low_ <<= 8;
        Carry_ = false;
    }

    TRangeDecoder::TRangeDecoder(const uint8_t* data, size_t size)
        : Data_(data)
        , Size_(size)
        , Position_(0)
        , Code_(0)
        , Range_(~uint64_t{0})
        , Step_(1)
    {
        for (int i = 0; i < STATE_BYTES; ++i) {
            Code_ = (Code_ << 8) | NextByte();
        }
    }

    uint64_t TRangeDecoder::GetTarget(uint64_t total) {
        Step_ = Range_ / total;
        const uint64_t target = Code_ / Step_;
        return target < total ? target : total - 1;
    }

    void TRangeDecoder::Consume(uint64_t start, uint64_t size) {
        Code_ -= Step_ * start;
        Range_ = Step_ * size;
        while (Range_ < TOP_VALUE) {
            Code_ = (Code_ << 8) | NextByte();
            Range_ <<= 8;
        }
    }

    uint64_t TRangeDecoder::DecodeBits(int bits) {
        if (bits == 0) {
            return 0;
        }
        const uint64_t value = GetTarget(uint64_t{1} << bits);
        Consume(value, 1);
        return value;
    }

    uint64_t TRangeDecoder::DecodeUniversal() {
        const int extraBits = static_cast<int>(DecodeBits(6));
        uint64_t shifted = 0;
        if (extraBits > 32) {
            const uint64_t high = DecodeBits(extraBits - 32);
            shifted = (uint64_t{1} << extraBits) | (high << 32) | DecodeBits(32);
        } else {
            shifted = (uint64_t{1} << extraBits) | DecodeBits(extraBits);
        }
        return shifted - 1;
    }

    uint8_t TRangeDecoder::NextByte() {
        return Position_ < Size_ ? Data_[Position_++] : 0;
    }
} // namespace NGraphCompressor
