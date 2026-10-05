#include "adaptive_model.h"

namespace NGraphCompressor {
    namespace {
        constexpr uint64_t INITIAL_COUNT = 1;
        constexpr uint64_t INCREMENT = 32;
    } // anonymous namespace

    TAdaptiveModel::TAdaptiveModel(int symbolCount)
        : SymbolCount_(symbolCount)
        , TopBit_(1)
        , Tree_(symbolCount + 1, 0)
        , Frequencies_(symbolCount, 0)
        , Total_(0)
    {
        while (TopBit_ * 2 <= SymbolCount_) {
            TopBit_ *= 2;
        }
        for (int symbol = 0; symbol < SymbolCount_; ++symbol) {
            Add(symbol, INITIAL_COUNT);
        }
    }

    void TAdaptiveModel::Encode(int symbol, TRangeEncoder* encoder) {
        encoder->Encode(Prefix(symbol), Frequencies_[symbol], Total_);
        Add(symbol, INCREMENT);
    }

    int TAdaptiveModel::Decode(TRangeDecoder* decoder) {
        const uint64_t target = decoder->GetTarget(Total_);
        // Finds the symbol whose cumulative interval contains `target`.
        int position = 0;
        uint64_t accumulated = 0;
        for (int step = TopBit_; step > 0; step >>= 1) {
            const int next = position + step;
            if (next <= SymbolCount_ && accumulated + Tree_[next] <= target) {
                position = next;
                accumulated += Tree_[next];
            }
        }
        decoder->Consume(accumulated, Frequencies_[position]);
        Add(position, INCREMENT);
        return position;
    }

    uint64_t TAdaptiveModel::Prefix(int end) const {
        uint64_t sum = 0;
        for (int i = end; i > 0; i -= i & -i) {
            sum += Tree_[i];
        }
        return sum;
    }

    void TAdaptiveModel::Add(int symbol, uint64_t delta) {
        Frequencies_[symbol] += delta;
        Total_ += delta;
        for (int i = symbol + 1; i <= SymbolCount_; i += i & -i) {
            Tree_[i] += delta;
        }
    }

    TEscapedIntModel::TEscapedIntModel(int directSymbolCount)
        : DirectSymbolCount_(directSymbolCount)
        , Model_(directSymbolCount + 1)
    {
    }

    void TEscapedIntModel::Encode(uint64_t value, TRangeEncoder* encoder) {
        if (value < static_cast<uint64_t>(DirectSymbolCount_)) {
            Model_.Encode(static_cast<int>(value), encoder);
            return;
        }
        Model_.Encode(DirectSymbolCount_, encoder);
        encoder->EncodeUniversal(value - DirectSymbolCount_);
    }

    uint64_t TEscapedIntModel::Decode(TRangeDecoder* decoder) {
        const int symbol = Model_.Decode(decoder);
        if (symbol < DirectSymbolCount_) {
            return symbol;
        }
        return DirectSymbolCount_ + decoder->DecodeUniversal();
    }
} // namespace NGraphCompressor
