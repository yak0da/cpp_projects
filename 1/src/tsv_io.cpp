#include "tsv_io.h"

#include "exception.h"

#include <cstdio>
#include <memory>

namespace NGraphCompressor {
    namespace {
        constexpr size_t IO_CHUNK_SIZE = 1 << 20;
        constexpr uint64_t MAX_ID = 0xFFFFFFFFULL;
        constexpr uint64_t MAX_WEIGHT = 255;

        class TFileCloser {
        public:
            void operator()(std::FILE* file) const {
                std::fclose(file);
            }
        };

        using TFilePtr = std::unique_ptr<std::FILE, TFileCloser>;

        TFilePtr OpenFile(const std::string& path, const char* mode) {
            TFilePtr file(std::fopen(path.c_str(), mode));
            if (!file) {
                throw TGraphCompressorException("cannot open " + path);
            }
            return file;
        }

        void WriteChunk(std::FILE* file, const char* data, size_t size, const std::string& path) {
            if (std::fwrite(data, 1, size, file) != size) {
                throw TGraphCompressorException("cannot write " + path);
            }
        }

        bool IsDigit(uint8_t symbol) {
            return symbol >= '0' && symbol <= '9';
        }

        bool IsBlank(uint8_t symbol) {
            return symbol == ' ' || symbol == '\t';
        }

        // Parses an unsigned number starting at `*position`, skipping leading
        // spaces and tabs. Returns false if there is no number up to `maxValue`.
        bool TryParseNumber(const std::vector<uint8_t>& text, size_t* position, uint64_t maxValue, uint64_t* value) {
            size_t i = *position;
            while (i < text.size() && IsBlank(text[i])) {
                ++i;
            }
            if (i >= text.size() || !IsDigit(text[i])) {
                return false;
            }
            uint64_t result = 0;
            while (i < text.size() && IsDigit(text[i])) {
                result = result * 10 + static_cast<uint64_t>(text[i] - '0');
                if (result > maxValue) {
                    return false;
                }
                ++i;
            }
            *position = i;
            *value = result;
            return true;
        }

        char* AppendNumber(uint32_t value, char* output) {
            char digits[10];
            int length = 0;
            do {
                digits[length++] = static_cast<char>('0' + value % 10);
                value /= 10;
            } while (value != 0);
            while (length > 0) {
                *output++ = digits[--length];
            }
            return output;
        }
    } // anonymous namespace

    std::vector<TEdge> ReadEdgesTsv(const std::string& path) {
        const std::vector<uint8_t> text = ReadBinaryFile(path);
        std::vector<TEdge> edges;
        edges.reserve(text.size() / 24);
        size_t position = 0;
        size_t line = 1;
        while (position < text.size()) {
            // Skips empty lines.
            if (text[position] == '\n' || text[position] == '\r') {
                if (text[position] == '\n') {
                    ++line;
                }
                ++position;
                continue;
            }
            uint64_t from = 0;
            uint64_t to = 0;
            uint64_t weight = 0;
            const bool parsed = TryParseNumber(text, &position, MAX_ID, &from) &&
                TryParseNumber(text, &position, MAX_ID, &to) &&
                TryParseNumber(text, &position, MAX_WEIGHT, &weight);
            while (position < text.size() && (IsBlank(text[position]) || text[position] == '\r')) {
                ++position;
            }
            if (!parsed || (position < text.size() && text[position] != '\n')) {
                throw TGraphCompressorException("malformed line " + std::to_string(line) + " in " + path);
            }
            edges.push_back({static_cast<uint32_t>(from), static_cast<uint32_t>(to), static_cast<uint8_t>(weight)});
        }
        return edges;
    }

    void WriteEdgesTsv(const std::string& path, const std::vector<TEdge>& edges) {
        const TFilePtr file = OpenFile(path, "wb");
        // One line takes at most 10 + 1 + 10 + 1 + 3 + 1 characters.
        std::vector<char> buffer(IO_CHUNK_SIZE + 32);
        char* output = buffer.data();
        for (const TEdge& edge : edges) {
            output = AppendNumber(edge.From, output);
            *output++ = '\t';
            output = AppendNumber(edge.To, output);
            *output++ = '\t';
            output = AppendNumber(edge.Weight, output);
            *output++ = '\n';
            if (static_cast<size_t>(output - buffer.data()) >= IO_CHUNK_SIZE) {
                WriteChunk(file.get(), buffer.data(), output - buffer.data(), path);
                output = buffer.data();
            }
        }
        WriteChunk(file.get(), buffer.data(), output - buffer.data(), path);
    }

    std::vector<uint8_t> ReadBinaryFile(const std::string& path) {
        const TFilePtr file = OpenFile(path, "rb");
        std::vector<uint8_t> data;
        std::vector<uint8_t> chunk(IO_CHUNK_SIZE);
        size_t readSize = 0;
        while ((readSize = std::fread(chunk.data(), 1, chunk.size(), file.get())) > 0) {
            data.insert(data.end(), chunk.begin(), chunk.begin() + readSize);
        }
        if (std::ferror(file.get())) {
            throw TGraphCompressorException("cannot read " + path);
        }
        return data;
    }

    void WriteBinaryFile(const std::string& path, const std::vector<uint8_t>& data) {
        const TFilePtr file = OpenFile(path, "wb");
        WriteChunk(file.get(), reinterpret_cast<const char*>(data.data()), data.size(), path);
    }
} // namespace NGraphCompressor
