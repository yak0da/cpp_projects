// Graph serializer.

#include "graph_codec.h"
#include "tsv_io.h"

#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
    enum class EMode {
        None,
        Serialize,
        Deserialize
    };

    struct TOptions {
        EMode Mode = EMode::None;
        std::string Input = {};
        std::string Output = {};
    };

    class TUsageException: public std::runtime_error {
    public:
        using std::runtime_error::runtime_error;
    };

    TOptions ParseOptions(int argc, char** argv) {
        TOptions options;
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument == "-s") {
                options.Mode = EMode::Serialize;
            } else if (argument == "-d") {
                options.Mode = EMode::Deserialize;
            } else if (argument == "-i" && i + 1 < argc) {
                options.Input = argv[++i];
            } else if (argument == "-o" && i + 1 < argc) {
                options.Output = argv[++i];
            } else {
                throw TUsageException("unexpected argument " + argument);
            }
        }
        if (options.Mode == EMode::None || options.Input.empty() || options.Output.empty()) {
            throw TUsageException("mode, input and output must be specified");
        }
        return options;
    }

    void Serialize(const TOptions& options) {
        const std::vector<NGraphCompressor::TEdge> edges = NGraphCompressor::ReadEdgesTsv(options.Input);
        NGraphCompressor::WriteBinaryFile(options.Output, NGraphCompressor::SerializeGraph(edges));
    }

    void Deserialize(const TOptions& options) {
        const std::vector<uint8_t> data = NGraphCompressor::ReadBinaryFile(options.Input);
        NGraphCompressor::WriteEdgesTsv(options.Output, NGraphCompressor::DeserializeGraph(data));
    }
} // anonymous namespace

int main(int argc, char** argv) {
    TOptions options;
    try {
        options = ParseOptions(argc, argv);
    } catch (const TUsageException& exception) {
        std::fprintf(stderr, "Error: %s\n", exception.what());
        std::fprintf(stderr, "Usage:\n");
        std::fprintf(stderr, "  %s -s -i input.tsv -o graph.bin\n", argv[0]);
        std::fprintf(stderr, "  %s -d -i graph.bin -o output.tsv\n", argv[0]);
        return 2;
    }
    try {
        if (options.Mode == EMode::Serialize) {
            Serialize(options);
        } else {
            Deserialize(options);
        }
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "Error: %s\n", exception.what());
        return 1;
    }
    return 0;
}
