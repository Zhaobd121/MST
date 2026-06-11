#include <iostream>
#include <fstream>
#include <filesystem>  // C++17
#include "FCI.h"
#include "Timer.h"
#include "utility.hpp"
#include <filesystem>

void write_clusters_to_file(const std::vector<Cluster>& clusters, const std::string& file_path) {
    std::ofstream out(file_path);
    if (!out) {
        std::cerr << "Failed to open file: " << file_path << std::endl;
        return;
    }

    for (const auto& cluster : clusters) {
        out << "Cluster rep: " << cluster.rep << "\n";
        out << "  Cores: ";
        for (int core : cluster.cores) {
            out << core << " ";
        }
        out << "\n  Non-Cores: ";
        for (int nc : cluster.nonCores) {
            out << nc << " ";
        }
        out << "\n\n";
    }

    out.close();
}

void write_edge_binary(std::ofstream& file, int u, int v, int num, int denom) {
    file.write(reinterpret_cast<const char*>(&u), sizeof(int));
    file.write(reinterpret_cast<const char*>(&v), sizeof(int));
    file.write(reinterpret_cast<const char*>(&num), sizeof(int));
    file.write(reinterpret_cast<const char*>(&denom), sizeof(int));
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " index <miu> <input_file>\n"
                  << "  " << argv[0] << " query <miu> <epsilon> <input_file>\n";
        return 1;
    }

    std::string mode = argv[1];

    if (mode == "index") {
        if (argc != 4) {
            std::cerr << "Usage: " << argv[0] << " index <miu> <input_file>\n";
            return 1;
        }
        int miu = std::stoi(argv[2]);
        std::string input_file = argv[3];
        std::string filename = std::filesystem::path(input_file).stem().string();
        std::string output_dir = "index/" + filename + "_" + std::to_string(miu);
        if (!std::filesystem::exists(output_dir)) {
            std::filesystem::create_directories(output_dir);
        }
        FCI fci;
        fci.read_dataset(input_file);
        std::cout << "finished reading graph" << std::endl;
        Timer t;
        fci.bulid_core_index(miu);
        std::cout << "Core index build time: " << Utility::integer_to_string(t.elapsed()) << std::endl;

        fci.build_spanning_tree();
        std::cout << "Spanning tree build time: " << Utility::integer_to_string(t.elapsed()) << std::endl;

        fci.build_NCI_density();
        std::cout << "NCI density build time: " << Utility::integer_to_string(t.elapsed()) << std::endl;

        // Save FCI
        {
            std::ofstream file(output_dir + "/FCI.bin", std::ios::binary);
            std::vector<Edge> edges;
            for (const auto& tree : fci.FCF) {
                edges.insert(edges.end(),
                    std::make_move_iterator(tree.begin()),
                    std::make_move_iterator(tree.end()));
            }
            for (const auto& edge : edges) {
                write_edge_binary(file, edge.id, edge.id2, edge.weight.num(), edge.weight.den());
            }
        }

        // Save NCI
        {
            std::ofstream file(output_dir + "/NCI.bin", std::ios::binary);
            for (const auto& edge : fci.NCI) {
                write_edge_binary(file, edge.id, edge.id2, edge.weight.num(), edge.weight.den());
            }
        }

        // Save DNCI
        {
            std::ofstream file(output_dir + "/DNCI.bin", std::ios::binary);
            for (const auto& edge : fci.NCIEdges) {
                write_edge_binary(file, edge.id, edge.id2, edge.weight.num(), edge.weight.den());
            }
        }

    }
    else if (mode == "query") {
        if (argc != 7) {
            std::cerr << "Usage: " << argv[0]
                    << " query <exact|density> <miu> <epsilon_num> <epsilon_denom> <input_file>\n";
            return 1;
        }

        std::string method = argv[2]; // "exact" or "density"
        int miu = std::stoi(argv[3]);
        int epsilon_num = std::stoi(argv[4]);
        int epsilon_denom = std::stoi(argv[5]);
        std::string input_file = argv[6];

        std::ifstream dataset_file(input_file);
        if (!dataset_file) {
            std::cerr << "Error: Cannot open input file " << input_file << "\n";
            return 1;
        }

        int n1, n2, m;
        dataset_file >> n1 >> n2 >> m;
        dataset_file.close();

        FCI fci;
        fci.n1 = n1;

        Fraction epsilon{epsilon_num, epsilon_denom};
        std::string filename = std::filesystem::path(input_file).stem().string();
        std::string output_filename = "result/" + filename + "_" + method + "_" + std::to_string(miu) + "_" + std::to_string(epsilon_num) + "_" + std::to_string(epsilon_denom) + ".txt";
        if (method == "exact") {
            auto result = fci.exactClustering(filename, miu, epsilon);
            write_clusters_to_file(result, output_filename);
        } else if (method == "density") {
            auto result = fci.densityClustering(filename, miu, epsilon);
            write_clusters_to_file(result, output_filename);
        } else {
            std::cerr << "Invalid method: use 'exact' or 'density'\n";
            return 1;
        }


    }
    else {
        std::cerr << "Unknown mode: " << mode << "\n";
        return 1;
    }

    return 0;
}