#include <iostream>
#include <fstream>
#include <filesystem>  // C++17
#include "FCI.h"
#include "Timer.h"
#include "utility.hpp"
#include <filesystem>

void write_clusters_to_file(const std::vector<Cluster>& clusters, const std::string& file_path) {
    std::filesystem::path path(file_path);
    if (!path.parent_path().empty() && !std::filesystem::exists(path.parent_path())) {
        std::filesystem::create_directories(path.parent_path());
        std::cout << "Created directory: " << path.parent_path() << std::endl;
    }

    std::ofstream fout(file_path);
    if (!fout) {
        std::cerr << "Failed to open file: " << file_path << std::endl;
        return;
    }

    if (clusters.empty()) {
        std::cout << "Warning: clusters vector is empty. File will be empty." << std::endl;
    }

    for (const auto& cluster : clusters) {
        fout << "Cluster rep: " << cluster.rep << "\n";

        fout << "  Cores: ";
        if (cluster.cores.empty()) fout << "(empty)";
        for (int core : cluster.cores) {
            fout << core << " ";
        }

        fout << "\n  Non-Cores: ";
        if (cluster.nonCores.empty()) fout << "(empty)";
        for (int nc : cluster.nonCores) {
            fout << nc << " ";
        }

        fout << "\n\n";
    }

    fout.flush();
    std::cout << "Successfully wrote " << clusters.size() << " clusters to " << file_path << std::endl;
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
        fci.mu = miu;
        fci.read_vectors(input_file);
        std::cout << "finished reading graph" << std::endl;
        Timer t;
        fci.bulid_core_index();
        std::cout << "Core index build time: " << Utility::integer_to_string(t.elapsed()) << std::endl;

        fci.build_spanning_tree();
        std::cout << "Spanning tree build time: " << Utility::integer_to_string(t.elapsed()) << std::endl;

        fci.build_NCI_density();
        std::cout << "NCI density build time: " << Utility::integer_to_string(t.elapsed()) << std::endl;

        // Save FCI
        {
            std::vector<Edge> all_edges;
            for (const auto& tree : fci.FCF) {
                all_edges.insert(all_edges.end(), tree.begin(), tree.end());
            }

            std::sort(all_edges.begin(), all_edges.end(), [](const Edge &a, const Edge &b) {
                return a.weight < b.weight;
            });

            std::ofstream file(output_dir + "/newFCI.txt");

            for (const auto& edge : all_edges) {
                file << edge.id << " " << edge.id2 << " " << edge.weight << "\n";
            }
            file.close();
        }
    
        {
            std::ofstream file(output_dir + "/newNCI.txt");
                for (const auto& edge : fci.NCI) {
                    int u = edge.id, v = edge.id2;
                    file << u << " " << v << " " << edge.weight << "\n";
                }
        }
    
        {
            std::ofstream file(output_dir + "/newDNCI.txt");
            for (const auto& edge : fci.NCIEdges) {
                int u = edge.id, v = edge.id2;
                file << u << " " << v << " " << edge.weight << "\n";
            }
        }

    }
    else if (mode == "query") {
        if (argc != 6) {
            std::cerr << "Usage: " << argv[0]
                    << " query <exact|density> <miu> <epsilon> <input_file>\n";
            return 1;
        }

        std::string method = argv[2]; // "exact" or "density"
        int miu = std::stoi(argv[3]);
        double epsilon = std::stod(argv[4]);
        std::string input_file = argv[5];

        std::ifstream dataset_file(input_file);
        if (!dataset_file) {
            std::cerr << "Error: Cannot open input file " << input_file << "\n";
            return 1;
        }

        FCI fci;
        fci.mu = miu;

        std::ifstream fin(input_file);
        if (!fin.is_open()) {
            std::cerr << "Failed to open file: " << input_file << std::endl;
            return 1;
        }
    
        std::string line;
    
        if (std::getline(fin, line)) {
            std::stringstream ss(line);
            ss >> fci.n1;
        }
        
        std::string filename = std::filesystem::path(input_file).stem().string();
        std::string output_filename = "result/" + filename + "_" + method + "_" + std::to_string(miu) + "_" + std::to_string(epsilon) + ".txt";
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