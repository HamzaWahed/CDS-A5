#include <cstdlib>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/sd_vector.hpp>

#include <string>
#include <vector>

#define CHAR_COUNT 256

int n = 0;
int r = 0;
int sigma = 0;

// all the characters in T or BWT(T)
std::vector<unsigned char> alphabet;

// to map a character to the index
// Ex. #:0, $:1, A:2, C:3, G:4, T:5
std::vector<int> char_to_index;

// The data structures needed for the assignment
sdsl::sd_vector B_F_sparse;
sdsl::sd_vector B_L_sparse;

std::vector<int> C;
std::vector<char> H_L;
std::vector<std::unique_ptr<sdsl::bit_vector>> B_x;

// Rank and Select data structures on B_F and B_L
sdsl::rank_support_sd<> rank_B_L;
sdsl::select_support_sd<> select_B_L;
sdsl::rank_support_sd<> rank_B_F;
sdsl::select_support_sd<> select_B_F;

// Rank data structure for B_x bit vectors
std::vector<std::unique_ptr<sdsl::rank_support_v<>>> B_x_ranks;


void deserialize_data(char *inputFileName) {
    std::ifstream in_file(inputFileName, std::ios::in | std::ios::binary);

    in_file.read(reinterpret_cast<char *>(&n), sizeof(n));
    in_file.read(reinterpret_cast<char *>(&r), sizeof(r));
    in_file.read(reinterpret_cast<char *>(&sigma), sizeof(sigma));
    std::cerr << "n: " << n << "\n";
    std::cerr << "r: " << r << "\n";
    std::cerr << "sigma: " << sigma << "\n";

    H_L.resize(r);
    in_file.read(reinterpret_cast<char *>(&H_L[0]), r * sizeof(H_L[0]));

    C.resize(sigma);
    in_file.read(reinterpret_cast<char *>(&C[0]), sigma * sizeof(C[0]));

    alphabet.resize(sigma);
    in_file.read(reinterpret_cast<char *>(&alphabet[0]),
                 sigma * sizeof(alphabet[0]));
    std::cerr << "\nThe alphabet in the BWT:\n";
    for (int i = 0; i < sigma; i++)
        std::cerr << "\t" << i << " -> " << static_cast<int>(alphabet[i]) << "("
                  << alphabet[i] << ")\n";
    std::cerr << "\n";
    char_to_index.resize(CHAR_COUNT);
    in_file.read(reinterpret_cast<char *>(&char_to_index[0]),
                 CHAR_COUNT * sizeof(char_to_index[0]));

    B_L_sparse.load(in_file);

    // Building the rank and select data structures for querying B_L
    rank_B_L = sdsl::rank_support_sd<>(
        &B_L_sparse); // usage example: rank_B_L(i) gives the rank result at
                      // index i on B_L
    select_B_L = sdsl::select_support_sd<>(
        &B_L_sparse); // usage example: select_B_L(i) gives the select result at
                      // index i on B_L

    B_F_sparse.load(in_file);

    // Building the rank and select data structures for querying B_L
    rank_B_F = sdsl::rank_support_sd<>(
        &B_F_sparse); // usage example: rank_B_L(i) gives the rank result at
                      // index i on B_L
    select_B_F = sdsl::select_support_sd<>(
        &B_F_sparse); // usage example: select_B_L(i) gives the select result at
                      // index i on B_L

    for (int i = 0; i < sigma; i++) {
        auto new_b_vector = std::make_unique<sdsl::bit_vector>();
        new_b_vector->load(in_file);
        B_x.push_back(std::move(new_b_vector));
    }

    // create the rank objects for the B_x bit vectors
    for (auto &B : B_x) {
        B_x_ranks.emplace_back(std::unique_ptr<sdsl::rank_support_v<>>(
            new sdsl::rank_support_v<>(B.get())));
    }
    // Example: code to perform rank query on B_2 at position 10:
    // std::cerr << (*B_x_ranks[2])(10) << "\n";

    in_file.close();
}

int main(int argc, char **argv) {
    std::cerr << "usage: load DATASET.ri\n";
    std::cerr << "The DATASET.ri is generated by the build script.\n\n";

    deserialize_data(argv[1]);
    std::cerr << "All the arrays and bit vectors are loaded.\n\n";

    size_t idx, run_idx, offset_idx, pred_run_head, run_head_idx_in_F;
    char run_head;

    std::cerr << "Enter position in L: ";
    while (std::cin >> (idx)) {
        run_idx = rank_B_L(idx + 1) - 1;
        offset_idx = idx - select_B_L(run_idx + 1);

        run_head = H_L[run_idx];
        pred_run_head = (*B_x_ranks[char_to_index[run_head]])(run_idx);

        run_head_idx_in_F = C[char_to_index[run_head]] + pred_run_head;

        // TODO: C array definition does not match the definition in the
        // assignment specs
        std::cout << "Position in F: "
                  << select_B_F(rank_B_F(C[char_to_index[run_head]] + 1) +
                                pred_run_head) +
                         offset_idx
                  << std::endl;
        std::cerr << "Enter position in L: ";
    }
    return 0;
}
