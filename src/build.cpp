#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/sd_vector.hpp>
#include <string>
#include <vector>

#define CHAR_COUNT 256

std::vector<int> chars(CHAR_COUNT, 0);

char *BWT;
int n = 0;
int r = 0;
int sigma = 0;

// to hold the BWT offset of every run head
std::vector<int> run_heads;

// for computing the LF mapping for building the B_F vector
std::vector<std::unique_ptr<sdsl::bit_vector>> occs;
std::vector<std::unique_ptr<sdsl::rank_support_v<>>> occs_rank;

// to keep the count of each character, used in computing C
std::vector<int> counts;

// all the characters in T or BWT(T)
std::vector<unsigned char> alphabet;

// to map a character to the index
// Ex. #:0, $:1, A:2, C:3, G:4, T:5
std::vector<int> char_to_index;

// The data structures needed for the assignment
sdsl::bit_vector B_F;
sdsl::bit_vector B_L;

sdsl::sd_vector B_F_sparse;
sdsl::sd_vector B_L_sparse;

std::vector<int> C;
std::vector<char> H_L;
std::vector<std::unique_ptr<sdsl::bit_vector>> B_x;

void build_BL(char *BWT, int n) {
    char last_char = ' ';
    int run_length = 0;

    for (int i = 0; i < n; i++) {

        char current_char = BWT[i];

        // increase count of the characters to find the characters that exist in
        // the BWT
        chars[current_char] += 1;

        if (last_char == ' ' or current_char != last_char) {
            // a new run is starting
            r += 1;

            // save the run head character in H_L
            H_L.push_back(current_char);

            // save the run head offset
            run_heads.push_back(i);

            // set the run head in B_L
            B_L[i] = 1;

            run_length = 1;

        } else {

            // a new row in the run
            run_length += 1;
        }

        last_char = current_char;
    }

    char_to_index.resize(CHAR_COUNT);
    std::fill(char_to_index.begin(), char_to_index.end(), CHAR_COUNT);

    sigma = 0;
    std::cerr << "The count of each character in the BWT:\n";
    for (int i = 0; i < CHAR_COUNT; i++) {
        if (chars[i] != 0) {
            std::cerr << "\t" << static_cast<char>(i) << ": " << chars[i]
                      << "\n";
            counts.push_back(chars[i]);
            alphabet.push_back(static_cast<char>(i));
            char_to_index[i] = sigma;
            sigma += 1;
            sdsl::bit_vector *new_b_vector = new sdsl::bit_vector(r, 0);
            B_x.emplace_back(std::unique_ptr<sdsl::bit_vector>(new_b_vector));

            sdsl::bit_vector *new_occ_vector = new sdsl::bit_vector(n, 0);
            occs.emplace_back(
                std::unique_ptr<sdsl::bit_vector>(new_occ_vector));
        }
    }
    std::cerr << "\n";
}

void build_C() {
    int count = 0;

    for (int i = 0; i < counts.size(); i++) {

        C.push_back(count);
        count += counts[i];
    }
}

void fill_occs_and_bs(int n, int r) {
    for (int i = 0; i < r; i++) {
        (*B_x[char_to_index[static_cast<int>(H_L[i])]])[i] = 1;
    }

    for (int i = 0; i < n; i++) {
        (*occs[char_to_index[static_cast<int>(BWT[i])]])[i] = 1;
    }
}

int LF(int bwt_row) {
    int lf = 0;

    int alphabet_index = char_to_index[static_cast<int>(BWT[bwt_row])];

    // add up the counts of the smaller characters
    lf += C[alphabet_index];

    // add the rank of the bwt_row
    auto &occ_rank = *occs_rank[alphabet_index];
    lf += static_cast<int>(occ_rank(bwt_row));

    return lf;
}

std::string reconstruct() {
    std::string original_T = "";

    int i = 0;
    do {
        original_T = BWT[i] + original_T;
        i = LF(i);
    } while (BWT[i] != alphabet[0]);

    return original_T;
}

void loadBWT(char *textFileName) {
    FILE *textFile = fopen(textFileName, "r");
    fseek(textFile, 0, SEEK_END);
    n = (int)ftell(textFile);
    rewind(textFile);

    BWT = (char *)malloc(n);
    fread(BWT, 1, n, textFile);
    fclose(textFile);
    return;
}

void serialize_data(char *outputFileName) {
    std::ofstream out_file(outputFileName, std::ios::out | std::ios::binary);

    out_file.write(reinterpret_cast<char *>(&n), sizeof(n));
    out_file.write(reinterpret_cast<char *>(&r), sizeof(r));
    out_file.write(reinterpret_cast<char *>(&sigma), sizeof(sigma));

    out_file.write(reinterpret_cast<char *>(&H_L[0]), r * sizeof(H_L[0]));

    out_file.write(reinterpret_cast<char *>(&C[0]), sigma * sizeof(C[0]));

    out_file.write(reinterpret_cast<char *>(&alphabet[0]),
                   sigma * sizeof(alphabet[0]));
    out_file.write(reinterpret_cast<char *>(&char_to_index[0]),
                   CHAR_COUNT * sizeof(char_to_index[0]));

    B_L_sparse.serialize(out_file);
    B_F_sparse.serialize(out_file);

    for (int i = 0; i < sigma; i++) {
        auto B_x_i = *B_x[i];
        B_x_i.serialize(out_file);
    }

    out_file.close();
}

int main(int argc, char **argv) {
    std::cerr << "usage: build DATASET.bwt DATASET.ri\n";
    std::cerr << "The DATASET.ri is the output of this program.\n\n";

    loadBWT(argv[1]);

    // initialize the B_L and B_F bit vectors
    B_L = sdsl::bit_vector(n, 0);
    B_F = sdsl::bit_vector(n, 0);

    // building the B_L bit vector and the H_L vector
    // and some other auxiliary dat structure
    build_BL(BWT, n);
    std::cerr << "B_L and H_L are created.\n";
    std::cerr << "n: " << n << "\n";
    std::cerr << "r: " << r << "\n";

    // build the B_x bit vectors and occs bit vectors
    fill_occs_and_bs(n, r);
    std::cerr << "B_x bit vectors are created.\n";

    // create the rank objects for the occurance bit vectors
    for (auto &occ : occs) {
        occs_rank.emplace_back(std::unique_ptr<sdsl::rank_support_v<>>(
            new sdsl::rank_support_v<>(occ.get())));
    }

    sdsl::rank_support_v<> rank_B_L = sdsl::rank_support_v<>(&B_L);

    build_C();
    std::cerr << "C array is created.\n";

    for (int i = 0; i < r; i++) {
        int lf_i = LF(run_heads[i]);
        B_F[lf_i] = 1;
    }
    std::cerr << "B_F is created.\n";

    // Create the sparse bit vectors for B_F and B_L
    B_F_sparse = sdsl::sd_vector<>(B_F);
    B_L_sparse = sdsl::sd_vector<>(B_L);
    std::cerr << "Sparse bit vectors are created.\n";

    // For inspecting how the vectors look on a small example
    /* std::cerr << reconstruct() << "\n";
    std::cerr << "B_L: " << B_L_sparse << "\n";
    std::cerr << "B_F: " << B_F_sparse << "\n";
    std::cerr << "H_L: ";
    for (int i = 0; i < r; i++)
        std::cerr << H_L[i];
    std::cerr << "\n";
    for (int i = 0; i < sigma; i++) {
        std::cerr << "B_" << i << ": " << *B_x[i] << "\n";
    } */

    serialize_data(argv[2]);
    std::cerr << "\nThe output is stored at " << argv[2] << "\n\n";
}
