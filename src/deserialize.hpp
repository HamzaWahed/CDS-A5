#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/bit_vectors.hpp>

#include <vector>
#include <string>

#define CHAR_COUNT 256

void deserialize_data(char *inputFileName);

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
std::vector<std::unique_ptr<sdsl::bit_vector> > B_x;

// Rank and Select data structures on B_F and B_L
sdsl::rank_support_sd<> rank_B_L;
sdsl::select_support_sd<> select_B_L;
sdsl::rank_support_sd<> rank_B_F;
sdsl::select_support_sd<> select_B_F;

// Rank data structure for B_x bit vectors
std::vector<std::unique_ptr<sdsl::rank_support_v< > > > B_x_ranks;