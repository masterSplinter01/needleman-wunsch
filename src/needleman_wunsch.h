#include <iostream>
#include <string>
#include <map>

class dna_comp
{
private:
    std::string _align1, _align2;
    std::string _dna1, _dna2;
    
    int _gap_penalty;
    int _mismatch_penalty;
    int _match_penalty;
    
    size_t _matrix_rows, _matrix_cols;
    
    int** _similarity_matrix;
public:
    dna_comp(std::string dna1, std::string dna2);
    ~dna_comp();
    
    void init_matrix();
    void print_matrix();
    int sim(size_t i, size_t j);
    void trace_back();
    void print_strings();
};
int max(int a, int b, int c);