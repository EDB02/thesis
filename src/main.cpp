#include "algorithm.hpp"
#include "matrix.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 3 && argc != 4)
    {
        cout << "Usage: ./reorder matrix_name block_size [output]" << endl;
        return 1;
    }
    Matrix m;

    m.read_from_file(argv[1]);

    int block_size = stoi(argv[2]);

    cout << "zero blocks before: " << m.count_zero_blocks(block_size) << endl;
    cout << "zero block ratio before: " << m.zero_block_ratio(block_size) << endl;

    
    Matrix m_final = reorder(m, block_size);

    cout << "zero blocks after: " << m_final.count_zero_blocks(block_size) << endl;
    cout << "zero block ratio after: " << m_final.zero_block_ratio(block_size) << endl;

    if(argc == 4)
        m_final.save_permutation(argv[3], block_size);
    else
        m_final.save_permutation("perm.txt", block_size);

    return 0;
}