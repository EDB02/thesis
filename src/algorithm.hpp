#include "matrix.hpp"

#include <vector>
#include <algorithm>
#include <set>
#include <array>
#include <iostream>

using namespace std;

#define n mat.n
#define m mat.m

int find_column(vector<vector<int>> &blocks, vector<vector<int>> &cols, 
    int block_size, int cur_x, int cur_y, int cur_col, int M, int N)
{
    for (int c = 0; c < M - M % block_size; c++)
    {
        if (!blocks[cur_x][c / block_size])
        {
            c += block_size - 1;
            continue;
        }
        if (c / block_size == cur_y || cols[c][cur_x] != 0)
            continue;
        bool ck = 1;
        for (int b = 0; b < N / block_size; b++)
        {
            if (!((blocks[b][c / block_size] != 0 || cols[cur_col][b] == 0) &&
                    (blocks[b][cur_y] != 0 || cols[c][b] == 0)))
            {
                ck = 0;
                break;
            }
        }
        if (ck)
            return c;
    }
    return -1;
}

int swap_col(vector<vector<int>> &blocks, vector<vector<int>> &rows, 
    vector<vector<int>> &cols, set<array<int, 3>> &order, Matrix &ret, 
    int block_size, int cur_y, int cur_col, int M, int N, int swap_c)
{
    int zero_blocks_cnt = 0;
    for (int b = 0; b < N / block_size; b++)
    {
        int before = blocks[b][cur_y];
        blocks[b][cur_y] += cols[swap_c][b] - cols[cur_col][b];

        if (before && blocks[b][cur_y] == 0)
            zero_blocks_cnt++;

        if (order.count({before, b, cur_y}))
        {
            order.erase({before, b, cur_y});
            if (blocks[b][cur_y])
                order.insert({blocks[b][cur_y], b, cur_y});
        }

        before = blocks[b][swap_c / block_size];
        blocks[b][swap_c / block_size] += cols[cur_col][b] - cols[swap_c][b];

        if (before && blocks[b][swap_c / block_size] == 0)
            zero_blocks_cnt++;

        if (order.count({before, b, swap_c / block_size}))
        {
            order.erase({before, b, swap_c / block_size});
            if (blocks[b][swap_c / block_size])
                order.insert({blocks[b][swap_c / block_size], b, swap_c / block_size});
        }
    }
    vector<int> &col1 = ret.get_col(cur_col);
    vector<int> &col2 = ret.get_col(swap_c);
    vector<int> &row_perm = ret.get_row_permutation();
    for (int &r : col1)
    {
        rows[row_perm[r]][cur_y]--;
        rows[row_perm[r]][swap_c / block_size]++;
    }
    for (int &r : col2)
    {
        rows[row_perm[r]][cur_y]++;
        rows[row_perm[r]][swap_c / block_size]--;
    }
    swap(cols[cur_col], cols[swap_c]);
    ret.swap_columns(cur_col, swap_c);
    return zero_blocks_cnt;
}

int find_row(vector<vector<int>> &blocks, vector<vector<int>> &rows, 
    int block_size, int cur_x, int cur_y, int cur_row, int M, int N)
{
    for (int r = 0; r < N - N % block_size; r++)
    {
        if (!blocks[r / block_size][cur_y])
        {
            r += block_size - 1;
            continue;
        }
        if (r / block_size == cur_x || rows[r][cur_y] != 0)
            continue;
        bool ck = 1;
        for (int b = 0; b < M / block_size; b++)
        {
            if (!((blocks[r / block_size][b] != 0 || rows[cur_row][b] == 0) &&
                    (blocks[cur_x][b] != 0 || rows[r][b] == 0)))
            {
                ck = 0;
                break;
            }
        }
        if (ck)
            return r;
    }
    return -1;
}

int swap_row(vector<vector<int>> &blocks, vector<vector<int>> &rows, 
    vector<vector<int>> &cols, set<array<int, 3>> &order, Matrix &ret, 
    int block_size, int cur_x, int cur_row, int M, int N, int swap_r)
{
    int zero_blocks_cnt = 0;
    for (int b = 0; b < M / block_size; b++)
    {
        int before = blocks[cur_x][b];
        blocks[cur_x][b] += rows[swap_r][b] - rows[cur_row][b];

        if (before && blocks[cur_x][b] == 0)
            zero_blocks_cnt++;

        if (order.count({before, cur_x, b}))
        {
            order.erase({before, cur_x, b});
            if (blocks[cur_x][b])
                order.insert({blocks[cur_x][b], cur_x, b});
        }

        before = blocks[swap_r / block_size][b];
        blocks[swap_r / block_size][b] += rows[cur_row][b] - rows[swap_r][b];

        if (before && blocks[swap_r / block_size][b] == 0)
            zero_blocks_cnt++;

        if (order.count({before, swap_r / block_size, b}))
        {
            order.erase({before, swap_r / block_size, b});
            if (blocks[swap_r / block_size][b])
                order.insert({blocks[swap_r / block_size][b], swap_r / block_size, b});
        }
    }
    vector<int> &row1 = ret.get_row(cur_row);
    vector<int> &row2 = ret.get_row(swap_r);
    vector<int> &col_perm = ret.get_col_permutation();
    for (int &c : row1)
    {
        cols[col_perm[c]][cur_x]--;
        cols[col_perm[c]][swap_r / block_size]++;
    }
    for (int &c : row2)
    {
        cols[col_perm[c]][cur_x]++;
        cols[col_perm[c]][swap_r / block_size]--;
    }
    swap(rows[cur_row], rows[swap_r]);
    ret.swap_rows(cur_row, swap_r);
    return zero_blocks_cnt;
}

Matrix reorder(Matrix &mat, int block_size)
{
    Matrix ret = mat.copy();
    vector<vector<int>> blocks(n / block_size, vector<int>(m / block_size, 0));
    vector<vector<int>> rows(n, vector<int>(m / block_size, 0));
    vector<vector<int>> cols(m, vector<int>(n / block_size, 0));

    int zero_blocks_cnt = (n / block_size) * (m / block_size);

    for (int i = 0; i < n - n % block_size; i++)
    {
        vector<int> &row = ret.get_row(i);
        vector<int> &col_perm = ret.get_col_permutation();
        for (int &j : row)
        {
            if (col_perm[j] < m - m % block_size)
            {
                if (blocks[i / block_size][col_perm[j] / block_size] == 0)
                    zero_blocks_cnt--;
                blocks[i / block_size][col_perm[j] / block_size]++;
                rows[i][col_perm[j] / block_size]++;
                cols[col_perm[j]][i / block_size]++;
            }
        }
    }

    set<array<int, 3>> order;

    for (int i = 0; i < n / block_size; i++)
    {
        for (int j = 0; j < m / block_size; j++)
        {
            if (blocks[i][j] != 0)
                order.insert({blocks[i][j], i, j});
        }
    }

    int prev_block_cnt = -1;
    int cnt = 0;
    int tot = 0;

    while (order.size())
    {
        if (prev_block_cnt != zero_blocks_cnt)
        {
            prev_block_cnt = zero_blocks_cnt;
            cnt = 0;
        }
        if (++cnt > n)      //soft limit: if the number of zero blocks does not change for n iteration it stops
            break;
        if (++tot > 4 * n)  //hard limit: the number of iterations is limited to k * n
        {
            cout << "hard limit" << endl;
            break;
        }
        auto top = *order.begin();
        int cur_x = top[1];
        int cur_y = top[2];

        order.erase(order.begin());

        if (blocks[cur_x][cur_y] == 0)
            continue;

        for (int cur = 0; cur < block_size; cur++)
        {
            int cur_row = cur + cur_x * block_size;
            if (rows[cur_row][cur_y] != 0)
            {
                int swap_r = find_row(blocks, rows, block_size, cur_x, cur_y, cur_row, m, n);
                
                if (swap_r != -1)
                {
                    zero_blocks_cnt += swap_row(blocks, rows, cols, order, ret, 
                                block_size, cur_x, cur_row, m, n, swap_r);
                }
            }

            int cur_col = cur + cur_y * block_size;
            if (cols[cur_col][cur_x] != 0)
            {
                int swap_c = find_column(blocks, cols, block_size, cur_x, cur_y, cur_col, m, n);
                
                if (swap_c != -1)
                {
                    zero_blocks_cnt += swap_col(blocks, rows, cols, order, ret, 
                                block_size, cur_y, cur_col, m, n, swap_c);
                }
            }
        }
    }
    return ret;
}


#undef n
#undef m