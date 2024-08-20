#pragma once

#include <string>
#include <vector>
#include <numeric>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <climits>

using namespace std;

class Matrix
{
public:
    int64_t n, m;
    Matrix() {};
    Matrix(int n, int m)
    {
        assert(n > 0 && m > 0);
        this->n = n;
        this->m = m;
        resize();
    }
    Matrix copy()
    {
        Matrix ret(this->n, this->m);
        ret.col = this->col;
        ret.row = this->row;
        ret.row_index_inverse = this->row_index_inverse;
        ret.col_index_inverse = this->col_index_inverse;
        ret.row_index = this->row_index;
        ret.col_index = this->col_index;
        return ret;
    }
    void swap_rows(int i, int j)
    {
        assert(i >= 0 && i < n && j >= 0 && j < n);
        row_index_inverse[row_index[i]] = j;
        row_index_inverse[row_index[j]] = i;
        swap(row_index[i], row_index[j]);
    }
    void swap_columns(int i, int j)
    {
        assert(i >= 0 && i < m && j >= 0 && j < m);
        col_index_inverse[col_index[i]] = j;
        col_index_inverse[col_index[j]] = i;
        swap(col_index[i], col_index[j]);
    }
    bool read_from_file(string path)
    {
        ifstream f(path);
        if(f.fail()) return 1;
        string line;
        while (getline(f, line))
        {
            if (line[0] != '%')
                break;
        }

        int cnt, a, b;
        std::istringstream iss(line);
        iss >> n >> m >> cnt;
        resize();

        while (cnt--)
        {
            f >> a >> b;
            f.ignore(UINT_MAX, '\n');
            row[a - 1].push_back(b - 1);
            col[b - 1].push_back(a - 1);
        }
        for (int i = 0; i < n; i++)
        {
            sort(row[i].begin(), row[i].end());
        }
        f.close();
        return 0;
    }
    int64_t count_zero_blocks(int64_t block_size)
    {
        int64_t cnt = this->total_blocks(block_size);
        vector<bool> is_zero((m + block_size - 1) / block_size);
        for (int i = 0; i < n; i++)
        {
            if (i % block_size == 0)
                fill(is_zero.begin(), is_zero.end(), 1);
            for (int &j : row[row_index[i]])
            {
                if (is_zero[col_index_inverse[j] / block_size])
                {
                    cnt--;
                    is_zero[col_index_inverse[j] / block_size] = 0;
                }
            }
        }
        return cnt;
    }
    int64_t total_blocks(int64_t block_size)
    {
        return ((n + block_size - 1) / block_size) * (int64_t)((m + block_size - 1) / block_size);
    }
    double zero_block_ratio(int64_t block_size)
    {
        return count_zero_blocks(block_size) / (double)total_blocks(block_size);
    }
    void reset_index_order()
    {
        iota(row_index.begin(), row_index.end(), 0);
        iota(col_index.begin(), col_index.end(), 0);
        iota(row_index_inverse.begin(), row_index_inverse.end(), 0);
        iota(col_index_inverse.begin(), col_index_inverse.end(), 0);
    }
    bool compare(Matrix &mat)
    {
        if (mat.n != this->n || mat.m != this->m)
            return 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < static_cast<int>(row[i].size()); j++)
            {
                if (mat.row[i][j] != this->row[i][j])
                    return 0;
            }
        }
        return 1;
    }
    void set_row_permutation(vector<int> permutation)
    {
        this->row_index = permutation;
        int ind = 0;
        for (int &i : permutation)
        {
            this->row_index_inverse[i] = ind++;
        }
    }
    void set_column_permutation(vector<int> permutation)
    {
        this->col_index = permutation;
        int ind = 0;
        for (int &i : permutation)
        {
            this->col_index_inverse[i] = ind++;
        }
    }
    vector<int> &get_row_permutation()
    {
        return this->row_index_inverse;
    }
    vector<int> &get_col_permutation()
    {
        return this->col_index_inverse;
    }
    vector<int> &get_row(int i)
    {
        return row[row_index[i]];
    }
    vector<int> &get_col(int i)
    {
        return col[col_index[i]];
    }
    void save_permutation(string path, int block_size)
    {
        ofstream out(path, ios_base::app);
        out << "Block size: " << block_size << "\n";
        for (int i : row_index)
            out << i << " ";
        out << "\n";
        for (int i : col_index)
            out << i << " ";
        out << endl;
        out.close();
    }
    void save_blocks_density(string path, int block_size)
    {
        ofstream out(path);
        vector<int> cnt((m + block_size - 1) / block_size);
        for (int i = 0; i < n; i++)
        {
            if (i % block_size == 0)
            {
                if (i != 0)
                {
                    for (int j = 0; j < m / block_size; j++)
                    {
                        if (cnt[j])
                            out << i / block_size << ";" << j << ";" << cnt[j] << ";\n";
                    }
                }
                fill(cnt.begin(), cnt.end(), 0);
            }
            for (int &j : row[row_index[i]])
            {
                cnt[col_index_inverse[j] / block_size]++;
            }
        }
    }
    void generate(double density)
    {
        int64_t cnt = (n * m) * density;
        set<pair<int, int>> s;
        while (cnt)
        {
            int x = rand() % n;
            int y = rand() % m;
            if (s.count({x, y}))
                continue;
            s.insert({x, y});
            row[x].push_back(y);
            col[y].push_back(x);
            cnt--;
        }
    }
    int64_t nonzero_count()
    {
        int64_t nonzero = 0;
        for (int i = 0; i < n; i++)
            nonzero += row[i].size();
        return nonzero;
    }
    double nonzero_density(int64_t block_size)
    {
        return nonzero_count() / (double)((n * m - count_zero_blocks(block_size)) * (block_size * block_size));
    }

private:
    void resize()
    {
        row.resize(n);
        col.resize(m);
        row_index.resize(n);
        col_index.resize(m);
        row_index_inverse.resize(n);
        col_index_inverse.resize(m);
        reset_index_order();
    }

    vector<vector<int>> row, col;

    vector<int> row_index, col_index;
    vector<int> row_index_inverse, col_index_inverse;
};