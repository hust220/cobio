#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>

struct Alignment {
    int gap0 = 10; // opening gap
    int gap = 1; // extension gap
    std::vector<std::vector<int>> sim {
        { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
        {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
        {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
        {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
        { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
        {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
        {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
        { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
        {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
        {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
        {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
        {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
        {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
        {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
        {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
        { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
        { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
        {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
        {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
        { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
        {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
        {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
        { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
        {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
    };
    std::vector<char> aa {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'};
    std::map<char, int> ind;

    Alignment() {
        for (int i = 0; i < int(aa.size()); i++) {
            ind[aa[i]] = i;
        }
    }

    std::vector<int> format_seq(const std::string &seq) {
        int n = seq.size();
        std::vector<int> v(n);
        for (int i = 0; i < n; i++) {
            if (ind.count(seq[i]) == 0) {
                v[i] = ind['*'];
            } else {
                v[i] = ind[seq[i]];
            }
        }
        return v;
    }

    /* 
     * Returns the Needleman-Wunsch score for the best alignment of a and b
     * and stores the aligned sequences in a_aligned and b_aligned
     */
    std::tuple<int, std::string, std::string> operator()(const std::string &a_, const std::string &b_) {
        auto &&a = format_seq(a_);
        auto &&b = format_seq(b_);

        int n = a.size();
        int m = b.size();

        std::vector<std::vector<int>> A(n+1, std::vector<int>(m+1, 0));
        std::vector<std::vector<int>> B(n+1, std::vector<int>(m+1, 0)); // 0: no gap, 1: gap

        for (int i = 0; i <= n; ++i) A[i][0] = (i == 0) ? 0 : (-gap0 - gap*(i-1));
        for (int i = 0; i <= m; ++i) A[0][i] = (i == 0) ? 0 : (-gap0 - gap*(i-1));

        int s1, s2, s3;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= m; ++j) {
                s1 = A[i-1][j-1] + sim[a[i-1]][b[j-1]];
                s2 = A[i-1][j] - (B[i-1][j] == 0 ? gap0 : gap);
                s3 = A[i][j-1] - (B[i][j-1] == 0 ? gap0 : gap);
                if (s1 >= s2 && s1 >= s3) {
                    A[i][j] = s1;
                    B[i][j] = 0;
                } else if (s2 >= s1 && s2 >= s3) {
                    A[i][j] = s2;
                    B[i][j] = 1;
                } else if (s3 >= s1 && s3 >= s2) {
                    A[i][j] = s3;
                    B[i][j] = 2;
                }
            }
        }

//        for (int i = 0; i <= n; i++) {
//            for (int j = 0; j <= m; j++) {
//                printf("%5d ", A[i][j]);
//            }
//            std::cout << std::endl;
//        }

        std::list<int> a_aligned;
        std::list<int> b_aligned;
        int i = n;
        int j = m;
        while (i > 0 && j > 0) {
            if (B[i][j] == 0) {
                a_aligned.push_front(--i);
                b_aligned.push_front(--j);
            } else if (B[i][j] == 1) {
                a_aligned.push_front(--i);
                b_aligned.push_front(-1);
            } else {
                a_aligned.push_front(-1);
                b_aligned.push_front(--j);
            }
        }
        
        while (i > 0) {
            a_aligned.push_front(--i);
            b_aligned.push_front(-1);
        }
        
        while (j > 0) {
            a_aligned.push_front(-1);
            b_aligned.push_front(--j);
        }

        int sz = a_aligned.size();
        std::string a1(sz, '-');
        std::string b1(sz, '-');
        std::transform(a_aligned.begin(), a_aligned.end(), a1.begin(), [&a_](int i){
            return (i == -1) ? '-' : a_[i];
        });
        std::transform(b_aligned.begin(), b_aligned.end(), b1.begin(), [&b_](int i){
            return (i == -1) ? '-' : b_[i];
        });
       
        return std::make_tuple(A[n][m], a1, b1);
    }
};



