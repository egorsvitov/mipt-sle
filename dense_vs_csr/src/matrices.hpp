#pragma once
#include <vector>
#include <iostream>
#include <numeric>
#include "vector_operations.hpp"

namespace mtrx {

    template<typename T>
    class dense {
        public:
            dense() : data(), m(0), n(0) {}
            dense(std::size_t m, std::size_t n, const std::vector<T>& data) : m(m), n(n), data(data) {}
            const T& operator()(std::size_t i, std::size_t j) const {
                return data[i * this->width() + j];
            }
            std::vector<T> operator*(const std::vector<T>& x) const {
                std::vector<T> res(m, 0);
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < n; j++) {
                        res[i] += (*this)(i, j)*x[j];
                    }
                }
                return res;
            }
            mtrx::dense<T> operator*(const T& rhs) const {
                std::vector<T> new_data(data.size());
                for (int i = 0; i < data.size(); i++) {
                    new_data[i] = data[i]*rhs;
                }
                mtrx::dense<T> res(m, n, new_data);
                return res;
            }
            mtrx::dense<T> operator-(const mtrx::dense<T>& rhs) const {
                std::vector<T> new_data = data - rhs.get_data();
                return mtrx::dense<T>(m, n, new_data);
            }
            std::size_t width() const {
                return n;
            }
            std::size_t height() const {
                return m;
            }
            void print() {
                for (int i = 0; i < m; i++) {
                    std::cout << "| ";
                    for (int j = 0; j < n; j++) {
                        std::cout << (*this)(i, j) << " ";
                    }
                    std::cout << "|\n";
                }
            }
            const std::vector<T>& get_data() const {
                return data;
            }
        private:
            std::vector<T> data;
            std::size_t m;
            std::size_t n;
    };

    template<typename T>
    mtrx::dense<T> operator*(const T& lhs, const mtrx::dense<T>& rhs) {
        return rhs*lhs;
    }

    template<typename T>
    class csr {
        public:
            csr() : values(), col_indxs(), row_indxs(), m(0), n(0) {}
            csr(std::size_t m, std::size_t n, const std::vector<T>& values, const std::vector<std::size_t>& col_indxs, const std::vector<std::size_t>& row_indxs) : m(m), n(n), values(values), col_indxs(col_indxs), row_indxs(row_indxs) {}
            csr(const mtrx::dense<T>& source) : m(source.height()), n(source.width()), values(), col_indxs(), row_indxs() {
                values.reserve(m*n);
                col_indxs.reserve(m*n);
                row_indxs.reserve(m+1);
                int count = 0;
                row_indxs.push_back(0);
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < n; j++) {
                        if (source(i, j) != 0) {
                            values.push_back(source(i, j));
                            col_indxs.push_back(j);
                            count++;
                        }
                    }
                    row_indxs.push_back(count);
                }
                values.shrink_to_fit();
                col_indxs.shrink_to_fit();
            }
            // T operator()(std::size_t i, std::size_t j) const {
            //     std::size_t row_start = row_indxs[i];
            //     std::size_t row_end = row_indxs[i+1];
            //     std::size_t true_k = 0;
            //     bool s = false;
            //     for (std::size_t k = row_start; k < row_end; k++) {
            //         if (col_indxs[k] == j) {
            //             true_k = k;
            //             s = true;
            //             break;
            //         }
            //     }
            //     return ((s) ? values[true_k] : 0);
            // }
            T operator()(std::size_t i, std::size_t j) const {
                std::size_t row_start = row_indxs[i];
                std::size_t row_end = row_indxs[i+1];
                for (std::size_t k = row_start; k < row_end; k++) {
                    if (col_indxs[k] == j) {
                        return values[k];
                    }
                }
                return 0;
            }
            // std::vector<T> operator*(const std::vector<T>& x) const {
            //     std::vector<T> res(m, 0);
            //     for (int i = 0; i < m; i++) {
            //         for (int j = 0; j < n; j++) {
            //             res[i] += (*this)(i, j)*x[j];
            //         }
            //     }
            //     return res;
            // }
            std::vector<T> operator*(const std::vector<T>& x) const {
                std::vector<T> res(m, 0);
                for (auto i = 0; i < m; i++) {
                    for (auto j = row_indxs[i]; j < row_indxs[i+1]; j++) {
                        res[i] += x[col_indxs[j]]*values[j];
                    }
                }
                return res;
            }
            mtrx::csr<T> operator*(const T& rhs) const {
                std::vector<T> new_values(values.size());
                for (int i = 0; i < values.size(); i++) {
                    new_values[i] = values[i]*rhs;
                }
                mtrx::csr<T> res(m, n, new_values, col_indxs, row_indxs);
                return res;
            }

            std::size_t width() const {
                return n;
            }
            std::size_t height() const {
                return m;
            }
            const std::vector<T>& get_raw_values() const {
                return values;
            }
            const std::vector<std::size_t>& get_raw_cols() const {
                return col_indxs;
            }
            const std::vector<std::size_t>& get_raw_rows() const {
                return row_indxs;
            }
            void print() {
                for (int i = 0; i < m; i++) {
                    std::cout << "| ";
                    for (int j = 0; j < n; j++) {
                        std::cout << (*this)(i, j) << " ";
                    }
                    std::cout << "|\n";
                }
            }
        private:
            std::vector<T> values;
            std::vector<std::size_t> col_indxs;
            std::vector<std::size_t> row_indxs;
            std::size_t m;
            std::size_t n;
    };

template<typename T>
mtrx::csr<T> operator*(const T& lhs, const mtrx::csr<T>& rhs) {
    return rhs*lhs;
}

template<typename T>
mtrx::csr<T> csr_ident(std::size_t n) {
    std::vector<std::size_t> cols(n, 0);
    std::vector<std::size_t> rows(n+1, 0);
    std::iota(cols.begin(), cols.end(), 0);
    std::iota(rows.begin(), rows.end(), 0);
    return mtrx::csr<T>(n, n, std::vector<T>(n, 1), cols, rows);
}

template<typename T>
mtrx::dense<T> dense_ident(std::size_t n) {
    std::vector<T> data(n*n, 0);
    for (int i = 0; i < n; i++) {
        data[i*n + i] = 1;
    }
    return mtrx::dense<T>(n, n, data);
}

}