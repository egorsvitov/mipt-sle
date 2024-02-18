#include <vector>

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
            T width() const {
                return this->n;
            }
            T height() const {
                return this->m;
            }
        private:
            std::vector<T> data;
            std::size_t m;
            std::size_t n;
    };



    template<typename T>
    class csr {
        public:
            csr() : values(), col_indxs(), row_indxs(), m(0), n(0) {}
            csr(std::size_t m, std::size_t n, const std::vector<T>& values, const std::vector<T>& col_indxs, const  std::vector<T>& row_indxs) : m(m), n(n), values(values), col_indxs(col_indxs), row_indxs(row_indxs) {}
            csr(const mtrx::dense<T>& source) : m(source.height()), n(source.width()), values(), col_indxs(), row_indxs() {
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
            }
            const T& operator()(std::size_t i, std::size_t j) const {
                unsigned int row_start = row_indxs[i];
                unsigned int row_end = row_indxs[i+1];
                unsigned int true_k = 0;
                for (unsigned int k = row_start; k < row_end; k++) {
                    if (col_indxs[k] == j) {
                        true_k = k;
                        break;
                    }
                }
                return values[true_k];
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
            T width() const {
                return this->n;
            }
            T height() const {
                return this->m;
            }
            std::vector<T> get_raw_values() const {
                return values;
            }
            std::vector<unsigned int> get_raw_cols() const {
                return col_indxs;
            }
            std::vector<unsigned int> get_raw_rows() const {
                return row_indxs;
            }
        private:
            std::vector<T> values;
            std::vector<unsigned int> col_indxs;
            std::vector<unsigned int> row_indxs;
            std::size_t m;
            std::size_t n;
    };
}

template<typename T>
std::vector<T> vec_add(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    std::vector<T> res(rhs.size());
    for (int i = 0; i < rhs.size(); i++) {
        res[i] = lhs[i] + rhs[i];       
    }
    return res;
}

template<typename T>
std::vector<T> vec_scal_prod(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    std::vector<T> res(rhs.size());
    for (int i = 0; i < rhs.size(); i++) {
        res[i] = lhs[i] * rhs[i];       
    }
    return res;
}