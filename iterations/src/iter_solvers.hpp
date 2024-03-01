#include "../../dense_vs_csr/src/vector_operations.hpp"
#include "../../dense_vs_csr/src/matrices.hpp"
#include <cmath>

namespace slv {

    template<typename T, typename M>
    std::vector<T> simple_solve(const M& A, const std::vector<T>& b, const std::vector<T> x0, T tolerance, T tau, int N_max) {
        std::vector<T> x = x0;
        T delta = tolerance + 1;
        for (int cycle = 0; (delta > tolerance) && (cycle < N_max); cycle++) {
            delta = 0;
            auto delta_x = tau*(A*x - b);
            for (int i = 0; i < delta_x.size(); i++) {
                delta += delta_x[i]*delta_x[i];
            }
            delta = std::sqrt(delta);
            x = x - delta_x;
            //std::cout << cycle << " " << delta << std::endl;
        }
        return x;
    }

    template<typename T, typename M>
    std::vector<T> simple_solve(const M& A, const std::vector<T>& b, const std::vector<T> x0, int N, T tau) {
        std::vector<T> x = x0;
        for (int cycle = 0; cycle < N; cycle++) {
            T delta = 0;
            auto delta_x = tau*(A*x - b);
            for (int i = 0; i < delta_x.size(); i++) {
                delta += delta_x[i]*delta_x[i];
            }
            x = x - delta_x;
        }
        return x;
    }

    template<typename T>
    std::vector<T> jacobi_solve(const mtrx::csr<T>& A, const std::vector<T>& b, const std::vector<T> x0, T tolerance, int N_max) {
        std::vector<T> x = x0;
        T delta = tolerance + 1;
        for (int cycle = 0; (delta > tolerance) && (cycle < N_max); cycle++) {
            std::vector<T> old_x = x;
            for (auto i = 0u; i < x.size(); i++) {
                T xi = 0;
                for (auto j = A.get_raw_rows()[i]; j < A.get_raw_rows()[i+1]; j++) {
                    if (i != A.get_raw_cols()[j]) {
                        xi += old_x[A.get_raw_cols()[j]]*A.get_raw_values()[j];
                    }
                }
                x[i] = xi;
            }
            x = b - x;
            delta = 0;
            for (auto i = 0; i < x.size(); i++) {
                x[i] /= A(i, i);
                delta += (old_x[i] - x[i])*(old_x[i] - x[i]);
            }
            delta = std::sqrt(delta);
        }
        return x;
    }

    template<typename T>
    std::vector<T> jacobi_solve(const mtrx::csr<T>& A, const std::vector<T>& b, const std::vector<T> x0, int N) {
        std::vector<T> x = x0;
        for (int cycle = 0; cycle < N; cycle++) {
            std::vector<T> old_x = x;
            for (auto i = 0u; i < x.size(); i++) {
                T xi = 0;
                for (auto j = A.get_raw_rows()[i]; j < A.get_raw_rows()[i+1]; j++) {
                    if (i != A.get_raw_cols()[j]) {
                        xi += old_x[A.get_raw_cols()[j]]*A.get_raw_values()[j];
                    }
                }
                x[i] = xi;
            }
            x = b - x;
            for (auto i = 0; i < x.size(); i++) {
                x[i] /= A(i, i);
            }
        }
        return x;
    }

    template<typename T>
    std::vector<T> gauss_seidel_solve(const mtrx::dense<T>& A, const std::vector<T>& b, const std::vector<T> x0, T tolerance, int N_max) {
        std::vector<T> x = x0;
        T delta = tolerance + 1;
        for (int cycle = 0; (delta > tolerance) && (cycle < N_max); cycle++) {
            std::vector<T> old_x = x;
            for (auto i = 0u; i < x.size(); i++) {
                T xi = 0;
                for (auto j = i + 1; j < A.width(); j++) {
                    xi += A(i, j)*x[j];
                }
                for (auto j = 0; j < i; j++) {
                    xi += A(i, j)*x[j];
                }
                x[i] = (b[i] - xi)/A(i, i);
            }
            delta = 0;
            for (auto i = 0u; i < x.size(); i++) {
                delta += (x[i]-old_x[i])*(x[i] - old_x[i]);
            }
            delta = std::sqrt(delta);
        }
        return x;
    }

    template<typename T>
    std::vector<T> gauss_seidel_solve(const mtrx::dense<T>& A, const std::vector<T>& b, const std::vector<T> x0, int N) {
        std::vector<T> x = x0;
        for (int cycle = 0; cycle < N; cycle++) {
            for (auto i = 0u; i < x.size(); i++) {
                T xi = 0;
                for (auto j = i + 1; j < A.width(); j++) {
                    xi += A(i, j)*x[j];
                }
                for (auto j = 0; j < i; j++) {
                    xi += A(i, j)*x[j];
                }
                x[i] = (b[i] - xi)/A(i, i);
            }
        }
        return x;
    }

    template<typename T>
    std::vector<T> gauss_seidel_solve(const mtrx::csr<T>& A, const std::vector<T>& b, const std::vector<T> x0, T tolerance, int N_max) {
        std::vector<T> x = x0;
        T delta = tolerance + 1;
        for (int cycle = 0; (delta > tolerance) && (cycle < N_max); cycle++) {
            std::vector<T> old_x = x;
            for (auto i = 0u; i < x.size(); i++) {
                T xi = 0;
                for (auto j = A.get_raw_rows()[i]; j < A.get_raw_rows()[i+1]; j++) {
                    if (A.get_raw_cols()[j] > i) {
                        xi += A.get_raw_values()[j]*x[A.get_raw_cols()[j]];
                    }
                }
                for (auto j = A.get_raw_rows()[i]; j < A.get_raw_rows()[i+1]; j++) {
                    if (A.get_raw_cols()[j] < i) {
                        xi += A.get_raw_values()[j]*x[A.get_raw_cols()[j]];
                    }
                }
                x[i] = (b[i] - xi)/A(i, i);
            }
            delta = 0;
            for (auto i = 0u; i < x.size(); i++) {
                delta += (x[i] - old_x[i])*(x[i] - old_x[i]);
            }
            delta = std::sqrt(delta);
        }
        return x;
    }

    template<typename T>
    std::vector<T> gauss_seidel_solve(const mtrx::csr<T>& A, const std::vector<T>& b, const std::vector<T> x0, int N) {
        std::vector<T> x = x0;
        for (int cycle = 0; cycle < N; cycle++) {
            for (auto i = 0u; i < x.size(); i++) {
                T xi = 0;
                for (auto j = A.get_raw_rows()[i]; j < A.get_raw_rows()[i+1]; j++) {
                    if (A.get_raw_cols()[j] > i) {
                        xi += A.get_raw_values()[j]*x[A.get_raw_cols()[j]];
                    }
                }
                for (auto j = A.get_raw_rows()[i]; j < A.get_raw_rows()[i+1]; j++) {
                    if (A.get_raw_cols()[j] < i) {
                        xi += A.get_raw_values()[j]*x[A.get_raw_cols()[j]];
                    }
                }
                x[i] = (b[i] - xi)/A(i, i);
            }
        }
        return x;
    }
}