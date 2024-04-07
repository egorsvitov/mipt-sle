#pragma once
#include "../../dense_vs_csr/src/matrices.hpp"
#include "../../dense_vs_csr/src/vector_operations.hpp"
#include <cmath>

namespace slv {
    template<typename T>
    std::vector<T>fastest_gradiend_descent(const mtrx::csr<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T tolerance, int N_max) {
        std::vector<T> x = x0;
        T delta = tolerance + 1;
        for (int cycle = 0; !(delta <= tolerance) && (cycle < N_max); cycle++) {
            std::vector<T> old_x = x;
            std::vector<T> r = A*x - b;
            T tau = (r*r)/(r*(A*r));
            auto delta_x = tau*r;
            x = x - delta_x;
            delta = std::sqrt((x-old_x)*(x - old_x));
            //std::cout << std::sqrt(x*x) << std::endl;
        }
        return x;
    }

    template<typename T>
    std::vector<T> sym_gauss_seidel(const mtrx::csr<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T tolerance, int N_max) {
        std::vector<T> x = x0;
        T delta = tolerance + 1;
        for (int cycle = 0; !(delta <= tolerance) && (cycle < N_max); cycle++) {
            auto old_x = x;
            for (int i = 0; i < x.size(); i++) {
                T sum = 0;
                for (int j = 0; j < A.height(); j++) {
                    if (i != j) {
                        sum += A(i, j)*x[j];
                    }
                }
                x[i] = (b[i] - sum)/A(i,i);
            }
            for (int i = x.size()-1; i >= 0; i--) {
                T sum = 0;
                for (int j = 0; j < A.height(); j++) {
                    if (i != j) {
                        sum += A(i, j)*x[j];
                    }
                }
                x[i] = (b[i] - sum)/A(i,i);
            }
            delta = std::sqrt((x - old_x)*(x - old_x));
        }
        return x;
    }
}