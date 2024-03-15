#pragma once
#include "../../dense_vs_csr/src/matrices.hpp"
#include "../../dense_vs_csr/src/vector_operations.hpp"
#include <cmath>
#include <vector>

namespace slv {

    std::vector<int> rearrange_indxs(int r) {
        int n = std::pow(2, r);
        std::vector<int> indxs(n);
        int step = n;
        indxs[0] = 0;
        indxs[step/2] = 1;
        int p = 2;
        for (int k = 1; k < r; k++) {
            step /= 2;
            p *= 2;
            for (int i = 0; i < n; i+=step) {
                indxs[i+step/2] = p - 1 - indxs[i];
            }
        }
        return indxs;
    }

    template<typename T, typename M>
    std::vector<T> chebyshev_simple_solve(const M& A, const std::vector<T>& b, const std::vector<T> x0, int r, T lambda_min, T lambda_max) {
        std::vector<T> x = x0;
        int N = std::pow(2, r);
        auto indxs = rearrange_indxs(r);
        std::vector<T> taus(N);
        for (int i = 0; i < N; i++) {
            taus[indxs[i]] = 1/((lambda_min + lambda_max)/2 + (lambda_max - lambda_min)*std::cos(M_PI*(2*i+1)/(2*N))/2);
        }
        // std::cout << taus << std::endl;
        for (int cycle = 0; cycle < N; cycle++) {
            T delta = 0;
            auto delta_x = taus[cycle]*(A*x - b);
            for (int i = 0; i < delta_x.size(); i++) {
                delta += delta_x[i]*delta_x[i];
            }
            x = x - delta_x;
        }
        return x;
    }

    template<typename T, typename M>
    std::vector<T> chebyshev_simple_solve(const M& A, const std::vector<T>& b, const std::vector<T> x0, T tolerance, int r_max, T lambda_min, T lambda_max) {
        std::vector<T> x = x0;
        int N = std::pow(2, r_max);
        auto indxs = rearrange_indxs(r_max);
        std::vector<T> taus(N);
        for (int i = 0; i < N; i++) {
            taus[indxs[i]] = 1/((lambda_min + lambda_max)/2 + (lambda_max - lambda_min)*std::cos(M_PI*(2*i+1)/(2*N))/2);
        }
        // std::cout << taus << std::endl;
        T delta = tolerance + 1;
        for (int cycle = 0; (cycle < N) && (delta > tolerance); cycle++) {
            delta = 0;
            auto delta_x = taus[cycle]*(A*x - b);
            for (int i = 0; i < delta_x.size(); i++) {
                delta += delta_x[i]*delta_x[i];
            }
            x = x - delta_x;
        }
        return x;
    }
}