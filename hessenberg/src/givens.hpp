#pragma once
#include "../../dense_vs_csr/src/matrices.hpp"
#include "../../dense_vs_csr/src/matrices.hpp"
#include "arnoldi.hpp"
#include <cmath>
#include <utility>

namespace givens {

    template<typename T>
    mtrx::csr<T> get_Gn_to_zero_last(const std::vector<std::vector<T>>& hs, const mtrx::dense<T>& Omega_n) {
        T rho = 0;
        T sigma = hs.back().back();
        for (int i = 0; i < Omega_n.width(); i++) {
            rho += Omega_n(Omega_n.height()-1, i) * hs.back()[i];
        }
        T c = rho/std::sqrt(rho*rho + sigma*sigma);
        T s = -sigma/std::sqrt(rho*rho + sigma*sigma);
        std::vector<T> vals;
        std::vector<std::size_t> cols;
        vals.reserve(Omega_n.height()+3);
        cols.reserve(Omega_n.height()+3);
        std::vector<std::size_t> rows(Omega_n.height()+2);
        std::size_t count = 0;
        rows[0] = 0;
        for (int i = 0; i < Omega_n.height() - 1; i++) {
            count++;
            vals.push_back(1);
            cols.push_back(i);
            rows[i+1] = count;
        }
        vals.push_back(c);
        cols.push_back(Omega_n.height()-1);
        vals.push_back(s);
        cols.push_back(Omega_n.height());
        count += 2;
        rows[Omega_n.height()] = count;
        vals.push_back(-s);
        cols.push_back(Omega_n.height()-1);
        vals.push_back(c);
        cols.push_back(Omega_n.height());
        count += 2;
        rows[Omega_n.height()+1] = count;
        std::cout << vals << std::endl;
        std::cout << cols << std::endl;
        std::cout << rows << std::endl;
        return mtrx::csr<T>(Omega_n.height()+1, Omega_n.height()+1, vals, cols, rows);
    }

    template<typename T>
    std::pair<T, T> get_c_s_to_zero_last(const std::vector<std::vector<T>>& hs, const mtrx::dense<T>& Omega_n) {
        T rho = 0;
        T sigma = hs.back().back();
        for (int i = 0; i < Omega_n.width(); i++) {
            rho += Omega_n(Omega_n.height()-1, i) * hs.back()[i];
        }
        T c = rho/std::sqrt(rho*rho + sigma*sigma);
        T s = -sigma/std::sqrt(rho*rho + sigma*sigma);
        
        return std::pair(c, s);
    }

    template<typename T>
    std::pair<T, T> zero_last(std::vector<std::vector<T>>& hs, const mtrx::dense<T>& Omega_n) {
        T rho = hs.back()[hs.back().size()-2];
        T sigma = hs.back().back();
        T c = rho/std::sqrt(rho*rho + sigma*sigma);
        T s = sigma/std::sqrt(rho*rho + sigma*sigma);
        hs.back()[Omega_n.height()] = 0;
        hs.back()[Omega_n.height()-1] = std::sqrt(rho*rho + sigma*sigma);

        return std::pair(c, s);
    }

}