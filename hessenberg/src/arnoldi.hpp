#pragma once
#include "../../dense_vs_csr/src/matrices.hpp"
#include "../../dense_vs_csr/src/matrices.hpp"
#include <cmath>
#include <utility>

namespace krylov {

    template<typename T>
    std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>> arnoldi_onb(const mtrx::csr<T>& A, const std::vector<T>& r0, std::size_t n) {
        std::vector<std::vector<T>> vs(n+1);
        vs.reserve(A.height()+1);
        std::vector<std::vector<T>> hs(n);
        hs.reserve(A.height());
        int count = 2;
        for (auto i : hs) {
            i.reserve(count);
            count++;
        }
        auto v = r0*(1/std::sqrt(r0*r0));
        vs[0] = v;
        for (int i = 0; i < n; i++) {
            auto t = A*v;
            for (int j = 0; j <= i; j++) {
                T h_ji = vs[j]*t;
                hs[i].push_back(h_ji);
                t = t - h_ji*vs[j];
            }
            T h_ipi = std::sqrt(t*t);
            hs[i].push_back(h_ipi);
            vs[i+1] = t*(1/h_ipi);
        }
        return std::pair(hs, vs);
    }
    template<typename T>
    void arnoldi_iter(const mtrx::csr<T>& A, std::vector<std::vector<T>>& vs, std::vector<std::vector<T>>& hs) {
        auto t = A*vs.back();
        hs.push_back(std::vector<T>(hs.back().size()+1));
        for (int i = 0; i < vs.size(); i++) {
            T h_ji = vs[i]*t;
            hs.back()[i] = h_ji;
            t = t - h_ji*vs[i];
        }
        T h_ipi = std::sqrt(t*t);
        hs.back()[vs.size()] = h_ipi;
        vs.push_back(t*(1/h_ipi));
    }

}

namespace mtrx {

    template<typename T>
    csr<T> form_H_CSR(const std::vector<std::vector<T>>& hs) {
        int w = hs.size();
        std::vector<T> values;
        values.reserve((w*w+w)/2);
        std::vector<std::size_t> cols;
        cols.reserve((w*w+w)/2);
        std::vector<std::size_t> rows(hs.back().size()+1, 0);
        for (int i = 0; i < hs.back().size(); i++) {
            rows[i+1]=rows[i];
            for (int j = ((i-2 < 0) ? 0 : i - 1); j < w; j++) {
                values.push_back(hs[j][i]);
                std::cout << hs[j][i] << " " << j << " " << i << std::endl;
                cols.push_back(j);
                rows[i+1]++;
            }
        }
        std::cout << values << std::endl;
        std::cout << cols << std::endl;
        std::cout << rows << std::endl;
        return csr<T>(hs.back().size(), hs.size(), values, cols, rows);
    }
    template<typename T>
    dense<T> form_H_dense(const std::vector<std::vector<T>>& hs) {
        std::vector<T> data(hs.back().size() * hs.size());
        for (int i = 0; i < hs.back().size(); i++) {
            for (int j = 0; j < hs.size(); j++) {
                if (i <= j + 1) {
                    data[i*hs.size()+j] = hs[j][i];
                } else {
                    data[i*hs.size()+j] = 0;
                }
            }
        }
        return dense<T>(hs.back().size(), hs.size(), data);
    }

}