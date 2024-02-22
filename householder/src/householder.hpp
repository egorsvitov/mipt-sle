#include "../../dense_vs_csr/src/matrices.hpp"
#include "../../dense_vs_csr/src/vector_operations.hpp"
#include <cmath>
#include <utility>
#include <algorithm>

namespace hldr {

    template<typename T>
    std::pair<mtrx::dense<T>, mtrx::dense<T>> hh_QR(const mtrx::dense<T>& A) {
        std::vector<std::vector<T>> vs(A.width());
        std::vector<T> R_data = A.get_data(); 
        for (int i = 0; i < A.width() - 1; i++) { 
            std::vector<T> x(A.height() - i);
            for (int l = 0; l < A.height() - i; l++) {
                x[l] = R_data[(l+i)*A.width()+i];
            }
            std::vector<T> v = x;
            v[0] += (x[0] >= 0) ? std::sqrt(x*x) : -std::sqrt(x*x); 
            vs[i] = v;
            T vx = v*x; 
            T vv = v*v;
            for (int l = 0; l < A.height() - i; l++) {
                R_data[(l+i)*A.width() + i] = x[l] - 2*(vx/vv)*v[l];
            }
            for (int j = i+1; j < A.width(); j++) {
                for (int l = 0; l < A.height() - i; l++) {
                    x[l] = R_data[(l+i)*A.width()+j];
                }
                vx = v*x;
                for (int l = 0; l < A.height() - i; l++) { 
                    R_data[(l+i)*A.width() + j] = x[l] - 2*(vx/vv)*v[l];
                }
            }
        }
        
        std::vector<T> Q_data(A.width()*A.height(), 0); 
        for (int i = 0; i < A.width(); i++) {
            Q_data[i*A.width() + i] = 1;
        }
        std::vector<T> diad_data(A.height()*A.height());
        for (int i = 0; i < A.height(); i++) {
            for (int j = 0; j < A.width(); j++) {
                diad_data[i*A.width()+j] = vs[0][i]*vs[0][j];
            }
        }
        Q_data = Q_data - 2/(vs[0]*vs[0])*diad_data;
        for (int i = 1; i < A.width() - 1; i++) { 
             for (int j = 0; j < A.height(); j++) {
                 std::vector<T> x(A.width() - i); 
                 for (int l = 0; l < A.width() - i; l++) {
                     x[l] = Q_data[j*A.width() + l+i];
                 }
                 T vx = vs[i]*x; 
                 T vv = vs[i]*vs[i];
                 for (int l = 0; l < A.height() - i; l++) {
                     Q_data[j*A.width() + l+i] = x[l] - 2*(vx/vv)*vs[i][l];
                } 
             }
        }
        return std::pair(mtrx::dense<T>(A.height(), A.width(), Q_data), mtrx::dense<T>(A.height(), A.width(), R_data));
    }

    template<typename T>
    std::vector<T> QR_solve(const mtrx::dense<T>& A, const std::vector<T>& b) {
        auto p = hldr::hh_QR(A);
        std::vector<T> Q_T_data(p.first.height()*p.first.width());
        for (int i = 0; i < p.first.height(); i++) {
            for (int j = 0; j < p.first.width(); j++) {
                Q_T_data[j*p.first.width() + i] = p.first.get_data()[i*p.first.width() + j];
            }
        }
        mtrx::dense<T> Q_T(p.first.width(), p.first.height(), Q_T_data);
        std::vector<T> x = Q_T*b;
        x[x.size() - 1] /= (p.second)(x.size() - 1, x.size() - 1);
        for (int i = x.size() - 2; i >= 0; i--) {
            for (int j = i + 1; j < x.size(); j++) {
                x[i] -= x[j]*(p.second)(i, j);

            }
            x[i] /= (p.second)(i, i);
        }
    return x;
    }
}