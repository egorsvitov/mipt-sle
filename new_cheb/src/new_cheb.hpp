#pragma once
#include "../../dense_vs_csr/src/matrices.hpp"
#include "../../dense_vs_csr/src/vector_operations.hpp"

namespace slv {

template<typename T>
std::vector<T> new_cheb_M(const mtrx::csr<T>& M, const std::vector<T>& b, const std::vector<T>& x0, T tolerance, int N_max, T tau, T rho) {
    
    std::vector<T> x = x0;
    T delta = tolerance + 1;

    auto ypp = x0;
    auto yp = M*x0 + tau*b;
    std::vector<T> y;
    T mup = 1;
    T mu = 1/rho;
    T mun;
    for (int cycle = 0; !(delta < tolerance) && (cycle < N_max); cycle++) {

        mun = 2*mu/rho - mup;

        y = (2*mu)/(mun*rho) * (M*yp + tau*b) - (mup/mun)*ypp;
        ypp = yp;
        yp = y;

        mup = mu;
        mu = mun;

        delta = std::sqrt((yp - ypp)*(yp - ypp));

    }
    return y;
}

template<typename T>
std::vector<T> new_cheb(const mtrx::csr<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T tolerance, int N_max, T tau, T rho) {

    std::vector<T> x = x0;
    T delta = tolerance + 1;

    std::vector<T> data(A.height()*A.width());
    for (int i = 0; i < A.height(); i++) {
        for (int j = 0; j < A.width(); j++) {
            data[i*A.width() + j] = A(i, j);
        }
    }
    mtrx::dense<T> B(A.height(), A.width(), data);

    mtrx::csr<T> M(mtrx::dense_ident<double>(B.height()) - tau*B);
    mtrx::csr<T> N(tau*mtrx::dense_ident<double>(B.height()));
    auto ypp = x0;
    auto yp = M*x0 + tau*b;
    std::vector<T> y;
    T mup = 1;
    T mu = 1/rho;
    T mun;
    for (int cycle = 0; !(delta < tolerance) && (cycle < N_max); cycle++) {

        // std::cout << mu << " " << mup << " " << mupp << std::endl;

        mun = 2*mu/rho - mup;

        y = (2*mu)/(mun*rho) * (M*yp + tau*b) - (mup/mun)*ypp;
        ypp = yp;
        yp = y;

        // std::cout << y << "|" << yp << "|" << ypp << std::endl;

        mup = mu;
        mu = mun;

    }
    return y;
}

}