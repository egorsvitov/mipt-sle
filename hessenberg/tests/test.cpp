#include "../src/arnoldi.hpp"
#include "../src/givens.hpp"

int main() {
    mtrx::dense<double> B(3, 3, {30, 2, 3, 
                                2, 20, 3, 
                                3, 3, 50});

    mtrx::csr<double> A(B);
    auto p = krylov::arnoldi_onb(A, {1, 1, 1}, 1);
    //krylov::arnoldi_iter(A, p.second, p.first);

    double r = std::sqrt(12.9185*12.9185 + 38.6667*38.6667);
    double c = 38.6667/r;
    double s = -12.9185/r;
    mtrx::dense<double> Omega_1(2, 2, {c, -s, s, c});

    givens::zero_last(p.first, Omega_1);

    for (auto i : p.second) {
        std::cout << i << std::endl;
    }
    std::cout << std::endl;
    for (auto i : p.first) {
        std::cout << i << std::endl;
    }
    // mtrx::form_H_CSR(p.first).print();
    // mtrx::form_H_dense(p.first).print();
    
    // std::cout << std::endl;
    // givens::get_Gn_to_zero_last(p.first, Omega_1).print();
}