#include "../../dense_vs_csr/src/matrices.hpp"
#include "../../dense_vs_csr/src/vector_operations.hpp"
#include "../src/new_cheb.hpp"
#include "../../iterations/src/iter_solvers.hpp"

int main() {
//    mtrx::csr<double>(3, 3, {1, 1, 1}, {0, 1, 2}, {0, 1, 2, 3}).print();
//    mtrx::csr_ident<double>(5).print();
    mtrx::dense<double> B(3, 3, {30, 2, 3, 
                                2, 20, 3, 
                                3, 3, 50});

    mtrx::csr<double> A(B);
    std::cout << slv::new_cheb(A, {1, 2, 3}, {0, 0, 0}, 0.00001, 100, 0.0285, 0.4467) << std::endl;
    std::cout << slv::simple_solve(A, {1, 2, 3}, {0, 0, 0}, 100, 0.0285) << std::endl;
}