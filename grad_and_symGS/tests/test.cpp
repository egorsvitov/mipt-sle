#include "../src/grad_and_symGS.hpp"
#include "../../iterations/src/iter_solvers.hpp"

int main() {
    mtrx::dense<double> B(3, 3, {30, 2, 3, 
                                2, 20, 3, 
                                3, 3, 50});

    mtrx::csr<double> A(B);
    std::cout << slv::gauss_seidel_solve(A, {1, 2, 3}, {1, 1, 1}, 0.00001, 100) << std::endl;
    std::cout << slv::sym_gauss_seidel(A, {1, 2, 3}, {1, 1, 1}, 0.00001, 100) << std::endl;
    std::cout << slv::fastest_gradiend_descent(A, {1, 2, 3}, {1, 1, 1}, 0.00001, 100) << std::endl;
}