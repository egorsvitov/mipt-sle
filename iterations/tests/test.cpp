#include "../src/iter_solvers.hpp"

int main() {
    mtrx::dense<double> B(3, 3, {30, 2, 3, 
                                2, 20, 3, 
                                3, 3, 50});

    mtrx::csr<double> A(B);
    std::cout << slv::jacobi_solve(A, {1, 2, 3}, {1, 1, 1}, 0.001, 100) << std::endl;
    std::cout << slv::simple_solve(A, {1, 2, 3}, {1, 1, 1}, 0.001, 0.0285, 100) << std::endl;
    std::cout << slv::gauss_seidel_solve(B, {1, 2, 3}, {1, 1, 1}, 0.001, 100) << std::endl;
    std::cout << slv::gauss_seidel_solve(A, {1, 2, 3}, {1, 1, 1}, 0.001, 100) << std::endl;
    std::cout << std::endl;
    std::cout << slv::jacobi_solve(A, {1, 2, 3}, {1, 1, 1}, 10) << std::endl;
    std::cout << slv::simple_solve(A, {1, 2, 3}, {1, 1, 1}, 10, 0.0285) << std::endl;
    std::cout << slv::gauss_seidel_solve(B, {1, 2, 3}, {1, 1, 1}, 10) << std::endl;
    std::cout << slv::gauss_seidel_solve(A, {1, 2, 3}, {1, 1, 1}, 10) << std::endl;
}