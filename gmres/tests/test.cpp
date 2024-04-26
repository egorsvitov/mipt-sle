#include "../src/gmres.hpp"
#include "../../iterations/src/iter_solvers.hpp"

int main() {
    mtrx::dense<double> B(3, 3, {30, 2, 3, 
                                2, 20, 3, 
                                3, 3, 50});

    mtrx::csr<double> A(B);

    std::cout << slv::gmres(A, {1, 2, 3}, {1, 1, 1}, 0.0001, 4) << std::endl;
    std::cout << slv::simple_solve(A, {1, 2, 3}, {0, 0, 0}, 100, 0.0285) << std::endl;

}