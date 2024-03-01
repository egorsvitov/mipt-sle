#include "../src/iter_solvers.hpp"
#include <chrono>
#include <cstdio>
#include <cmath>

int main() {

    std::freopen("by_iterations.txt", "w", stdout);

    mtrx::dense<double> B(3, 3, {30, 2, 3, 
                                2, 20, 3, 
                                3, 3, 50});
    
    mtrx::csr<double> A(B);
    std::vector<double> b = {1, 2, 3};
    std::vector<double> x0 = {1, 1, 1};
    double tau = 0.0285;

    for (int n = 1; n < 20; n++) {
            std::cout << n << " ";
            std::vector<double> x;
            x = slv::simple_solve<double>(A, b, x0, n, tau);
            std::cout << std::sqrt(x*x) << " ";
            x = slv::jacobi_solve<double>(A, b, x0, n);
            std::cout << std::sqrt(x*x) << " ";
            x = slv::gauss_seidel_solve<double>(A, b, x0, n);
            std::cout << std::sqrt(x*x) << std::endl;
    }

    std::freopen("by_time.txt", "w", stdout);
    for (double tolerance = 0.0001; tolerance > 0.00000001; tolerance *= 0.8) {
        std::cout << tolerance << " ";
        std::vector<double> x;
        auto start = std::chrono::high_resolution_clock::now();
        x = slv::simple_solve<double>(A, b, x0, tolerance, tau, 100);;
        auto end = std::chrono::high_resolution_clock::now();
        //std::cout << x << std::endl;
        std::cout << std::to_string((end-start).count()) << " ";

        start = std::chrono::high_resolution_clock::now();
        x = slv::jacobi_solve<double>(A, b, x0, tolerance, 100);;
        end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string((end-start).count()) << " ";

        start = std::chrono::high_resolution_clock::now();
        x = slv::gauss_seidel_solve<double>(A, b, x0, tolerance, 100);;
        end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string((end-start).count()) << std::endl;

    }
}