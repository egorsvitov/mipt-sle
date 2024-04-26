#include "../../iterations/src/iter_solvers.hpp"
#include "../../new_cheb/src/new_cheb.hpp"
#include "../../chebyshev/src/chebyshev.hpp"
#include "../../grad_and_symGS/src/grad_and_symGS.hpp"
#include "../../conjugate_gradient/src/conjugate_gradient.hpp"
#include "../src/gmres.hpp"
#include <chrono>
#include <cstdio>
#include <cmath>

int main() {
    
    // mtrx::dense<double> B(3, 3, {19.423, 0, 0, 
    //                             0, 20, 0, 
    //                             0, 0, 50.786});
    mtrx::dense<double> B(3, 3, {30, 2, 3, 
                                2, 20, 3, 
                                3, 3, 50});

    mtrx::csr<double> A(B);
    std::vector<double> b = {1, 2, 3};
    std::vector<double> x0 = {1, 1, 1};
    double tau = 0.0285;

    std::vector<double> data(A.height()*A.width());
    for (int i = 0; i < A.height(); i++) {
        for (int j = 0; j < A.width(); j++) {
            data[i*A.width() + j] = A(i, j);
        }
    }

    mtrx::dense<double> D(A.height(), A.width(), data);

    mtrx::csr<double> M(mtrx::dense_ident<double>(D.height()) - tau*D);


    std::freopen("by_time.txt", "w", stdout);
    for (double tolerance = 0.0001; tolerance > 0.00000001; tolerance *= 0.8) {
        std::cout << tolerance << " ";
        std::vector<double> x;
        auto start = std::chrono::high_resolution_clock::now();
        x = slv::simple_solve<double>(A, b, x0, tolerance, tau, 1000);;
        auto end = std::chrono::high_resolution_clock::now();
        //std::cout << x << std::endl;
        std::cout << std::to_string((end-start).count()) << " ";

        start = std::chrono::high_resolution_clock::now();
        x = slv::jacobi_solve<double>(A, b, x0, tolerance, 1000);;
        end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string((end-start).count()) << " ";

        start = std::chrono::high_resolution_clock::now();
        x = slv::gauss_seidel_solve<double>(A, b, x0, tolerance, 1000);;
        end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string((end-start).count()) << " ";

        start = std::chrono::high_resolution_clock::now();
        x = slv::chebyshev_simple_solve<double>(A, b, x0, 10, tolerance, 19.423, 50.786);
        end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string((end-start).count()) << " ";

        start = std::chrono::high_resolution_clock::now();
        x = slv::new_cheb_M<double>(M, b, x0, tolerance, 1000, tau, 0.4467);
        end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string((end-start).count()) << " ";

        start = std::chrono::high_resolution_clock::now();
        x = slv::fastest_gradiend_descent<double>(A, b, x0, tolerance, 100);
        end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string((end-start).count()) << " ";

        start = std::chrono::high_resolution_clock::now();
        x = slv::sym_gauss_seidel<double>(A, b, x0, tolerance, 1000);
        end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string((end-start).count()) << " ";

        start = std::chrono::high_resolution_clock::now();
        x = slv::conj_grad<double>(A, b, x0, tolerance, 1000);
        end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string((end-start).count()) << " ";

        start = std::chrono::high_resolution_clock::now();
        x = slv::gmres<double>(A, b, x0, -1.0, 4);
        end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string((end-start).count()) << std::endl;
    }
}