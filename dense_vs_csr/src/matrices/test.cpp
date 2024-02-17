#include <iostream>
#include "matrices.hpp"

int main() {
    mtrx::dense<float> A(3, 3, {1, 2, 3, 
                                2, 1, 3, 
                                3, 1, 2});
    std::vector<float> b = {1, 2, 3};
    std::vector<float> x = A*b;
    for (int i = 0; i < A.height(); i++) {
        std::cout << x[i] << ' ';
    }
    std::cout << '\n';
    mtrx::csr<float> B(A);
    for (int i =0; i < 3; i++) {
        for (int j = 0; j< 3; j++) {
            std::cout << B(i, j) << ' ';
        }
        std::cout << std::endl;
    }
    for (float i : B.get_raw_values()) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
    for (unsigned int i : B.get_raw_cols()) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
    for (unsigned int i : B.get_raw_rows()) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
}