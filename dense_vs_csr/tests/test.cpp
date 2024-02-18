#include <iostream>
#include "../src/matrices.hpp"
#include "../src/vector_operations.hpp"

int main() {
    mtrx::dense<float> A(3, 3, {1, 2, 3, 
                                2, 0, 3, 
                                3, 1, 2});
    std::vector<float> b = {1, 2, 3};
    std::vector<float> x = A*b;
    for (int i = 0; i < A.height(); i++) {
        std::cout << x[i] << ' ';
    }
    std::cout << '\n';
    mtrx::csr<float> B(A);
    B.print();
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
    mtrx::csr<float> C = 6.0f*B;
    C.print();

    std::vector<float> q = {1, 2, 3};
    std::vector<float> r = {4, 5, 6};
    std::cout << q << "| " << r << std::endl;
    std::cout << q*r << ' ' << r*q << std::endl;
    std::cout << q + r << std::endl;
    std::cout << q*3.0f << ' ' << 3.0f*q << std::endl;
}