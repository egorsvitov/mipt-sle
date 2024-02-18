#include <iostream>
#include <random>
#include "../src/matrices.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>

int main()
{
    std::ofstream n;
	n.open("output.txt");
	n << '\n';
	n.close();

    std::ofstream out;
    out.open("output.txt", std::ios::app);
    for (double p = 1; p > 0.00005; p*=0.8)  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> distrib({1 - p, p});
    for (int n = 100; n <= 1000; n+=100) {
        std::vector<double> data(n*n);       
        for (int i = 0; i < n*n; i++) {
            data[i] = 3.0*distrib(gen);
        }
        mtrx::dense<double> A(n, n, data);
        std::vector<double> b(n);
        std::iota(b.begin(), b.end(), 1);
        
        std::cout << n << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        auto c = A*b;
        auto end = std::chrono::high_resolution_clock::now();
        out << std::to_string(p) << " " << std::to_string(n) << " " << std::to_string((end-start).count()) << " ";
        mtrx::csr<double> B(A);
        start = std::chrono::high_resolution_clock::now();
        c = B*b;
        end = std::chrono::high_resolution_clock::now();
        out << std::to_string((end-start).count()) << std::endl;
    }
    }
    out.close();
}