#include "../../dense_vs_csr/src/matrices.hpp"
#include "../../dense_vs_csr/src/vector_operations.hpp"
#include <cmath>

namespace slv {
    template<typename T>
    std::vector<T> conj_grad(const mtrx::csr<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T tolerance, int N_max) {
        std::vector<T> x = x0;
        std::vector<T> r = A*x0 - b;
        std::vector<T> d = r;
        T alpha = 0;
        T beta = 0;
        for (int i = 0; !(std::sqrt(d*d) < tolerance) && (i < N_max); i++) { // сходится за число итераций, равное размерности задачи
            alpha = (r*r)/(d*(A*d)); // n + n + n = O(n), если матрица разреженная
            x = x - alpha*d; // n + n = O(n)
            std::vector old_r = r; // O(n)  
            r = A*x - b; // n + n = O(n)
            beta = (r*r)/(old_r*old_r); // n + n = O(n)
            d = r + beta*d; // n + n = O(n)
        }
        return x;
    }
}