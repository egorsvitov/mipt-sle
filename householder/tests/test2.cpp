#include <vector>
#include <iostream>

template<typename T>
T operator*(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    T res = 0;
    for (int i = 0; i < rhs.size(); i++) {
        res += lhs[i] * rhs[i];       
    }
    return res;
}

int main() {
    std::vector<float> x = {1, 2, 3};
    std::vector<float> y = {3, 4, 5};
    std::cout << x*y << std::endl;
}