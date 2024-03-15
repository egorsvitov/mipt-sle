#pragma once
#include <vector>
#include <iostream>
#include <iterator>

template<typename T>
std::vector<T> operator+(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    std::vector<T> res(rhs.size());
    for (int i = 0; i < rhs.size(); i++) {
        res[i] = lhs[i] + rhs[i];       
    }
    return res;
}

template<typename T>
T operator*(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    T res = 0;
    for (int i = 0; i < rhs.size(); i++) {
        res += lhs[i] * rhs[i];       
    }
    return res;
}
template<typename T>
std::vector<T> operator*(const T& lhs, const std::vector<T>& rhs) {
    std::vector<T> res(rhs.size());
    for (int i = 0; i < rhs.size(); i++) {
        res[i] = lhs * rhs[i];       
    }
    return res;
}
template<typename T>
std::vector<T> operator*(const std::vector<T>& lhs, const T& rhs) {
    return rhs*lhs;
}

template<typename T>
std::vector<T> operator-(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    std::vector<T> res(rhs.size());
    for (int i = 0; i < rhs.size(); i++) {
        res[i] = lhs[i] - rhs[i];       
    }
    return res;
}
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}