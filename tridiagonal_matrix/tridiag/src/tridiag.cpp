#include <vector>
#include <array>
#include <iostream>
#include "../include/tridiag.h"

tridiag::tridiag() : s(0) {}

tridiag::tridiag(unsigned int n) : s(n) {}

tridiag::tridiag(std::vector<double> a, std::vector<double> b, std::vector<double> c) : s(b.size()) {
    data[0] = a;
    data[1] = b;
    data[2] = c;
}

unsigned int tridiag::size() const {return s;}

std::vector<double> tridiag::solve(std::vector<double> f) {
    std::vector<double> p(this->s);
    std::vector<double> q(this->s);
    std::vector<double> x(this->s);

    p[0] = -this->data[2][0]/this->data[1][0];
    q[0] = f[0]/this->data[1][0];

    for (int i = 1; i < this->s; i++) {
        p[i] = - this->data[2][i-1]/(this->data[0][i-1]*p[i-1] + this->data[1][i-1]);
        q[i] = (f[i-1] - this->data[0][i-2]*q[i-1])/(this->data[0][i-2]*p[i-1] + this->data[1][i-1]);
    }

    x[this->s-1] = (f[this->s-1] - this->data[0][this->s-2]*q[this->s-1])/(this->data[0][this->s-2]*p[this->s-1] + this->data[1][this->s-1]);
    
    for (int i = this->s-2; i >= 0; i--) {
        x[i] = p[i+1]*x[i+1] + q[i+1];
    }

    return x;
}