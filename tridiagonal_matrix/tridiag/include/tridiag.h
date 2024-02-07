#include <vector>
#include <array>

class tridiag {
    unsigned int s;
    std::array<std::vector<double>, 3> data;
    public:
    tridiag();
    tridiag(unsigned int);
    tridiag(std::vector<double>, std::vector<double>, std::vector<double>);
    unsigned int size() const;
    std::vector<double> solve(std::vector<double>);
};