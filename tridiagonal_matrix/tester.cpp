#include "tridiag/include/tridiag.h"
#include <iostream>
#include <gtest/gtest.h>

TEST(test1, satisfies_conditions) {
    std::vector<double> f = {1, 2, 3};
    tridiag t({1, 1}, {2, 2, 2}, {1, 1});
    std::vector<double> x = t.solve(f);
    std::vector<double> true_x = {0.5, 0, 1.5};
    for (int i = 0; i < t.size(); i++)  {
        ASSERT_NEAR(true_x[i], x[i], 0.00001);
    }
}

TEST(test2, doesnt_satisfy_conditions) {
    std::vector<double> f = {1, 2, 3};
    tridiag t({4, 5}, {1, 1, 1}, {4, 4});
    std::vector<double> x = t.solve(f);
    std::vector<double> true_x = {-0.6, 0.4, 1};
    for (int i = 0; i < t.size(); i++)  {
        ASSERT_NEAR(true_x[i], x[i], 0.00001);
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS();
}
