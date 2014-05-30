#include "../src/mathUtils.h"

#include "gtest/gtest.h"

class TestVectorOps : public testing::Test
{
public:
    TestVectorOps()
    {
        size_t len = 20;
        data.resize(len);
        for (size_t i=0; i<len; i++) {
            data[i] = i;
        }
    }

    dvector data;
    dvector baddata;
};

TEST_F(TestVectorOps, maxval)
{
    EXPECT_EQ(19, mathUtils::maxval(data));
    EXPECT_EQ(19, mathUtils::maxval(data, 0, 19));
    EXPECT_EQ(10, mathUtils::maxval(data, 5, 10));

    data[3] = 100;
    EXPECT_EQ(100, mathUtils::maxval(data));
    EXPECT_EQ(100, mathUtils::maxval(data, 0, 10));
    EXPECT_EQ(19, mathUtils::maxval(data, 5, 19));
}

TEST_F(TestVectorOps, minval)
{
    EXPECT_EQ(0, mathUtils::minval(data));
    EXPECT_EQ(0, mathUtils::minval(data, 0, 19));
    EXPECT_EQ(5, mathUtils::minval(data, 5, 10));

    data[3] = -100;
    EXPECT_EQ(-100, mathUtils::minval(data));
    EXPECT_EQ(-100, mathUtils::minval(data, 0, 10));
    EXPECT_EQ(5, mathUtils::minval(data, 5, 19));
}

TEST(TestSplines, integrate_eigen)
{
    dvec x = dvec::LinSpaced(31, 0.0, 3.0);
    dvec y = x.pow(2);
    EXPECT_NEAR(mathUtils::integrate(x, y), 9.0, 1e-4);
}

TEST(TestSplines, integrate_stdvector)
{
    std::vector<double> x(31), y(31);
    for (size_t i=0; i < 31; i++) {
        x[i] = 0.1 * i;
        y[i] = pow(x[i], 2);
    }
    EXPECT_NEAR(mathUtils::integrate(x, y), 9.0, 1e-4);
}
