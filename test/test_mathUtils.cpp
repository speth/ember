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
