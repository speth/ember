#include <UnitTest++.h>
#include "../src/mathUtils.h"

SUITE(MathUtils)
{

struct TestVectorOps
{
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

TEST_FIXTURE(TestVectorOps, maxval)
{
    CHECK_EQUAL(19, mathUtils::maxval(data));
    CHECK_EQUAL(19, mathUtils::maxval(data, 0, 19));
    CHECK_EQUAL(10, mathUtils::maxval(data, 5, 10));

    data[3] = 100;
    CHECK_EQUAL(100, mathUtils::maxval(data));
    CHECK_EQUAL(100, mathUtils::maxval(data, 0, 10));
    CHECK_EQUAL(19, mathUtils::maxval(data, 5, 19));
}

TEST_FIXTURE(TestVectorOps, minval) {
    CHECK_EQUAL(0, mathUtils::minval(data));
    CHECK_EQUAL(0, mathUtils::minval(data, 0, 19));
    CHECK_EQUAL(5, mathUtils::minval(data, 5, 10));

    data[3] = -100;
    CHECK_EQUAL(-100, mathUtils::minval(data));
    CHECK_EQUAL(-100, mathUtils::minval(data, 0, 10));
    CHECK_EQUAL(5, mathUtils::minval(data, 5, 19));
}

}
