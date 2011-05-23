#include <UnitTest++.h>
#include <iostream>

int main(int argc, char** argv)
{
    std::cout << "****** Running unit tests *******" << std::endl;
    int err = UnitTest::RunAllTests();
    std::cout << "****** Finished unit tests ******" << std::endl;
    return err;
}
