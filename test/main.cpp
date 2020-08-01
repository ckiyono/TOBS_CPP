#include "gtest/gtest.h"
#include <string>
#include <sstream>
#include <cstdio>
#include <iostream>

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
	int ret = RUN_ALL_TESTS();

	std::cin.ignore();
	return ret;
}