// to compile g++-11 -std=c++23 -O2 -Wall -Wextra -Wpedantic test_file.cpp -lgtest -lgtest_main

#include <iostream>
#include <ranges>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <concepts>

#include <boost/tokenizer.hpp>
#include <boost/type_index.hpp>

#include <gtest/gtest.h>

template<typename T, typename U>
concept Addable = requires (T a, U b) { a + b; };


TEST(test1, try_something) {
    std::stringstream s;
    std::string str = "aaa-bbb-ccc/ddd;eee";
    boost::char_separator<char> sep("-;/");
    boost::tokenizer<boost::char_separator<char>> tokens(str, sep);
    //std::copy(tokens.begin(), tokens.end(), std::ostream_iterator<std::string>(s, " "));
    std::ranges::copy(tokens, std::ostream_iterator<std::string>(s, " "));
    EXPECT_EQ(s.str(), "aaa bbb ccc ddd eee ");
}

TEST(test2, try_something) {
    auto f = []<typename T, typename U> (T x, U y)
        requires Addable<T, U> { return x + y; };
    int res = f(2, 3);
    EXPECT_EQ(res, 5);
    std::cout << std::endl << boost::typeindex::type_id_with_cvr<decltype(f)>().pretty_name();
}
