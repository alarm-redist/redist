/*
 * BigInteger.hpp
 *
 * Copyright (c) 2012 -- 2019 Jun Kawahara
 */


#ifndef BIGINTEGER_HPP
#define BIGINTEGER_HPP

#include <iostream>
#include <vector>
#include <string>
#include <climits>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <sstream>

class BigDecimal {
private:
    std::vector<int> place_list_;
public:
    BigDecimal()
    {
        place_list_.resize(1);
        place_list_[0] = 0;
    }

    BigDecimal& operator=(unsigned long long int integer)
    {
        place_list_.clear();
        while (integer > 0) {
            place_list_.push_back(integer % 10);
            integer /= 10;
        }
        return *this;
    }

    BigDecimal& operator+=(const BigDecimal& integer)
    {
        if (place_list_.size() < integer.place_list_.size()) {
            place_list_.resize(integer.place_list_.size());
        }
        unsigned long long int c = 0;
        for (unsigned int i = 0; i < place_list_.size(); ++i) {
            unsigned long long int q1 = place_list_[i];
            unsigned long long int q2 = (i < integer.place_list_.size() ? integer.place_list_[i] : 0);
            place_list_[i] = (q1 + q2 + c) % 10;
            c = (q1 + q2 + c) / 10;
        }
        if (c > 0) {
            place_list_.push_back(c);
        }
        return *this;
    }

    std::string GetString() const
    {
        std::ostringstream oss;
        for (int i = static_cast<int>(place_list_.size()) - 1; i >= 0; --i) {
            oss << place_list_[i];
        }
        return oss.str();
    }
};

namespace BigIntegerUtility {

// http://developer.cybozu.co.jp/takesako/2006/11/binary_hacks.html
int CountBits(unsigned long long int bit) {
    bit = ((bit & 0xaaaaaaaaaaaaaaaaUL) >> 1)
        + (bit & 0x5555555555555555UL);
    bit = ((bit & 0xccccccccccccccccUL) >> 2)
        + (bit & 0x3333333333333333UL); 
    bit = ((bit & 0xf0f0f0f0f0f0f0f0UL) >> 4)
        + (bit & 0x0f0f0f0f0f0f0f0fUL);
    bit = ((bit & 0xff00ff00ff00ff00UL) >> 8)
        + (bit & 0x00ff00ff00ff00ffUL);
    bit = ((bit & 0xffff0000ffff0000UL) >> 16)
        + (bit & 0x0000ffff0000ffffUL);
    bit = ((bit & 0xffffffff00000000UL) >> 32)
        + (bit & 0x00000000ffffffffUL);
    return bit;
}

// http://www.nminoru.jp/~nminoru/programming/bitcount.html
int GetNumberOfLeadingZero(unsigned long long int bit) {
    bit |= (bit >> 1);
    bit |= (bit >> 2);
    bit |= (bit >> 4);
    bit |= (bit >> 8);
    bit |= (bit >> 16);
    bit |= (bit >> 32);
    return CountBits(~bit);
}

}


class BigInteger {
private:
    std::vector<unsigned long long int> place_list_;
public:
    BigInteger()
    {
        place_list_.resize(1);
        place_list_[0] = 0;
    }

    BigInteger(unsigned long long int integer)
    {
        place_list_.resize(1);
        place_list_[0] = integer;
    }

    BigInteger& operator=(unsigned long long int integer)
    {
        place_list_.resize(1);
        place_list_[0] = integer;
        return *this;
    }

    BigInteger& operator+=(const BigInteger& integer)
    {
        if (place_list_.size() < integer.place_list_.size()) {
            place_list_.resize(integer.place_list_.size());
        }
        unsigned long long int c = 0;
        for (unsigned int i = 0; i < place_list_.size(); ++i) {
            unsigned long long int q1 = place_list_[i];
            unsigned long long int q2 = (i < integer.place_list_.size() ? integer.place_list_[i] : 0);
            place_list_[i] = static_cast<unsigned long long int>(q1 + q2 + c);
            if (q1 > ULLONG_MAX - q2) { // check overflow (whether q1 + q2 is larger than 2^64 - 1)
                c = 1;
            } else if (c > ULLONG_MAX - (q1 + q2)) {
                c = 1;
            } else {
                c = 0;
            }
        }
        if (c > 0) {
            place_list_.push_back(c);
       }
       return *this;
    }

    // If this < integer, throw exception
    BigInteger& operator-=(const BigInteger& integer)
    {
        int large_digit1, large_digit2, small_digit1, small_digit2;

        GetDigit(&large_digit1, &small_digit1);
        integer.GetDigit(&large_digit2, &small_digit2);

        if (IsSmaller(integer, large_digit1, small_digit1,
                      large_digit2, small_digit2)) {
            throw std::range_error("integer1 must not be smaller than integer2!");
        }

        BigInteger bi = integer;
        bi.place_list_.resize(large_digit1);
        for (size_t i = 0; i < large_digit1; ++i) {
            bi.place_list_[i] = ~bi.place_list_[i];
        }
        *this += bi;
        *this += 1;
        place_list_.resize(large_digit1);
        return *this;
    }

    operator double() const
    {
        double d = 0.0;
        for (int i = static_cast<int>(place_list_.size()) - 1; i >= 0; --i) {
            d *= 18446744073709551616.0; // d *= 2^64
            d += place_list_[i];
        }
        return d;
    }

    bool IsZero() const
    {
        for (size_t i = 0; i < place_list_.size(); ++i) {
            if (place_list_[i] != 0) {
                return false;
            }
        }
        return true;
    }

    // Check whether "value" is equal to this integer.
    // If this integer is larger than 2^64 - 1, always return false.
    bool Equals(unsigned long long int value) const
    {
        // check if the upper digits are 0
        for (size_t i = 1; i < place_list_.size(); ++i) {
            if (place_list_[i] != 0) {
                return false;
            }
        }
        return place_list_[0] == value;
    }

    void GetDigit(int* large_digit, int* small_digit) const
    {
        int pos = static_cast<int>(place_list_.size()) - 1;
        for (; pos >= 1; --pos) {
            if (place_list_[pos] != 0) {
                break;
            }
        }
        *large_digit = pos + 1;
        *small_digit = 64 - BigIntegerUtility::GetNumberOfLeadingZero(place_list_[pos]);
    }
    
    void Set(int place, unsigned long long int value)
    {
        if (place >= static_cast<int>(place_list_.size())) {
            place_list_.resize(place + 1);
        }
        place_list_[place] = value;
    }

    std::string GetString() const
    {
        BigDecimal d;
        for (int i = static_cast<int>(place_list_.size()) - 1; i >= 0; --i) {
            for (int j = 0; j < 64; ++j) { // compute d * 2^64
                BigDecimal d2 = d;
                d += d2;
            }
            BigDecimal d3;
            d3 = place_list_[i];
            d += d3;
        }
        return d.GetString();
    }

    bool IsSmaller(const BigInteger& integer2,
                   int large_digit1, int small_digit1,
                   int large_digit2, int small_digit2) const
    {
        if (large_digit1 < large_digit2) {
            return true;
        } else if (large_digit1 > large_digit2) {
            return false;
        } else {
            for (int i = large_digit1 - 1; i >= 0; --i) {
                if (place_list_[i] < integer2.place_list_[i]) {
                    return true;
                } else if (place_list_[i] > integer2.place_list_[i]) {
                    return false;
                }
            }
            return false; // this == integer2
        }
    }

    friend void TestBigInteger();
    friend bool operator<(const BigInteger& integer1, const BigInteger& integer2);
};

inline std::ostream& operator<<(std::ostream& ost, const BigInteger& integer)
{
    ost << integer.GetString();
    return ost;
}

inline BigInteger operator+(const BigInteger& integer1, const BigInteger& integer2)
{
    BigInteger integer3(integer1);
    integer3 += integer2;
    return integer3;
}

inline BigInteger operator-(const BigInteger& integer1, const BigInteger& integer2)
{
    BigInteger integer3(integer1);
    integer3 -= integer2;
    return integer3;
}

inline bool operator<(const BigInteger& integer1, const BigInteger& integer2)
{
    int large_digit1, large_digit2, small_digit1, small_digit2;

    integer1.GetDigit(&large_digit1, &small_digit1);
    integer2.GetDigit(&large_digit2, &small_digit2);

    return integer1.IsSmaller(integer2, large_digit1, small_digit1,
                              large_digit2, small_digit2);
}

class BigIntegerRandom {

//extern "C" {
//#include "mt19937ar.h"
//}

private:
    std::mt19937 engine_;
    std::uniform_int_distribution<unsigned long long int> dis_;

public:

    BigIntegerRandom()
    {
        std::random_device rd;
        std::array<int, std::mt19937::state_size> seed_data;
        std::generate_n(seed_data.data(), seed_data.size(), std::ref(rd));
        std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
        engine_ = std::mt19937(seq);
        dis_ = std::uniform_int_distribution<unsigned long long int>(0, 0xffffffffffffffffll);
    }

    //void MTInit(unsigned long int s)
    //{
        //init_genrand(s);
    //}

    //void MTInitByArray(unsigned long int init_key[], int key_length)
    //{
    //    //init_by_array(init_key, key_length);
    //}

    unsigned long long int GenRandUInt64()
    {
        unsigned long long int x = dis_(engine_);
        //std::cerr << x << std::endl;
        return x;
        //return ((static_cast<unsigned long long int>(genrand_int32()) << 32) |
        //         static_cast<unsigned long long int>(genrand_int32()));
    }

    // get rand between [0,r)
    BigInteger GetRand(const BigInteger& r)
    {
        int large_digit;
        int small_digit;
        r.GetDigit(&large_digit, &small_digit);

        assert(small_digit > 0);

        //std::cout << "l" << large_digit << " " << small_digit << std::endl;

        BigInteger ret(0);
        do {
            for (int i = 0; i < large_digit - 1; ++i) {
                unsigned long long int x = GenRandUInt64();
                ret.Set(i, x);
            }
            unsigned long long int y = GenRandUInt64();
            if (small_digit < 64) {
                y &= ((1ull << small_digit) - 1);
            }
            ret.Set(large_digit - 1, y);
        } while (!(ret < r));

        return ret;
    }
};

#include <cassert>

void TestBigInteger()
{
    BigInteger zero(0);
    BigInteger one(1);
    BigInteger bigint(100);

    assert(zero.IsZero());
    assert(!bigint.IsZero());

    assert(bigint.GetString() == std::string("100"));

    int large_digit, small_digit;
    bigint.GetDigit(&large_digit, &small_digit);
    assert(large_digit == 1);
    assert(small_digit == 7);

    bigint.Set(1, 123);
    assert(bigint.GetString() == std::string("2268949521066274848868"));

    bigint.GetDigit(&large_digit, &small_digit);
    assert(large_digit == 2);
    assert(small_digit == 7);

    BigInteger bigint1(0xffffffffffffffffull);
    BigInteger bigint2(0);
    bigint2.Set(1, 1);
    assert(bigint1.GetString() == std::string("18446744073709551615"));
    assert(bigint2.GetString() == std::string("18446744073709551616"));
    assert(bigint1 < bigint2);
    assert(!(bigint2 < bigint1));
    assert(zero < one);
    assert(!(one < zero));

    assert((BigInteger(9) - BigInteger(3)).GetString() == std::string("6"));
    assert((bigint2 - bigint1).GetString() == std::string("1"));

    BigInteger bigint3(0);
    bigint3.Set(10, 1);
    BigInteger bigint4(0);
    bigint4.Set(5, 123);
    bigint4.Set(4, 456);

    assert((bigint3 - bigint4).GetString() == std::string("456244061762219521"
        "86411716057002913248932285072485599305791925178992751672086773865059"
        "12811317108673373224037633407018448137633626166156972838831343037636"
        "957482329550887253504122441468437594112"));

    //BigIntegerUtility::MTInit(12345678);
    BigIntegerRandom random;
    BigInteger hundred(100);
    for (int i = 0; i < 1000; ++i) {
        BigInteger c = random.GetRand(hundred);
        //std::cerr << c << std::endl;
        assert(!(c < zero));
        assert(c < hundred);
    }
    for (int i = 0; i < 1000; ++i) {
        BigInteger c = random.GetRand(bigint1);
        //std::cerr << c << std::endl;
        assert(!(c < zero));
        assert(c < bigint1);
    }
    for (int i = 0; i < 1000; ++i) {
        BigInteger c = random.GetRand(bigint2);
        //std::cerr << c << std::endl;
        assert(!(c < zero));
        assert(c < bigint2);
    }

    std::cerr << "pass!" << std::endl;
}
/*
void TestBigInteger2() {

    for (int i = 0; i < 100; ++i) {
        std::cout << i << "\t" << BigIntegerUtility::CountBits(i) << "\t"
                  << BigIntegerUtility::GetNumberOfLeadingZero(i) << "\n";
    }
    std::cout.flush();
    unsigned long long int c = 0xffffffffffffffffull;
    std::cout << c << "\t" << BigIntegerUtility::CountBits(c) << "\t"
                  << BigIntegerUtility::GetNumberOfLeadingZero(c) << "\n";
    --c;
    std::cout << c << "\t" << BigIntegerUtility::CountBits(c) << "\t"
                  << BigIntegerUtility::GetNumberOfLeadingZero(c) << "\n";
    c = 0x7fffffffffffffffull;
    std::cout << c << "\t" << BigIntegerUtility::CountBits(c) << "\t"
                  << BigIntegerUtility::GetNumberOfLeadingZero(c) << "\n";

    BigInteger bi(0);
    int large_digit;
    int small_digit;
    bi.GetDigit(&large_digit, &small_digit);
    std::cout << bi << "\t" << large_digit << "\t" << small_digit << "\n";
    bi += 1;
    bi.GetDigit(&large_digit, &small_digit);
    std::cout << bi << "\t" << large_digit << "\t" << small_digit << "\n";
    bi += 1;
    bi.GetDigit(&large_digit, &small_digit);
    std::cout << bi << "\t" << large_digit << "\t" << small_digit << "\n";
    bi = BigInteger(0xffffffffffffffffull);
    bi.GetDigit(&large_digit, &small_digit);
    std::cout << bi << "\t" << large_digit << "\t" << small_digit << "\n";
    bi.Set(0, 0);
    bi.Set(1, 1);
    bi.GetDigit(&large_digit, &small_digit);
    std::cout << bi << "\t" << large_digit << "\t" << small_digit << "\n";
    bi += 1;
    bi.GetDigit(&large_digit, &small_digit);
    std::cout << bi << "\t" << large_digit << "\t" << small_digit << "\n";

    BigInteger bi10;
    bi10.Set(10, 1);

    for (int i = 0; i < 20; ++i) {
        std::cout << BigIntegerUtility::GetRand(bi10).GetString() << "\n";
    }
    std::cout << RAND_MAX << std::endl;
}

void TestBigInteger3()
{
    std::ifstream ifs("bigintegerdata.txt");
    if (!ifs) {
        std::cerr << "file bigintegerdata.txt cannot be opened" << std::endl;
        exit(1);
    }
    std::vector<BigInteger> bi_vec;
    for (int i = 0; i < 1000; ++i) {
        std::string bignum;
        std::string str;

        std::getline(ifs, bignum);
        //std::cout << bignum << std::endl;;
        std::getline(ifs, str);
        int n = atoi(str.c_str());
        BigInteger bi;
        for (int j = 0; j < n; ++j) {
            std::getline(ifs, str);
            std::istringstream iss(str);
            unsigned long long int c;
            iss >> c;
            bi.Set(j, c);
        }

        if (bi.GetString() != bignum) {
            std::cout << bi.GetString() << std::endl;
            std::cout << bignum << std::endl;
        }
        assert(bi.GetString() == bignum);
        bi_vec.push_back(bi);
        //std::cout << "i = " << i << " OK" << std::endl;
    }
    BigInteger sum;
    for (int i = 0; i < 1000; ++i) {
        sum += bi_vec[i];
        std::string bnum;
        std::getline(ifs, bnum);
        assert(sum.GetString() == bnum);
    }
    for (int i = 0; i < 1000; ++i) {
        sum -= bi_vec[i];
        std::string bnum;
        std::getline(ifs, bnum);
        assert(sum.GetString() == bnum);
    }
    std::string bstr;
    std::getline(ifs, bstr);
    for (int i = 0; i < 1000 - 1; ++i) {
        if (bi_vec[i] < bi_vec[i + 1]) {
            assert(bstr[i] == 't');
        } else {
            assert(bstr[i] == 'f');
        }
    }
}
*/

#endif // BIGINTEGER_HPP
