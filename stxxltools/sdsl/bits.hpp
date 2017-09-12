/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file bits.hpp
    \brief bits.hpp contains the sdsl::bits class.
	\author Simon Gog
*/
// --- EXCERPT ---
#ifndef INCLUDED_SDSL_BITS
#define INCLUDED_SDSL_BITS

#include <stdint.h> // for uint64_t uint32_t declaration
#include <iostream>// for cerr
#include <cassert>
#ifdef __BMI2__
#include <immintrin.h>
#endif
#ifdef __SSE4_2__
#include <xmmintrin.h>
#endif

//! Namespace for the succinct data structure library.
namespace sdsl{
namespace bits{

//! Lookup table for select on bytes.
/*! Entry at idx = 256*j + i equals the position of the
    (j+1)-th set bit in byte i. Positions lie in the range \f$[0..7]\f$.
 */
constexpr uint8_t lt_sel[] = {
    0, 1, 1, 2, 1, 2, 2, 3,
    1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7,
    5, 6, 6, 7, 6, 7, 7, 8
};


//! Use to help to decide if a prefix sum stored in a byte overflows.
constexpr uint64_t ps_overflow[] = {
    0x8080808080808080ULL,
    0x7f7f7f7f7f7f7f7fULL,
    0x7e7e7e7e7e7e7e7eULL,
    0x7d7d7d7d7d7d7d7dULL,
    0x7c7c7c7c7c7c7c7cULL,
    0x7b7b7b7b7b7b7b7bULL,
    0x7a7a7a7a7a7a7a7aULL,
    0x7979797979797979ULL,
    0x7878787878787878ULL,
    0x7777777777777777ULL,
    0x7676767676767676ULL,
    0x7575757575757575ULL,
    0x7474747474747474ULL,
    0x7373737373737373ULL,
    0x7272727272727272ULL,
    0x7171717171717171ULL,
    0x7070707070707070ULL,
    0x6f6f6f6f6f6f6f6fULL,
    0x6e6e6e6e6e6e6e6eULL,
    0x6d6d6d6d6d6d6d6dULL,
    0x6c6c6c6c6c6c6c6cULL,
    0x6b6b6b6b6b6b6b6bULL,
    0x6a6a6a6a6a6a6a6aULL,
    0x6969696969696969ULL,
    0x6868686868686868ULL,
    0x6767676767676767ULL,
    0x6666666666666666ULL,
    0x6565656565656565ULL,
    0x6464646464646464ULL,
    0x6363636363636363ULL,
    0x6262626262626262ULL,
    0x6161616161616161ULL,
    0x6060606060606060ULL,
    0x5f5f5f5f5f5f5f5fULL,
    0x5e5e5e5e5e5e5e5eULL,
    0x5d5d5d5d5d5d5d5dULL,
    0x5c5c5c5c5c5c5c5cULL,
    0x5b5b5b5b5b5b5b5bULL,
    0x5a5a5a5a5a5a5a5aULL,
    0x5959595959595959ULL,
    0x5858585858585858ULL,
    0x5757575757575757ULL,
    0x5656565656565656ULL,
    0x5555555555555555ULL,
    0x5454545454545454ULL,
    0x5353535353535353ULL,
    0x5252525252525252ULL,
    0x5151515151515151ULL,
    0x5050505050505050ULL,
    0x4f4f4f4f4f4f4f4fULL,
    0x4e4e4e4e4e4e4e4eULL,
    0x4d4d4d4d4d4d4d4dULL,
    0x4c4c4c4c4c4c4c4cULL,
    0x4b4b4b4b4b4b4b4bULL,
    0x4a4a4a4a4a4a4a4aULL,
    0x4949494949494949ULL,
    0x4848484848484848ULL,
    0x4747474747474747ULL,
    0x4646464646464646ULL,
    0x4545454545454545ULL,
    0x4444444444444444ULL,
    0x4343434343434343ULL,
    0x4242424242424242ULL,
    0x4141414141414141ULL,
    0x4040404040404040ULL
};

// ============= inline - implementations ================

// see page 11, Knuth TAOCP Vol 4 F1A
inline uint64_t cnt(uint64_t x)
{
    return __builtin_popcountll(x);
}

inline uint32_t sel(uint64_t x, uint32_t i)
{
    uint64_t s = x, b;
    s = s-((s>>1) & 0x5555555555555555ULL);
    s = (s & 0x3333333333333333ULL) + ((s >> 2) & 0x3333333333333333ULL);
    s = (s + (s >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    s = 0x0101010101010101ULL*s;
// now s contains 8 bytes s[7],...,s[0]; s[j] contains the cumulative sum
// of (j+1)*8 least significant bits of s
    b = (s+ps_overflow[i]) & 0x8080808080808080ULL;
// ps_overflow contains a bit mask x consisting of 8 bytes
// x[7],...,x[0] and x[j] is set to 128-j
// => a byte b[j] in b is >= 128 if cum sum >= j

// __builtin_ctzll returns the number of trailing zeros, if b!=0
    int  byte_nr = __builtin_ctzll(b) >> 3;   // byte nr in [0..7]
    s <<= 8;
    i -= (s >> (byte_nr<<3)) & 0xFFULL;
    return (byte_nr << 3) + lt_sel[((i-1) << 8) + ((x>>(byte_nr<<3))&0xFFULL) ];
}

}} // end namespace sdsl

#endif

