/**
 * @file    src/lcpscan_src/tuples.h
 * @section LICENCE
 *
 * This file is part of LCPscan v0.2.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2016
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __LCPSCAN_SRC_TUPLES_H_INCLUDED
#define __LCPSCAN_SRC_TUPLES_H_INCLUDED

#include <algorithm>
#include <limits>

#include "./uint40.h"
#include "./radixsort.h"
#include "./radixsort_stxxl.h"


namespace lcpscan_private {

struct Triple {
  Triple() {}
  Triple(uint40 f, uint40 s, uint40 t)
    : first(f), second(s), third(t) {}
  uint40 first, second, third;
};

struct Pair {
  Pair() {}
  Pair(uint40 f, uint40 s)
    : first(f), second(s) {}
  uint40 first, second;
};

struct CmpTriple {
  inline bool operator() (const Triple &a, const Triple &b) const {
    return a.first < b.first;
  }

  Triple min_value() const {
    return Triple(std::numeric_limits<uint40>::min(), 0, 0);
  }

  Triple max_value() const {
    return Triple(std::numeric_limits<uint40>::max(), 0, 0);
  }
};

struct CmpPair {
  inline bool operator() (const Pair &a, const Pair &b) const {
    return a.first < b.first;
  }

  Pair min_value() const {
    return Pair(std::numeric_limits<uint40>::min(), 0);
  }

  Pair max_value() const {
    return Pair(std::numeric_limits<uint40>::max(), 0);
  }
};

}  // namespace lcpscan_private

namespace stxxl {
namespace algorithm {

template <typename Iterator>
struct Sort<Iterator, lcpscan_private::CmpPair> {
  static inline void sort(Iterator begin, Iterator end, const lcpscan_private::CmpPair& cmp) {
    size_t depth = sizeof(uint40) - 1;
    radixsort8msb_page(begin, end, cmp, depth);
  }
};

template <typename Iterator>
struct Sort<Iterator, lcpscan_private::CmpTriple> {
  static inline void sort(Iterator begin, Iterator end, const lcpscan_private::CmpTriple& cmp) {
    size_t depth = sizeof(uint40) - 1;
    radixsort8msb_page(begin, end, cmp, depth);
  }
};


}  // namespace algorithm
}  // namespace stxxl

#endif  // __LCPSCAN_SRC_TUPLES_H_INCLUDED
