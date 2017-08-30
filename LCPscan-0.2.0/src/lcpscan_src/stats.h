/**
 * @file    src/lcpscan_src/stats.h
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

#ifndef __LCPSCAN_SRC_STATS_H_INCLUDED
#define __LCPSCAN_SRC_STATS_H_INCLUDED


namespace lcpscan_private {

struct stats_t {
  stats_t() {
    max_lcp = 0;
    lcp_sum = 0;
    irr_lcp_count = 0;
    irr_lcp_sum = 0;
  }

  stats_t& operator = (const stats_t &s) {
    max_lcp = s.max_lcp;
    lcp_sum = s.lcp_sum;
    irr_lcp_count = s.irr_lcp_count;
    irr_lcp_sum = s.irr_lcp_sum;

    return *this;
  }

  long max_lcp;
  long lcp_sum;
  long irr_lcp_count;
  long irr_lcp_sum;
};

}  // namespace lcpscan_private

#endif  // __LCPSCAN_SRC_STATS_H_INCLUDED
