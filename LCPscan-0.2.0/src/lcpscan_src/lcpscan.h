/**
 * @file    src/lcpscan_src/lcpscan.h
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

#ifndef __LCPSCAN_SRC_LCPSCAN_H_INCLUDED
#define __LCPSCAN_SRC_LCPSCAN_H_INCLUDED

#include <string>

#include "sort.h"


template<long disk_block_size>
void construct_lcp(
    std::string text_filename,
    std::string sa_filename,
    std::string out_filename,
    long ram_use,
    long n_sblock,
    bool dry_run,
    bool distribute_text,
    bool distribute_sa,
    bool disable_cache) {
  lcpscan_private::construct_lcp<disk_block_size>(text_filename,
      sa_filename, out_filename, ram_use, n_sblock, dry_run,
      distribute_text, distribute_sa, disable_cache);
}

#endif  // __LCPSCAN_SRC_LCPSCAN_H_INCLUDED
