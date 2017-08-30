/**
 * @file    src/lcpscan_src/sort.cc
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

#include "./sort.h"

#include <cstdio>
#include <algorithm>
#include <string>

#include "./utils.h"
#include "./stream.h"
#include "./uint40.h"
#include "./stats.h"


namespace lcpscan_private {

//==============================================================================
// Return length of lcp between suffixes starting at positions i and j.
//==============================================================================
long naive_lcp(long i, long j, long lcp, std::FILE *f_text, long text_length) {
  static const long bufsize = (1 << 20);
  unsigned char *b1 = new unsigned char[bufsize];
  unsigned char *b2 = new unsigned char[bufsize];

  while (true) {
    long toread = std::min(bufsize, text_length - std::max(i, j) - lcp);
    if (!toread) break;

    utils::read_at_offset(b1, i + lcp, toread, f_text);
    utils::read_at_offset(b2, j + lcp, toread, f_text);

    long lcp2 = 0L;
    while (lcp2 < toread && b1[lcp2] == b2[lcp2])
      ++lcp2;
    lcp += lcp2;

    if (lcp2 < toread)
      break;
  }
  delete[] b1;
  delete[] b2;

  return lcp;
}

//==============================================================================
// Build LCP array in internal memory.
//==============================================================================
void im_laca(std::string text_fname, std::string sa_fname,
    std::string out_fname, long text_length, bool dry_run, stats_t &stats) {
  long *tab = new long[text_length];

  // Compute Phi and store in tab.
  stream_reader<uint40> *sa_reader = new stream_reader<uint40>(sa_fname);
  for (long i = 0, prev = text_length; i < text_length; ++i) {
    long sa_i = sa_reader->read();
    tab[sa_i] = prev;
    prev = sa_i;
  }
  delete sa_reader;

  // Read text.
  unsigned char *text = new unsigned char[text_length];
  utils::read_from_file(text, text_length, text_fname);

  // Compute PLCP (overwrite Phi).
  stats_t st;
  for (long i = 0, match = 0, prev_phi = text_length; i < text_length; ++i) {
    long j = tab[i]; // j = Phi[i]
    while (i + match < text_length && j + match < text_length &&
        text[i + match] == text[j + match]) ++match;
    tab[i] = match;
    st.max_lcp = std::max(st.max_lcp, match);
    st.lcp_sum += match;
    match = std::max(match - 1, 0L);

    // Update stats.
    if (i == 0 || j != prev_phi + 1 || tab[i - 1] == 0) {
      st.irr_lcp_count += 1;
      st.irr_lcp_sum += tab[i];
    }
    prev_phi = j;
  }
  delete[] text;
  stats = st;

  // Permute PLCP into LCP and write to disk.
  if (!dry_run) {
    sa_reader = new stream_reader<uint40>(sa_fname);
    stream_writer<uint40> *lcp_writer = new stream_writer<uint40>(out_fname);
    for (long i = 0; i < text_length; ++i) {
      long sa_i = sa_reader->read();
      long lcp_i = tab[sa_i];
      lcp_writer->write(lcp_i);
    }
    delete sa_reader;
    delete lcp_writer;
  }

  delete[] tab;
}

}  // namespace lcpscan_private
