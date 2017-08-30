/**
 * @file    em_succinct_irreducible_src/compute_lcp_array.hpp
 * @section LICENCE
 *
 * This file is part of EM-SuccinctIrreducible v0.1.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2016
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

#ifndef __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_LCP_ARRAY_HPP_INCLUDED
#define __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_LCP_ARRAY_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <string>
#include <limits>
#include <algorithm>
#include <unistd.h>

#include "compute_plcp_bitvector.hpp"
#include "compute_lcp_from_plcp.hpp"
#include "utils.hpp"


namespace em_succinct_irreducible_private {

template<typename text_offset_type, typename ext_text_offset_type>
void compute_lcp_array(std::uint64_t text_length, std::uint64_t ram_use,
    std::string text_filename, std::string sa_filename, std::string bwt_filename,
    std::string output_filename, std::uint64_t &max_lcp, std::uint64_t &lcp_sum,
    std::uint64_t &n_irreducible_lcps, std::uint64_t &sum_irreducible_lcps,
    std::uint64_t &total_io_volume) {
  long double text_to_ram_ratio = (long double)text_length / (long double)ram_use;
  if (text_to_ram_ratio > 4.0L) {
    // Not enough RAM to hold B in RAM.
    std::string B_filename = output_filename + ".plcp." + utils::random_string_hash();
    compute_plcp_bitvector_small_ram<text_offset_type, ext_text_offset_type>(text_length, ram_use, text_filename,
        sa_filename, bwt_filename, B_filename, n_irreducible_lcps, sum_irreducible_lcps, total_io_volume);

    compute_lcp_from_plcp<text_offset_type>(text_length, ram_use, sa_filename,
        output_filename, B_filename, total_io_volume, max_lcp, lcp_sum);
  } else {
    // Enough RAM to hold B in RAM.
    std::uint64_t *B = compute_plcp_bitvector_large_ram<text_offset_type, ext_text_offset_type>(text_length,
        ram_use, text_filename, sa_filename, bwt_filename, output_filename, n_irreducible_lcps,
        sum_irreducible_lcps, total_io_volume);

    compute_lcp_from_plcp<text_offset_type>(text_length, ram_use, B, sa_filename,
        output_filename, total_io_volume, max_lcp, lcp_sum);
  }
}

template<typename text_offset_type, typename ext_text_offset_type>
void compute_lcp_array(std::string text_filename, std::string sa_filename,
    std::string bwt_filename, std::string output_filename, std::uint64_t ram_use) {
  srand(time(0) + getpid());
  utils::drop_disk_pages(text_filename);
  utils::drop_disk_pages(sa_filename);
  utils::drop_disk_pages(bwt_filename);
  long double global_start = utils::wclock();

  // Initialize basic parameters.
  std::uint64_t text_length = utils::file_size(text_filename);
  std::uint64_t lcp_sum = 0;
  std::uint64_t max_lcp = 0;
  std::uint64_t n_irreducible_lcps = 0;
  std::uint64_t sum_irreducible_lcps = 0;
  std::uint64_t total_io_volume = 0;
  long double text_to_ram_ratio = (long double)text_length / (long double)ram_use;

  if (text_length == 0) {
    fprintf(stderr, "Error: the input file is empty!\n");
    std::exit(EXIT_FAILURE);
  }

  // Turn paths absolute.
  text_filename = utils::absolute_path(text_filename);
  sa_filename = utils::absolute_path(sa_filename);
  bwt_filename = utils::absolute_path(bwt_filename);
  output_filename = utils::absolute_path(output_filename);

  // Print summary of basic parameters.
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "SA filename = %s\n", sa_filename.c_str());
  fprintf(stderr, "BWT filename = %s\n", bwt_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n", text_length, 1.L * text_length / (1 << 20));
  fprintf(stderr, "Text size / ram_use = %.2Lf\n", text_to_ram_ratio);
  fprintf(stderr, "RAM use = %lu (%.2LfMiB)\n", ram_use, ram_use / (1024.L * 1024));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n", sizeof(text_offset_type));
  fprintf(stderr, "sizeof(ext_text_offset_type) = %lu\n", sizeof(ext_text_offset_type));
#ifdef _OPENMP
  fprintf(stderr, "Max number of threads = %d\n", omp_get_max_threads());
#endif
  fprintf(stderr, "\n");

  compute_lcp_array<text_offset_type, ext_text_offset_type>(text_length, ram_use, text_filename,
      sa_filename, bwt_filename, output_filename, max_lcp, lcp_sum,
      n_irreducible_lcps, sum_irreducible_lcps, total_io_volume);

  // Print summary.
  long double total_time = utils::wclock() - global_start;
  long double avg_lcp = (long double)lcp_sum / text_length;
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  elapsed time = %.2Lfs (%.3Lfs/MiB of text)\n", total_time, total_time / (1.L * text_length / (1L << 20)));
  fprintf(stderr, "  speed = %.2LfMiB of text/s\n", (1.L * text_length / (1L << 20)) / total_time);
  fprintf(stderr, "  I/O volume = %lu (%.2Lfbytes/input symbol)\n", total_io_volume, (1.L * total_io_volume) / text_length);
  fprintf(stderr, "  number of irreducible LCPs = %lu\n", n_irreducible_lcps);
  fprintf(stderr, "  sum of irreducible LCPs = %lu\n", sum_irreducible_lcps);
  fprintf(stderr, "  sum of all LCPs = %lu\n", lcp_sum);
  fprintf(stderr, "  average LCP = %.2Lf\n", avg_lcp);
  fprintf(stderr, "  maximal LCP = %lu\n", max_lcp);
}

}  // namespace em_succinct_irreducible_private

#endif  // __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_LCP_ARRAY_HPP_INCLUDED
