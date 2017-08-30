/**
 * @file    em_succinct_irreducible_src/compute_plcp_bitvector.hpp
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

#ifndef __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_PLCP_BITVECTOR_HPP_INCLUDED
#define __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_PLCP_BITVECTOR_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <ctime>
#include <string>
#include <algorithm>
#include <omp.h>
#include <unistd.h>

#include "utils.hpp"
#include "distribute_pairs_and_compute_C.hpp"
#include "process_halfsegment_pairs.hpp"
#include "compute_B.hpp"


namespace em_succinct_irreducible_private {

// A version that returns the B bitvector as a file on disk.
template<typename text_offset_type, typename ext_text_offset_type>
void compute_plcp_bitvector_small_ram(std::uint64_t text_length, std::uint64_t ram_use,
    std::string text_filename, std::string sa_filename, std::string bwt_filename,
    std::string B_filename, std::uint64_t &n_irreducible_lcps,
    std::uint64_t &sum_irreducible_lcps, std::uint64_t &total_io_volume) {
  fprintf(stderr, "Compute PLCP bitvector (dest = EM):\n");
  long double compute_plcp_bitvector_start = utils::wclock();

  // Initialize basic parameters.
  static const std::uint64_t max_overflow_size = (1UL << 20);
  std::uint64_t max_halfsegment_size = std::max(1UL, ram_use / 2);
  std::uint64_t n_halfsegments = (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t n_different_halfsegment_pairs = (n_halfsegments * (n_halfsegments + 1)) / 2;
  std::uint64_t io_volume = 0;
  std::uint64_t max_block_size_B = std::max(64UL, (((ram_use * 8UL) >> 6) << 6));
  std::uint64_t n_blocks_B = (2UL * text_length + max_block_size_B - 1) / max_block_size_B;
  long double text_to_ram_ratio = (long double)text_length / (long double)ram_use;

  // Print info about halfsegments.
  fprintf(stderr, "  Max halfsegment size = %lu (%.2LfMiB)\n", max_halfsegment_size, (1.L * max_halfsegment_size / (1UL << 20)));
  fprintf(stderr, "  Number of halfsegments = %lu\n", n_halfsegments);
  fprintf(stderr, "  Number of halfsegment pairs = %lu\n", n_different_halfsegment_pairs);

  // Initialize file names with halfsegment pairs.
  std::string **pairs_filenames = new std::string*[n_halfsegments];
  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    pairs_filenames[i] = new std::string[n_halfsegments];
    for (std::uint64_t j = i; j < n_halfsegments; ++j) {
      std::string filename = B_filename + ".pairs." + utils::intToStr(i) + "_" + utils::intToStr(j);
      pairs_filenames[i][j] = filename;
    }
  }

  // Distribute pairs (i, Phi[i]) such that PLCP[i] is irreducible
  // into files corresponding to different halfsegment pairs and
  // compute the C bitvector.
  std::string C_filename = B_filename + ".irreducible_positions_bv";
  std::uint64_t phi_undefined_position = 0;
  if (text_to_ram_ratio > 8.0L) {
    // Distribute pairs.
    phi_undefined_position = distribute_pairs<text_offset_type>(text_length, max_halfsegment_size,
        ram_use, sa_filename, bwt_filename, pairs_filenames, n_irreducible_lcps, io_volume);

    // Compute C.
    compute_C<text_offset_type>(text_length, max_halfsegment_size, ram_use, phi_undefined_position,
        pairs_filenames, sa_filename, bwt_filename, C_filename, io_volume);
  } else {
    // Distribute pairs and compute C.
    phi_undefined_position = distribute_pairs_and_compute_C<text_offset_type>(text_length,
        max_halfsegment_size, ram_use, sa_filename, bwt_filename, C_filename,
        pairs_filenames, n_irreducible_lcps, io_volume);
  }

  std::string *irreducible_bits_filenames = new std::string[n_blocks_B];
  for (std::uint64_t block_id = 0; block_id < n_blocks_B; ++block_id) {
    std::string filename = B_filename + ".irreducible_bits_bv." + utils::intToStr(block_id);
    irreducible_bits_filenames[block_id] = filename;
  }

  // Process all pairs of halfsegments.
  sum_irreducible_lcps = process_halfsegment_pairs<text_offset_type, ext_text_offset_type>(text_filename,
      text_length, max_block_size_B, max_halfsegment_size, max_overflow_size,
      pairs_filenames, irreducible_bits_filenames, io_volume);

  // Clean up.
  for (std::uint64_t halfseg_id = 0; halfseg_id < n_halfsegments; ++halfseg_id)
    delete[] pairs_filenames[halfseg_id];
  delete[] pairs_filenames;

  // Compute B.
  compute_B<ext_text_offset_type>(text_length, max_block_size_B, phi_undefined_position,
      B_filename, C_filename, irreducible_bits_filenames, io_volume);

  // Update I/O volume.
  total_io_volume += io_volume;

  // Clean up.
  delete[] irreducible_bits_filenames;

  // Print summary.
  long double compute_plcp_bitvector_time = utils::wclock() - compute_plcp_bitvector_start;
  fprintf(stderr, "Summary: time = %.2Lfs, total I/O vol = %.2Lfn\n\n",
      compute_plcp_bitvector_time, (1.L * io_volume) / text_length);
}

// A version, that returns a pointer to B bitvector. Requires at least 2n bits of RAM.
template<typename text_offset_type, typename ext_text_offset_type>
std::uint64_t* compute_plcp_bitvector_large_ram(std::uint64_t text_length, std::uint64_t ram_use,
    std::string text_filename, std::string sa_filename, std::string bwt_filename,
    std::string output_filename, std::uint64_t &n_irreducible_lcps,
    std::uint64_t &sum_irreducible_lcps, std::uint64_t &total_io_volume) {
  fprintf(stderr, "Compute PLCP bitvector (dest = RAM):\n");
  long double compute_plcp_bitvector_start = utils::wclock();

  // Initialize basic parameters.
  long double ram_to_text_ratio = (long double)ram_use / (long double)text_length;
  std::uint64_t *B = NULL;
  std::uint64_t io_volume = 0;

  if (ram_to_text_ratio < 1.375L) {
    // Initialize basic parameters.
    static const std::uint64_t max_overflow_size = (1UL << 20);
    std::uint64_t max_halfsegment_size = std::max(1UL, ram_use / 2);
    std::uint64_t n_halfsegments = (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
    std::uint64_t n_different_halfsegment_pairs = (n_halfsegments * (n_halfsegments + 1)) / 2;

    // Print info about halfsegments.
    fprintf(stderr, "  Max halfsegment size = %lu (%.2LfMiB)\n", max_halfsegment_size, (1.L * max_halfsegment_size / (1UL << 20)));
    fprintf(stderr, "  Number of halfsegments = %lu\n", n_halfsegments);
    fprintf(stderr, "  Number of halfsegment pairs = %lu\n", n_different_halfsegment_pairs);

    // Initialize file names with halfsegment pairs.
    std::string **pairs_filenames = new std::string*[n_halfsegments];
    for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
      pairs_filenames[i] = new std::string[n_halfsegments];
      for (std::uint64_t j = i; j < n_halfsegments; ++j) {
        std::string filename = output_filename + ".pairs." + utils::intToStr(i) + "_" + utils::intToStr(j);
        pairs_filenames[i][j] = filename;
      }
    }

    // Distribute pairs (i, Phi[i]) such that PLCP[i] is irreducible
    // into files corresponding to different halfsegment pairs and
    // compute the C bitvector.
    std::string C_filename = output_filename + ".irreducible_positions_bv";
    std::uint64_t phi_undefined_position = distribute_pairs_and_compute_C<text_offset_type>(text_length,
        max_halfsegment_size, ram_use, sa_filename, bwt_filename, C_filename, pairs_filenames,
        n_irreducible_lcps, io_volume);

    // Process all pairs of halfsegments.
    std::string irreducible_bits_filename = output_filename + ".irreducible_bits";
    sum_irreducible_lcps = process_halfsegment_pairs<text_offset_type, ext_text_offset_type>(text_filename,
        text_length, max_halfsegment_size, max_overflow_size, pairs_filenames,
        irreducible_bits_filename, io_volume);

    // Clean up.
    for (std::uint64_t halfseg_id = 0; halfseg_id < n_halfsegments; ++halfseg_id)
      delete[] pairs_filenames[halfseg_id];
    delete[] pairs_filenames;

    // Allocate B.
    std::uint64_t B_size_in_words = (2UL * text_length + 63) / 64;
    B = new std::uint64_t[B_size_in_words];
    std::fill(B, B + B_size_in_words, 0UL);

    // Compute B.
    compute_B<ext_text_offset_type>(text_length, B, irreducible_bits_filename,
        C_filename, phi_undefined_position, io_volume);
  } else {
    // Compute B.
    B = compute_B<text_offset_type>(text_length, text_filename, sa_filename,
        n_irreducible_lcps, sum_irreducible_lcps, io_volume);
  }

  // Update I/O volume.
  total_io_volume += io_volume;

  // Print summary.
  long double compute_plcp_bitvector_time = utils::wclock() - compute_plcp_bitvector_start;
  fprintf(stderr, "Summary: time = %.2Lfs, total I/O vol = %.2Lfn\n\n",
      compute_plcp_bitvector_time, (1.L * io_volume) / text_length);

  // Return pointer to B.
  return B;
}

template<typename text_offset_type, typename ext_text_offset_type>
void compute_plcp_bitvector(std::uint64_t text_length, std::uint64_t ram_use,
    std::string text_filename, std::string sa_filename, std::string bwt_filename,
    std::string output_filename, std::uint64_t &n_irreducible_lcps,
    std::uint64_t &sum_irreducible_lcps, std::uint64_t &total_io_volume) {
  long double text_to_ram_ratio = (long double)text_length / (long double)ram_use;
  if (text_to_ram_ratio > 4.0L) {
    // Not enough RAM to hold B in RAM.
    compute_plcp_bitvector_small_ram<text_offset_type, ext_text_offset_type>(text_length, ram_use, text_filename,
        sa_filename, bwt_filename, output_filename, n_irreducible_lcps, sum_irreducible_lcps, total_io_volume);
  } else {
    // Enough RAM to hold B in RAM.
    std::uint64_t *B = compute_plcp_bitvector_large_ram<text_offset_type, ext_text_offset_type>(text_length,
        ram_use, text_filename, sa_filename, bwt_filename, output_filename, n_irreducible_lcps,
        sum_irreducible_lcps, total_io_volume);

    // Write B to disk.
    {
      // Start the timer.
      fprintf(stderr, "Write PLCP bitvector to disk: ");
      long double write_plcp_start = utils::wclock();
      std::uint64_t io_volume = 0;

      // Write the data.
      std::uint64_t length_of_B_in_words = (2UL * text_length + 63) / 64;
      utils::write_to_file(B, length_of_B_in_words, output_filename);

      // Update I/O volume.
      io_volume += length_of_B_in_words * sizeof(std::uint64_t);
      total_io_volume += io_volume;

      // Print summary.
      long double write_plcp_time = utils::wclock() - write_plcp_start;
      fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, I/O vol = %.2Lfn\n\n", write_plcp_time,
          ((1.L * io_volume) / (1L << 20)) / write_plcp_time, (1.L * io_volume) / text_length);
    }
    delete[] B;
  }
}

template<typename text_offset_type, typename ext_text_offset_type>
void compute_plcp_bitvector(std::string text_filename, std::string sa_filename,
    std::string bwt_filename, std::string output_filename, std::uint64_t ram_use) {
  utils::drop_disk_pages(text_filename);
  utils::drop_disk_pages(sa_filename);
  utils::drop_disk_pages(bwt_filename);
  srand(time(0) + getpid());
  long double global_start = utils::wclock();
  std::uint64_t total_io_volume = 0;

  // Compute basic parameters.
  std::uint64_t text_length = utils::file_size(text_filename);
  std::uint64_t n_irreducible_lcps = 0;
  std::uint64_t sum_irreducible_lcps = 0;
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

  compute_plcp_bitvector<text_offset_type, ext_text_offset_type>(text_length, ram_use, text_filename,
      sa_filename, bwt_filename, output_filename, n_irreducible_lcps,
      sum_irreducible_lcps, total_io_volume);

  // Print summary.
  long double total_time = utils::wclock() - global_start;
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  elapsed time = %.2Lfs (%.3Lfs/MiB of text)\n", total_time, total_time / (1.L * text_length / (1L << 20)));
  fprintf(stderr, "  speed = %.2LfMiB of text/s\n", (1.L * text_length / (1L << 20)) / total_time);
  fprintf(stderr, "  I/O volume = %lu (%.2Lfbytes/input symbol)\n", total_io_volume, (1.L * total_io_volume) / text_length);
  fprintf(stderr, "  number of irreducible LCPs = %lu\n", n_irreducible_lcps);
  fprintf(stderr, "  sum of irreducible LCPs = %lu\n", sum_irreducible_lcps);
}

}  // namespace em_succinct_irreducible_private

#endif  // __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_PLCP_BITVECTOR_HPP_INCLUDED
