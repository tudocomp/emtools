/**
 * @file    em_succinct_irreducible_src/distribute_pairs_and_compute_C.hpp
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

#ifndef __EM_SUCCINCT_IRREDUCIBLE_SRC_DISTRIBUTE_PAIRS_AND_COMPUTE_C_HPP_INCLUDED
#define __EM_SUCCINCT_IRREDUCIBLE_SRC_DISTRIBUTE_PAIRS_AND_COMPUTE_C_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <string>
#include <algorithm>
#include <omp.h>

#include "io/async_stream_reader.hpp"
#include "io/async_multi_stream_writer.hpp"
#include "set_bits.hpp"
#include "utils.hpp"


namespace em_succinct_irreducible_private {

template<typename text_offset_type>
std::uint64_t distribute_pairs(std::uint64_t text_length, std::uint64_t max_halfsegment_size,
    std::uint64_t ram_use, std::string sa_filename, std::string bwt_filename, std::string **pairs_filenames,
    std::uint64_t &n_irreducible_lcps, std::uint64_t &total_io_volume) {
  std::uint64_t n_halfsegments = (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t n_irreducible = 0;
  std::uint64_t phi_undefined_position = 0;

  fprintf(stderr, "  Distribute irreducible (i, Phi[i]) pairs: ");
  long double start = utils::wclock();

  // Create a map from used halfsegment pairs to a contiguous
  // range of integers. This is needed to use multifile writer.
  std::uint64_t **halfseg_ids_to_file_id = new std::uint64_t*[n_halfsegments];
  {
    for (std::uint64_t i = 0; i < n_halfsegments; ++i)
      halfseg_ids_to_file_id[i] = new std::uint64_t[n_halfsegments];

    std::uint64_t file_counter = 0;
    for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
      for (std::uint64_t j = i; j < n_halfsegments; ++j) {
        halfseg_ids_to_file_id[i][j] = file_counter;
        halfseg_ids_to_file_id[j][i] = file_counter;
        ++file_counter;
      }
    }
  }

  // Initialize multifile writer of (i, Phi[i]) pairs.
  static const std::uint64_t n_free_buffers = 4;
  std::uint64_t halfseg_buffers_ram = ram_use;
  std::uint64_t n_different_halfseg_pairs = (n_halfsegments * (n_halfsegments + 1)) / 2;
  std::uint64_t buffer_size = std::max(1UL, halfseg_buffers_ram / (n_different_halfseg_pairs + n_free_buffers));
  typedef async_multi_stream_writer<text_offset_type> pair_multiwriter_type;
  pair_multiwriter_type *pair_multiwriter = new pair_multiwriter_type(buffer_size, n_free_buffers);
  for (std::uint64_t i = 0; i < n_halfsegments; ++i)
    for (std::uint64_t j = i; j < n_halfsegments; ++j)
      pair_multiwriter->add_file(pairs_filenames[i][j]);

  // Initialize suffix array reader.
  typedef async_stream_reader<text_offset_type> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(sa_filename);

  // Initialize BWT reader.
  typedef async_stream_reader<std::uint8_t> bwt_reader_type;
  bwt_reader_type *bwt_reader = new bwt_reader_type(bwt_filename);

  // Distribution follows.
  std::uint8_t prev_bwt = 0;
  std::uint64_t prev_sa = 0;
  std::uint64_t prev_halfseg_id = 0;
  for (std::uint64_t i = 0; i < text_length; ++i) {
    std::uint64_t cur_sa = sa_reader->read();
    std::uint64_t cur_halfseg_id = cur_sa / max_halfsegment_size;
    std::uint8_t cur_bwt = bwt_reader->read();

    if (i == 0 || cur_sa == 0 || prev_sa == 0 || cur_bwt != prev_bwt) {
      // PLCP[cur_sa] is irreducible. Write (i, Phi[i]) to appropriate file.
      ++n_irreducible;
      if (i > 0) {
        std::uint64_t file_id = halfseg_ids_to_file_id[cur_halfseg_id][prev_halfseg_id];
        pair_multiwriter->write_to_ith_file(file_id, (text_offset_type)cur_sa);
        pair_multiwriter->write_to_ith_file(file_id, (text_offset_type)prev_sa);
      } else phi_undefined_position = cur_sa;
    }

    prev_halfseg_id = cur_halfseg_id;
    prev_sa = cur_sa;
    prev_bwt = cur_bwt;
  }

  // Print summary.
  long double elapsed = utils::wclock() - start;
  std::uint64_t io_volume = sa_reader->bytes_read() + bwt_reader->bytes_read() + pair_multiwriter->bytes_written();
  total_io_volume += io_volume;
  fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
      elapsed, ((1.L * io_volume) / (1L << 20)) / elapsed, (1.L * total_io_volume) / text_length);

  // Clean up.
  delete bwt_reader;
  delete sa_reader;
  delete pair_multiwriter;
  for (std::uint64_t i = 0; i < n_halfsegments; ++i)
    delete[] halfseg_ids_to_file_id[i];
  delete[] halfseg_ids_to_file_id;

  // Return undefined Phi position.
  n_irreducible_lcps = n_irreducible;
  return phi_undefined_position;
}

template<typename text_offset_type>
void compute_C(std::uint64_t text_length, std::uint64_t max_halfsegment_size, std::uint64_t ram_use,
    std::uint64_t phi_undefined_position, std::string **pairs_filenames, std::string sa_filename,
    std::string bwt_filename, std::string C_filename, std::uint64_t &total_io_volume) {
  std::uint64_t n_halfsegments = (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t max_block_size = 8UL * ram_use;
  while (max_block_size & 63UL)
    ++max_block_size;

  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;
  std::uint64_t io_vol_scan_sa = (1 + sizeof(text_offset_type)) * text_length * n_blocks;

  std::uint64_t io_vol_scan_pairs = 0;
  for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
    std::uint64_t block_beg = block_id * max_block_size;
    std::uint64_t block_end = std::min(block_beg + max_block_size, text_length);
    for (std::uint64_t left_halfseg_id = 0; left_halfseg_id < n_halfsegments; ++left_halfseg_id) {
      std::uint64_t left_halfseg_beg = left_halfseg_id * max_halfsegment_size;
      std::uint64_t left_halfseg_end = std::min(left_halfseg_beg + max_halfsegment_size, text_length);
      for (std::uint64_t right_halfseg_id = left_halfseg_id; right_halfseg_id < n_halfsegments; ++right_halfseg_id) {
        std::uint64_t right_halfseg_beg = right_halfseg_id * max_halfsegment_size;
        std::uint64_t right_halfseg_end = std::min(right_halfseg_beg + max_halfsegment_size, text_length);
        if ((left_halfseg_end > block_beg && block_end > left_halfseg_beg) ||
            (right_halfseg_end > block_beg && block_end > right_halfseg_beg))
          io_vol_scan_pairs += utils::file_size(pairs_filenames[left_halfseg_id][right_halfseg_id]);
      }
    }
  }

  if (io_vol_scan_sa <= io_vol_scan_pairs) {
    fprintf(stderr, "  Compute bitvector C (method I): ");
    long double start = utils::wclock();
    std::uint64_t io_vol = 0;

    // Allocate the array holding the block of C.
    std::uint64_t max_block_size_in_words = max_block_size / 64;
    std::uint64_t *C = new std::uint64_t[max_block_size_in_words];
    std::FILE *f = utils::file_open(C_filename, "w");

    // Initialize the buffer.
    static const std::uint64_t buffer_size = (1UL << 20);
    std::uint64_t *buf = new std::uint64_t[buffer_size];
#ifdef _OPENMP
    std::uint64_t *tempbuf = new std::uint64_t[buffer_size];
#endif

    for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
      std::uint64_t block_beg = block_id * max_block_size;
      std::uint64_t block_end = std::min(block_beg + max_block_size, text_length);
      std::uint64_t block_size = block_end - block_beg;
      std::uint64_t block_size_in_words = (block_size + 63) / 64;

      // Zero-initialize the block of C.
      std::fill(C, C + block_size_in_words, 0UL);

      // Initialize suffix array reader.
      typedef async_stream_reader<text_offset_type> sa_reader_type;
      sa_reader_type *sa_reader = new sa_reader_type(sa_filename);

      // Initialize BWT reader.
      typedef async_stream_reader<std::uint8_t> bwt_reader_type;
      bwt_reader_type *bwt_reader = new bwt_reader_type(bwt_filename);

      // Scan SA and BWT left to right.
      std::uint64_t filled = 0;
      std::uint8_t prev_bwt = 0;
      std::uint64_t prev_sa = 0;
      for (std::uint64_t i = 0; i < text_length; ++i) {
        std::uint64_t cur_sa = sa_reader->read();
        std::uint8_t cur_bwt = bwt_reader->read();

        if (block_beg <= cur_sa && cur_sa < block_end &&
            (i == 0 || cur_sa == 0 || prev_sa == 0 || cur_bwt != prev_bwt)) {
          // PLCP[cur_sa] is irreducible.
          std::uint64_t offset = cur_sa - block_beg;
          buf[filled++] = offset;
          if (filled == buffer_size) {
#ifdef _OPENMP
            set_bits(C, block_size, buf, filled, tempbuf);
#else
            set_bits(C, buf, filled);
#endif
            filled = 0;
          }
        }

        prev_sa = cur_sa;
        prev_bwt = cur_bwt;
      }

      // Flush the remaining items in the buffer.
      if (filled > 0) {
#ifdef _OPENMP
        set_bits(C, block_size, buf, filled, tempbuf);
#else
        set_bits(C, buf, filled);
#endif
        filled = 0;
      }

      // Write current block of C to file.
      utils::write_to_file(C, block_size_in_words, f);

      // Update I/O volume.
      io_vol += sa_reader->bytes_read() + bwt_reader->bytes_read() + block_size_in_words * sizeof(std::uint64_t);

      // Clean up.
      delete sa_reader;
      delete bwt_reader;
    }

    // Clean up.
#ifdef _OPENMP
    delete[] tempbuf;
#endif
    delete[] buf;
    delete[] C;
    std::fclose(f);

    // Update I/O volume.
    total_io_volume += io_vol;

    // Print summary.
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n", elapsed,
        ((1.L * io_vol) / (1L << 20)) / elapsed, (1.L * total_io_volume) / text_length);
  } else {
    fprintf(stderr, "  Compute bitvector C (method II): ");
    long double start = utils::wclock();
    std::uint64_t io_vol = 0;

    // Allocate the array holding the block of C.
    std::uint64_t max_block_size_in_words = max_block_size / 64;
    std::uint64_t *C = new std::uint64_t[max_block_size_in_words];
    std::FILE *f = utils::file_open(C_filename, "w");

    // Initialize the buffer.
    static const std::uint64_t buffer_size = (1UL << 20);
    std::uint64_t *buf = new std::uint64_t[buffer_size];
#ifdef _OPENMP
    std::uint64_t *tempbuf = new std::uint64_t[buffer_size];
#endif

    // Process blocks of C left to right.
    for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
      std::uint64_t block_beg = block_id * max_block_size;
      std::uint64_t block_end = std::min(block_beg + max_block_size, text_length);
      std::uint64_t block_size = block_end - block_beg;
      std::uint64_t block_size_in_words = (block_size + 63) / 64;

      // Zero-initialize the block of C.
      std::fill(C, C + block_size_in_words, 0UL);

      // Iterate through all pairs of halfsegments.
      std::uint64_t filled = 0;
      for (std::uint64_t left_halfseg_id = 0; left_halfseg_id < n_halfsegments; ++left_halfseg_id) {
        std::uint64_t left_halfseg_beg = left_halfseg_id * max_halfsegment_size;
        std::uint64_t left_halfseg_end = std::min(left_halfseg_beg + max_halfsegment_size, text_length);

        for (std::uint64_t right_halfseg_id = left_halfseg_id; right_halfseg_id < n_halfsegments; ++right_halfseg_id) {
          std::uint64_t right_halfseg_beg = right_halfseg_id * max_halfsegment_size;
          std::uint64_t right_halfseg_end = std::min(right_halfseg_beg + max_halfsegment_size, text_length);

          if ((left_halfseg_end > block_beg && block_end > left_halfseg_beg) ||
              (right_halfseg_end > block_beg && block_end > right_halfseg_beg)) {
            // Initialize reading of pairs.
            typedef async_stream_reader<text_offset_type> pair_reader_type;
            pair_reader_type *pair_reader = new pair_reader_type(pairs_filenames[left_halfseg_id][right_halfseg_id]);

            while (pair_reader->empty() == false) {
              std::uint64_t i = pair_reader->read();
              pair_reader->read();  // Skip Phi[i].

              if (block_beg <= i && i < block_end) {
                std::uint64_t offset = i - block_beg;
                buf[filled++] = offset;
                if (filled == buffer_size) {
#ifdef _OPENMP
                  set_bits(C, block_size, buf, filled, tempbuf);
#else
                  set_bits(C, buf, filled);
#endif
                  filled = 0;
                }
              }
            }

            // Update I/O volume.
            io_vol += pair_reader->bytes_read();

            // Clean up.
            delete pair_reader;
          }
        }
      }

      // Flush the remaining items in the buffer.
      if (filled > 0) {
#ifdef _OPENMP
        set_bits(C, block_size, buf, filled, tempbuf);
#else
        set_bits(C, buf, filled);
#endif
        filled = 0;
      }

      // Special case.
      if (block_beg <= phi_undefined_position && phi_undefined_position < block_end) {
        std::uint64_t offset = phi_undefined_position - block_beg;
        C[offset >> 6] |= (1UL << (offset & 63));
      }

      // Write current block of C to file.
      utils::write_to_file(C, block_size_in_words, f);

      // Update I/O volume.
      io_vol += block_size_in_words * sizeof(std::uint64_t);
    }

    // Clean up.
#ifdef _OPENMP
    delete[] tempbuf;
#endif
    delete[] buf;
    delete[] C;
    std::fclose(f);


    // Update I/O volume.
    total_io_volume += io_vol;

    // Print summary.
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n", elapsed,
        ((1.L * io_vol) / (1L << 20)) / elapsed, (1.L * total_io_volume) / text_length);
  }
}

template<typename text_offset_type>
std::uint64_t distribute_pairs_and_compute_C(std::uint64_t text_length,
    std::uint64_t max_halfsegment_size, std::uint64_t ram_use, std::string sa_filename,
    std::string bwt_filename, std::string C_filename, std::string **pairs_filenames,
    std::uint64_t &n_irreducible_lcps, std::uint64_t &total_io_volume) {
  fprintf(stderr, "  Distribute irreducible (i, Phi[i]) pairs and compute bitvector C: ");
  long double start = utils::wclock();

  // Initialize basic parameters.
  std::uint64_t n_halfsegments = (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t n_different_halfseg_pairs = (n_halfsegments * (n_halfsegments + 1)) / 2;
  std::uint64_t io_volume = 0;
  std::uint64_t n_irreducible = 0;
  std::uint64_t phi_undefined_position = 0;

  // Allocate bitvector C.
  std::uint64_t C_size_in_words = (text_length + 63) / 64;
  std::uint64_t C_size_in_bytes = (text_length + 7) / 8;
  std::uint64_t *C = new std::uint64_t[C_size_in_words];
  std::fill(C, C + C_size_in_words, 0UL);

  // Create a map from used halfsegment pairs to a contiguous
  // range of integers. This is needed to use multifile writer.
  std::uint64_t **halfseg_ids_to_file_id = new std::uint64_t*[n_halfsegments];
  {
    for (std::uint64_t i = 0; i < n_halfsegments; ++i)
      halfseg_ids_to_file_id[i] = new std::uint64_t[n_halfsegments];

    std::uint64_t file_counter = 0;
    for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
      for (std::uint64_t j = i; j < n_halfsegments; ++j) {
        halfseg_ids_to_file_id[i][j] = file_counter;
        halfseg_ids_to_file_id[j][i] = file_counter;
        ++file_counter;
      }
    }
  }

  // Initialize multifile writer of (i, Phi[i]) pairs.
  static const std::uint64_t n_free_buffers = 4;
  std::uint64_t halfseg_buffers_ram = ram_use - C_size_in_bytes;
  std::uint64_t buffer_size = std::max((1UL << 20), halfseg_buffers_ram / (n_different_halfseg_pairs + n_free_buffers));
  typedef async_multi_stream_writer<text_offset_type> pair_multiwriter_type;
  pair_multiwriter_type *pair_multiwriter = new pair_multiwriter_type(buffer_size, n_free_buffers);
  for (std::uint64_t i = 0; i < n_halfsegments; ++i)
    for (std::uint64_t j = i; j < n_halfsegments; ++j)
      pair_multiwriter->add_file(pairs_filenames[i][j]);

  // Initialize suffix array reader.
  typedef async_stream_reader<text_offset_type> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(sa_filename);

  // Initialize BWT reader.
  typedef async_stream_reader<std::uint8_t> bwt_reader_type;
  bwt_reader_type *bwt_reader = new bwt_reader_type(bwt_filename);

  // Initialize the buffer.
  static const std::uint64_t local_buffer_size = (1UL << 20);
  std::uint64_t *buf = new std::uint64_t[local_buffer_size];
#ifdef _OPENMP
  std::uint64_t *tempbuf = new std::uint64_t[local_buffer_size];
#endif

  // Distribution follows.
  std::uint64_t filled = 0;
  std::uint8_t prev_bwt = 0;
  std::uint64_t prev_sa = 0;
  std::uint64_t prev_halfseg_id = 0;
  for (std::uint64_t i = 0; i < text_length; ++i) {
    std::uint64_t cur_sa = sa_reader->read();
    std::uint64_t cur_halfseg_id = cur_sa / max_halfsegment_size;
    std::uint8_t cur_bwt = bwt_reader->read();

    if (i == 0 || cur_sa == 0 || prev_sa == 0 || cur_bwt != prev_bwt) {
      // PLCP[cur_sa] is irreducible. Write (i, Phi[i]) to appropriate file.
      ++n_irreducible;
      buf[filled++] = cur_sa;
      if (filled == local_buffer_size) {
#ifdef _OPENMP
        set_bits(C, text_length, buf, filled, tempbuf);
#else
        set_bits(C, buf, filled);
#endif
        filled = 0;
      }

      if (i > 0) {
        std::uint64_t file_id = halfseg_ids_to_file_id[cur_halfseg_id][prev_halfseg_id];
        pair_multiwriter->write_to_ith_file(file_id, (text_offset_type)cur_sa);
        pair_multiwriter->write_to_ith_file(file_id, (text_offset_type)prev_sa);
      } else phi_undefined_position = cur_sa;
    }

    prev_halfseg_id = cur_halfseg_id;
    prev_sa = cur_sa;
    prev_bwt = cur_bwt;
  }

  // Flush the remaining items in the buffer.
  if (filled > 0) {
#ifdef _OPENMP
    set_bits(C, text_length, buf, filled, tempbuf);
#else
    set_bits(C, buf, filled);
#endif
    filled = 0;
  }

  // Write C to disk.
  utils::write_to_file(C, C_size_in_words, C_filename);

  // Update I/O volume.
  io_volume += sa_reader->bytes_read() + bwt_reader->bytes_read() +
      pair_multiwriter->bytes_written() + C_size_in_words * sizeof(std::uint64_t);
  total_io_volume += io_volume;

  // Clean up.
#ifdef _OPENMP
  delete[] tempbuf;
#endif
  delete[] buf;
  delete[] C;
  delete bwt_reader;
  delete sa_reader;
  delete pair_multiwriter;
  for (std::uint64_t i = 0; i < n_halfsegments; ++i)
    delete[] halfseg_ids_to_file_id[i];
  delete[] halfseg_ids_to_file_id;

  // Print summary.
  long double elapsed = utils::wclock() - start;
  fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
      elapsed, ((1.L * io_volume) / (1L << 20)) / elapsed, (1.L * total_io_volume) / text_length);

  // Update reference variables.
  n_irreducible_lcps = n_irreducible;

  // Return undefined Phi position.
  return phi_undefined_position;
}

}  // namespace em_succinct_irreducible_private

#endif  // __EM_SUCCINCT_IRREDUCIBLE_SRC_DISTRIBUTE_PAIRS_AND_COMPUTE_C_HPP_INCLUDED
