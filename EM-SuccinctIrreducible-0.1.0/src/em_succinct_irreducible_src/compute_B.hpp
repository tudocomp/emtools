/**
 * @file    em_succinct_irreducible_src/compute_B.hpp
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

#ifndef __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_B_HPP_INCLUDED
#define __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_B_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "io/async_stream_reader.hpp"

#include "set_bits.hpp"
#include "utils.hpp"


namespace em_succinct_irreducible_private {

template<typename ext_text_offset_type>
void compute_B(std::uint64_t text_length, std::uint64_t *B,
    std::string irreducible_bits_filename, std::string C_filename,
    std::uint64_t phi_undefined_position, std::uint64_t &total_io_volume) {
  fprintf(stderr, "  Compute bitvector encoding of PLCP array: ");
  long double start = utils::wclock();
  std::uint64_t io_volume = 0;

  // Fill in the bits in B corresponding
  // to irreducible lcp values. 
  {
    // Initialize reader of irreducible positions.
    typedef async_stream_reader<ext_text_offset_type> irreducible_bits_reader_type;
    irreducible_bits_reader_type *irreducible_bits_reader =
      new irreducible_bits_reader_type(irreducible_bits_filename);

    // Allocate the buffer.
    static const std::uint64_t buffer_size = (1UL << 20);
    ext_text_offset_type *buf = new ext_text_offset_type[buffer_size];
#ifdef _OPENMP
    ext_text_offset_type *tempbuf = new ext_text_offset_type[buffer_size];
#endif

    // Stream and set bits inside B.
    std::uint64_t count = utils::file_size(irreducible_bits_filename) / sizeof(ext_text_offset_type);
    {
      std::uint64_t items_processed = 0;
      while (items_processed < count) {
        std::uint64_t filled = std::min(count - items_processed, buffer_size);
        irreducible_bits_reader->read(buf, filled);
#ifdef _OPENMP
        set_bits(B, 2UL * text_length, buf, filled, tempbuf);
#else
        for (std::uint64_t j = 0; j < filled; ++j) {
          std::uint64_t idx = buf[j];
          B[idx >> 6] |= (1UL << (idx & 63));
        }
#endif
        items_processed += filled;
      }
    }

    // Special case.
    {
      std::uint64_t idx = 2 * phi_undefined_position;
      B[idx >> 6] |= (1UL << (idx & 63));
    }

    // Update I/O volume.
    io_volume += irreducible_bits_reader->bytes_read();

    // Clean up.
#ifdef _OPENMP
    delete[] tempbuf;
#endif
    delete[] buf;
    delete irreducible_bits_reader;
    utils::file_delete(irreducible_bits_filename);
  }

  // Fill in reducible LCP values.
  {
    // Initialize reader of C.
    typedef async_stream_reader<std::uint64_t> C_reader_type;
    C_reader_type *C_reader = new C_reader_type(C_filename);

    // Initialize the bit-buffer for reader of C.
    std::uint64_t bitbuf = C_reader->read();
    std::uint64_t bitbuf_pos = 0;
    bool C_bit = (bitbuf & (1UL << (bitbuf_pos++)));

    // Add reducible bits.
    for (std::uint64_t j = 0; j < 2UL * text_length; ++j) {
      // Set the bit in B.
      if (C_bit == 0)
        B[j >> 6] |= (1UL << (j & 63));

      // Read the next bit from C.
      if (B[j >> 6] & (1UL << (j & 63))) {
        if (bitbuf_pos < 64 || C_reader->empty() == false) {
          if (bitbuf_pos == 64) {
            bitbuf = C_reader->read();
            bitbuf_pos = 0;
          }
          C_bit = (bitbuf & (1UL << (bitbuf_pos++)));
        }
      }
    }

    // Update I/O volume.
    io_volume += C_reader->bytes_read();

    // Clean up.
    delete C_reader;
    utils::file_delete(C_filename);
  }

  // Update I/O volume.
  total_io_volume += io_volume;

  // Print summary.
  long double elapsed = utils::wclock() - start;
  fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
      elapsed, ((1.L * io_volume) / (1L << 20)) / elapsed, (1.L * total_io_volume) / text_length);
}

template<typename ext_text_offset_type>
void compute_B(std::uint64_t text_length, std::uint64_t max_block_size_B,
    std::uint64_t phi_undefined_position, std::string B_filename, std::string C_filename,
    std::string *irreducible_bits_filenames, std::uint64_t &total_io_volume) {
  std::uint64_t n_blocks_B = (2UL * text_length + max_block_size_B - 1) / max_block_size_B;

  fprintf(stderr, "  Compute bitvector encoding of PLCP array: ");
  long double start = utils::wclock();

  // Initialize reader of C.
  typedef async_stream_reader<std::uint64_t> C_reader_type;
  C_reader_type *C_reader = new C_reader_type(C_filename);

  // Initialize the bit-buffer for reader of C.
  std::uint64_t bitbuf = C_reader->read();
  std::uint64_t bitbuf_pos = 0;
  bool C_bit = (bitbuf & (1UL << (bitbuf_pos++)));

  std::uint64_t io_vol = 0;
  std::uint64_t max_block_size_B_in_words = max_block_size_B / 64;
  std::uint64_t *B = new std::uint64_t[max_block_size_B_in_words];
  std::FILE *f = utils::file_open(B_filename, "w");

  // Allocate the buffer.
  static const std::uint64_t buffer_size = (1UL << 20);
  ext_text_offset_type *buf = new ext_text_offset_type[buffer_size];
#ifdef _OPENMP
  ext_text_offset_type *tempbuf = new ext_text_offset_type[buffer_size];
#endif

  for (std::uint64_t block_id = 0; block_id < n_blocks_B; ++block_id) {
    std::uint64_t block_beg = block_id * max_block_size_B;
    std::uint64_t block_end = std::min(block_beg + max_block_size_B, 2 * text_length);
    std::uint64_t block_size = block_end - block_beg;
    std::uint64_t block_size_in_words = (block_size + 63) / 64;

    // Zero-initialize the block of B.
    std::fill(B, B + block_size_in_words, 0UL);

    // Initialize the reader of irreducible positions.
    typedef async_stream_reader<ext_text_offset_type> irreducible_bits_reader_type;
    irreducible_bits_reader_type *irreducible_bits_reader =
      new irreducible_bits_reader_type(irreducible_bits_filenames[block_id]);

    // Read and set the bits in the block of B.
    std::uint64_t count = utils::file_size(irreducible_bits_filenames[block_id]) / sizeof(ext_text_offset_type);
    {
      std::uint64_t items_processed = 0;
      while (items_processed < count) {
        std::uint64_t filled = std::min(count - items_processed, buffer_size);
        irreducible_bits_reader->read(buf, filled);
#ifdef _OPENMP
        #pragma omp parallel for
        for (std::uint64_t j = 0; j < filled; ++j)
          buf[j] = (std::uint64_t)buf[j] - block_beg;

        set_bits(B, block_size, buf, filled, tempbuf);
#else
        for (std::uint64_t j = 0; j < filled; ++j) {
          std::uint64_t idx = buf[j];
          std::uint64_t offset = idx - block_beg;
          B[offset >> 6] |= (1UL << (offset & 63));
        }
#endif

        items_processed += filled;
      }
    }

    // Special case for 1-bit corresponding to PLCP[SA[0]].
    if (block_beg <= 2 * phi_undefined_position && 2 * phi_undefined_position < block_end) {
      std::uint64_t offset = 2 * phi_undefined_position - block_beg;
      B[offset >> 6] |= (1UL << (offset & 63));
    }

    // Add reducible bits.
    for (std::uint64_t j = 0; j < block_size; ++j) {
      // Set the bit in B.
      if (C_bit == 0)
        B[j >> 6] |= (1UL << (j & 63));

      // Read the next bit from C.
      if (B[j >> 6] & (1UL << (j & 63))) {
        if (bitbuf_pos < 64 || C_reader->empty() == false) {
          if (bitbuf_pos == 64) {
            bitbuf = C_reader->read();
            bitbuf_pos = 0;
          }
          C_bit = (bitbuf & (1UL << (bitbuf_pos++)));
        }
      }
    }

    // Write current block of B to file.
    utils::write_to_file(B, block_size_in_words, f);

    // Update I/O volume.
    io_vol += irreducible_bits_reader->bytes_read() + block_size_in_words * sizeof(std::uint64_t);

    // Clean up.
    delete irreducible_bits_reader;
    utils::file_delete(irreducible_bits_filenames[block_id]);
  }

  // Update I/O volume.
  io_vol += C_reader->bytes_read();
  total_io_volume += io_vol;

  // Clean up.
#ifdef _OPENMP
  delete[] tempbuf;
#endif
  delete[] buf;
  delete[] B;
  delete C_reader;
  std::fclose(f);
  utils::file_delete(C_filename);

  // Print summary.
  long double elapsed = utils::wclock() - start;
  fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n", elapsed,
      ((1.L * io_vol) / (1L << 20)) / elapsed, (1.L * total_io_volume) / text_length);
}

template<typename text_offset_type>
std::uint64_t *compute_B(std::uint64_t text_length, std::string text_filename,
    std::string sa_filename, std::uint64_t &n_irreducible_lcps,
    std::uint64_t &sum_irreducible_lcps, std::uint64_t &total_io_volume) {
  // Initialize basic parameters.
  std::uint64_t local_n_irreducible_lcps = 0;
  std::uint64_t local_sum_irreducible_lcps = 0;

  // Allocate bitvectors.
  std::uint64_t C_size_in_words = (text_length + 63) / 64;
  std::uint64_t B_size_in_words = (2UL * text_length + 63) / 64;
  std::uint64_t *C = new std::uint64_t[C_size_in_words];
  std::uint64_t *B = new std::uint64_t[B_size_in_words];
  std::fill(C, C + C_size_in_words, 0UL);
  std::fill(B, B + B_size_in_words, 0UL);

  // Read text.
  std::uint8_t *text = new std::uint8_t[text_length];
  {
    // Start the timer.
    fprintf(stderr, "  Read text: ");
    long double read_start = utils::wclock();
    std::uint64_t io_volume = 0;

    // Read data.
    utils::read_from_file(text, text_length, text_filename);

    // Update I/O volume.
    io_volume += text_length;
    total_io_volume += io_volume;

    // Print summary.
    long double read_time = utils::wclock() - read_start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n", read_time,
        ((1.L * io_volume) / (1L << 20)) / read_time, (1.L * total_io_volume) / text_length);
  }

  // Compute irreducible lcp values.
  {
    // Start the timer.
    fprintf(stderr, "  Compute irreducible LCP values: ");
    long double compute_irr_lcp_start = utils::wclock();

    // Initialize basic statistics.
    std::uint64_t io_volume = 0;

    // Initialize SA reader.
    typedef async_stream_reader<text_offset_type> sa_reader_type;
    sa_reader_type *sa_reader = new sa_reader_type(sa_filename);

    // Allocate buffers.
    static const std::uint64_t buf_size = (1UL << 20);
    text_offset_type *sa_buf = new text_offset_type[buf_size];
    std::uint8_t *bwt_buf = new std::uint8_t[buf_size];
#ifdef _OPENMP
    std::uint64_t *pair_buf = new std::uint64_t[buf_size * 2];
    std::uint64_t *ans_buf_B = new std::uint64_t[buf_size];
    std::uint64_t *ans_buf_C = new std::uint64_t[buf_size];
#endif

    // Processing of SA follows.
    std::uint64_t sa_items_read = 0;
    std::uint64_t prev_sa = text_length;
    std::uint8_t prev_bwt = 0;
    while (sa_items_read < text_length) {
      std::uint64_t buf_filled = std::min(buf_size, text_length - sa_items_read);
      sa_reader->read(sa_buf, buf_filled);

      // Compute BWT buffer.
#ifdef _OPENMP
      #pragma omp parallel for
      for (std::uint64_t j = 0; j < buf_filled; ++j) {
        std::uint64_t addr = (std::uint64_t)sa_buf[j];
        if (addr > 0) bwt_buf[j] = text[addr - 1];
      }
#else
      for (std::uint64_t j = 0; j < buf_filled; ++j) {
        std::uint64_t addr = (std::uint64_t)sa_buf[j];
        if (addr > 0) bwt_buf[j] = text[addr - 1];
      }
#endif

      // Process buffer.
#ifdef _OPENMP
      {
        // Bring the irreducible pairs together.
        std::uint64_t buf_irr_filled = 0;
        for (std::uint64_t j = 0; j < buf_filled; ++j) {
          std::uint64_t cur_sa = (std::uint64_t)sa_buf[j];
          std::uint8_t cur_bwt = bwt_buf[j];
          if ((sa_items_read == 0 && j == 0) || (cur_sa == 0) || (prev_sa == 0) || (cur_bwt != prev_bwt)) {
            pair_buf[2 * buf_irr_filled] = cur_sa;
            pair_buf[2 * buf_irr_filled + 1] = prev_sa;
            ++buf_irr_filled;
          }
          prev_sa = cur_sa;
          prev_bwt = cur_bwt;
        }

        // Update statistics.
        local_n_irreducible_lcps += buf_irr_filled;

        if (buf_irr_filled > 0) {
          // Compute lcp values in parallel.
          #pragma omp parallel
          {
            std::uint64_t thread_sum_irreducible_lcps = 0;

            #pragma omp for nowait
            for (std::uint64_t j = 0; j < buf_irr_filled; ++j) {
              std::uint64_t i = pair_buf[2 * j];
              std::uint64_t phi_i = pair_buf[2 * j + 1];
              std::uint64_t lcp = 0;
              while (i + lcp < text_length && phi_i + lcp < text_length &&
                  text[i + lcp] == text[phi_i + lcp]) ++lcp;
              thread_sum_irreducible_lcps += lcp;
              ans_buf_C[j] = i;
              ans_buf_B[j] = 2 * i + lcp;
            }

            #pragma omp critical
            {
              local_sum_irreducible_lcps += thread_sum_irreducible_lcps;
            }
          }

          // Set the bits in B and C in parallel.
          set_bits(B, 2UL * text_length, ans_buf_B, buf_irr_filled, pair_buf);
          set_bits(C, 1UL * text_length, ans_buf_C, buf_irr_filled, pair_buf);
        }
      }
#else
      for (std::uint64_t j = 0; j < buf_filled; ++j) {
        std::uint64_t cur_sa = (std::uint64_t)sa_buf[j];
        std::uint8_t cur_bwt = bwt_buf[j];
        if ((sa_items_read == 0 && j == 0) || (cur_sa == 0) || (prev_sa == 0) || (cur_bwt != prev_bwt)) {
          // Compute irreducible lcp(cur_sa, prev_sa) naively.
          std::uint64_t lcp = 0;
          while (cur_sa + lcp < text_length && prev_sa + lcp < text_length &&
              text[cur_sa + lcp] == text[prev_sa + lcp]) ++lcp;

          // Set the corresponding bits in the B and C.
          std::uint64_t bv_idx = 2UL * cur_sa + lcp;
          B[bv_idx >> 6] |= (1UL << (bv_idx & 63));
          C[cur_sa >> 6] |= (1UL << (cur_sa & 63));

          // Update statistics.
          ++local_n_irreducible_lcps;
          local_sum_irreducible_lcps += lcp;
        }

        prev_sa = cur_sa;
        prev_bwt = cur_bwt;
      }
#endif

      sa_items_read += buf_filled;
    }

    // Update I/O volume.
    io_volume += sa_reader->bytes_read();
    total_io_volume += io_volume;

    // Clean up.
    delete[] sa_buf;
    delete[] bwt_buf;
    delete sa_reader;
#ifdef _OPENMP
    delete[] pair_buf;
    delete[] ans_buf_B;
    delete[] ans_buf_C;
#endif

    // Print summary.
    long double compute_irr_lcp_time = utils::wclock() - compute_irr_lcp_start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB, total I/O vol = %.2Lfn\n", compute_irr_lcp_time,
        ((1.L * io_volume) / (1L << 20)) / compute_irr_lcp_time, (1.L * total_io_volume) / text_length);
  }

  // Clean up.
  delete[] text;

  // Fill in reducible LCP values.
  {
    fprintf(stderr, "  Fill missing reducible LCP values: ");
    long double fill_in_reduc_start = utils::wclock();

    std::uint64_t B_ptr = 0;
    for (std::uint64_t j = 0; j < text_length; ++j) {
      if ((C[j >> 6] & (1UL << (j & 63))) == 0) {
        // Mark the 1-bit corresponding to reducible LCP value.
        B[B_ptr >> 6] |= (1UL << (B_ptr & 63));
      } else {
        // Find the next 1-bit in B.
        while ((B[B_ptr >> 6] & (1UL << (B_ptr & 63))) == 0)
          ++B_ptr;
      }
      ++B_ptr;
    }

    // Print summary.
    long double fill_in_reduc_time = utils::wclock() - fill_in_reduc_start;
    fprintf(stderr, "time = %.2Lfs\n", fill_in_reduc_time);
  }

  // Clean up.
  delete[] C;

  // Update reference variables.
  n_irreducible_lcps = local_n_irreducible_lcps;
  sum_irreducible_lcps = local_sum_irreducible_lcps;

  // Return the pointer to B.
  return B;
}

}  // namespace em_succinct_irreducible_private

#endif  // __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_B_HPP_INCLUDED
