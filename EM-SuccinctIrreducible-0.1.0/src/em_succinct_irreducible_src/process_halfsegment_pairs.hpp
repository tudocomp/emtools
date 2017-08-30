/**
 * @file    em_succinct_irreducible_src/process_halfsegment_pairs.hpp
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

#ifndef __EM_SUCCINCT_IRREDUCIBLE_SRC_PROCESS_HALFSEGMENT_PAIRS_HPP_INCLUDED
#define __EM_SUCCINCT_IRREDUCIBLE_SRC_PROCESS_HALFSEGMENT_PAIRS_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"
#include "io/async_multi_stream_writer.hpp"
#include "utils.hpp"


namespace em_succinct_irreducible_private {

std::uint64_t naive_lcp(std::uint64_t i, std::uint64_t j, std::uint64_t lcp,
    std::FILE *f_text, std::uint64_t text_length, std::uint64_t &io_volume) {
  static const std::uint64_t bufsize = (1L << 20);
  std::uint8_t *b1 = new std::uint8_t[bufsize];
  std::uint8_t *b2 = new std::uint8_t[bufsize];
  std::uint64_t io_vol = 0;
  while (true) {
    std::uint64_t toread = std::min(bufsize, text_length - std::max(i, j) - lcp);
    if (!toread) break;
    utils::read_at_offset(b1, i + lcp, toread, f_text);
    utils::read_at_offset(b2, j + lcp, toread, f_text);
    io_vol += 2UL * toread;
    std::uint64_t lcp_delta = 0;
    while (lcp_delta < toread && b1[lcp_delta] == b2[lcp_delta])
      ++lcp_delta;
    lcp += lcp_delta;
    if (lcp_delta < toread)
      break;
  }
  delete[] b1;
  delete[] b2;
  io_volume += io_vol;
  return lcp;
}

struct buf_item_ext {
  std::uint64_t m_left_idx;
  std::uint64_t m_right_idx;
  std::uint64_t m_ans;
  std::uint64_t m_block_id;
};

template<typename text_offset_type, typename ext_text_offset_type>
std::uint64_t process_halfsegment_pairs(std::string text_filename,
    std::uint64_t text_length, std::uint64_t max_block_size_B,
    std::uint64_t max_halfsegment_size, std::uint64_t max_overflow_size,
    std::string **pairs_filenames, std::string *irreducible_bits_filenames,
    std::uint64_t &total_io_volume) {
  fprintf(stderr, "  Compute irreducible LCP values:\n");
  long double start = utils::wclock();

  // Initialize basic parameters.
  std::uint64_t n_blocks_B = (2UL * text_length + max_block_size_B - 1) / max_block_size_B;
  std::uint64_t n_halfsegments = (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t sum_irreducible_lcps = 0;

  // Open file with text.
  std::FILE *f_text = utils::file_open(text_filename, "r");

  // Initialize multiwriter of values 2i + PLCP[i].
  typedef async_multi_stream_writer<ext_text_offset_type> lcp_multiwriter_type;
  lcp_multiwriter_type *lcp_multiwriter = NULL;
  {
    static const std::uint64_t n_free_buffers = 4;
    std::uint64_t buffer_size = (1UL << 20);
    lcp_multiwriter = new lcp_multiwriter_type(buffer_size, n_free_buffers);
    for (std::uint64_t block_id = 0; block_id < n_blocks_B; ++block_id)
      lcp_multiwriter->add_file(irreducible_bits_filenames[block_id]);
  }

  // Allocate halfsegments.
  std::uint8_t *left_halfsegment = (std::uint8_t *)malloc(max_halfsegment_size + max_overflow_size);
  std::uint8_t *right_halfsegment = (std::uint8_t *)malloc(max_halfsegment_size + max_overflow_size);

  // Allocate buffers.
  static const std::uint64_t local_buf_size = (1UL << 20);
  text_offset_type *idx_buf = new text_offset_type[local_buf_size * 2];
#ifdef _OPENMP
  buf_item_ext *ans_buf = new buf_item_ext[local_buf_size];
#endif

  // Processing of halfsegment pairs follows.
  for (std::uint64_t left_halfsegment_id = 0; left_halfsegment_id < n_halfsegments; ++left_halfsegment_id) {
    std::uint64_t left_halfsegment_beg = left_halfsegment_id * max_halfsegment_size;
    std::uint64_t left_halfsegment_end = std::min(left_halfsegment_beg + max_halfsegment_size, text_length);
    std::uint64_t left_halfsegment_ext_end = std::min(left_halfsegment_end + max_overflow_size, text_length);
    std::uint64_t left_halfsegment_ext_size = left_halfsegment_ext_end - left_halfsegment_beg;
    bool left_halfsegment_loaded = false;

    // Scan all halfsegments to the right of left_halfsegment_id.
    for (std::uint64_t right_halfsegment_id = left_halfsegment_id; right_halfsegment_id < n_halfsegments; right_halfsegment_id++) {
      std::uint64_t right_halfsegment_beg = right_halfsegment_id * max_halfsegment_size;
      std::uint64_t right_halfsegment_end = std::min(right_halfsegment_beg + max_halfsegment_size, text_length);
      std::uint64_t right_halfsegment_ext_end = std::min(right_halfsegment_end + max_overflow_size, text_length);
      std::uint64_t right_halfsegment_ext_size = right_halfsegment_ext_end - right_halfsegment_beg;

      // Check if that pair of halfsegments has any associated pairs.
      std::string pairs_filename = pairs_filenames[left_halfsegment_id][right_halfsegment_id];
      if (utils::file_exists(pairs_filename) == false || utils::file_size(pairs_filename) == 0) {
        if (utils::file_exists(pairs_filename))
          utils::file_delete(pairs_filename);
        continue;
      }

      // Print initial progress message.
      fprintf(stderr, "    Process halfsegments %lu and %lu: ", left_halfsegment_id, right_halfsegment_id);
      long double halfsegment_process_start = utils::wclock();
      std::uint64_t local_lcp_sum = 0;
      std::uint64_t extra_io = 0;
      std::uint64_t io_vol = 0;

      // Initialize reading from file associated with current pair of halfsegments.
      typedef async_stream_reader<text_offset_type> pair_reader_type;
      std::uint64_t n_pairs = utils::file_size(pairs_filename) / (2 * sizeof(text_offset_type));
      pair_reader_type *pair_reader = new pair_reader_type(pairs_filename);

      // Read left halfsegment from disk (if it wasn't already)
      if (left_halfsegment_loaded == false) {
        utils::read_at_offset(left_halfsegment, left_halfsegment_beg, left_halfsegment_ext_size, text_filename);
        left_halfsegment_loaded = true;
        extra_io += left_halfsegment_ext_size;
      }

      // Read right halfsegment from disk.
      std::uint8_t *right_halfsegment_ptr = right_halfsegment;
      if (right_halfsegment_id != left_halfsegment_id) {
        utils::read_at_offset(right_halfsegment, right_halfsegment_beg, right_halfsegment_ext_size, text_filename);
        extra_io += right_halfsegment_ext_size;
      } else right_halfsegment_ptr = left_halfsegment;

      std::uint64_t pairs_processed = 0;
      while (pairs_processed < n_pairs) {
        std::uint64_t filled = std::min(n_pairs - pairs_processed, local_buf_size);
        pair_reader->read(idx_buf, filled * 2);

#ifdef _OPENMP
        std::vector<std::uint64_t> long_lcps;
        std::uint64_t max_threads = omp_get_max_threads();
        std::uint64_t max_block_size = (filled + max_threads - 1) / max_threads;
        std::uint64_t n_threads = (filled + max_block_size - 1) / max_block_size;
        #pragma omp parallel num_threads(n_threads)
        {
          std::uint64_t thread_id = omp_get_thread_num();
          std::uint64_t block_beg = thread_id * max_block_size;
          std::uint64_t block_end = std::min(block_beg + max_block_size, filled);
          std::vector<std::uint64_t> local_long_lcps;
          std::uint64_t thread_lcp_sum = 0;

          for (std::uint64_t j = block_beg; j < block_end; ++j) {
            std::uint64_t i = idx_buf[2 * j];
            std::uint64_t phi_i = idx_buf[2 * j + 1];
            std::uint64_t left_idx = i;
            std::uint64_t right_idx = phi_i;
            if (!(left_halfsegment_beg <= left_idx && left_idx < left_halfsegment_end &&
                  right_halfsegment_beg <= right_idx && right_idx < right_halfsegment_end))
              std::swap(left_idx, right_idx);

            // Compute LCP value.
            std::uint64_t lcp = 0;
            while (left_idx + lcp < left_halfsegment_ext_end &&
                right_idx + lcp < right_halfsegment_ext_end &&
                left_halfsegment[left_idx - left_halfsegment_beg + lcp] ==
                right_halfsegment_ptr[right_idx - right_halfsegment_beg + lcp])
              ++lcp;

            // If the LCP computation cannot be completed, add it to the list of unfinished LCPs.
            if ((left_idx + lcp == left_halfsegment_ext_end && left_halfsegment_ext_end < text_length) ||
                (right_idx + lcp == right_halfsegment_ext_end && right_halfsegment_ext_end < text_length)) {
              ans_buf[j].m_left_idx = left_idx;
              ans_buf[j].m_right_idx = right_idx;
              ans_buf[j].m_ans = lcp;
              local_long_lcps.push_back(j);
            } else {
              std::uint64_t pos_in_B = 2UL * i + lcp;
              std::uint64_t block_id = pos_in_B / max_block_size_B;
              ans_buf[j].m_ans = pos_in_B;
              ans_buf[j].m_block_id = block_id;
              thread_lcp_sum += lcp;
            }
          }

          #pragma omp critical
          {
            // Concatenate the list of long LCP processed by a given thread with a global list.
            long_lcps.insert(long_lcps.end(), local_long_lcps.begin(), local_long_lcps.end());
            local_lcp_sum += thread_lcp_sum;
          }
        }

        // Finish the computatino of long LCPs using naive method.
        for (std::uint64_t j = 0; j < long_lcps.size(); ++j) {
          std::uint64_t which = long_lcps[j];

          // Retreive indexes from the buffer.
          std::uint64_t i = idx_buf[2 * which];
          std::uint64_t left_idx = ans_buf[which].m_left_idx;
          std::uint64_t right_idx = ans_buf[which].m_right_idx;
          std::uint64_t lcp = ans_buf[which].m_ans;

          // Compute LCP.
          lcp = naive_lcp(left_idx, right_idx, lcp, f_text, text_length, io_vol);

          // Compute answer.
          std::uint64_t pos_in_B = 2UL * i + lcp;
          std::uint64_t block_id = pos_in_B / max_block_size_B;

          // Write answer to buffer.
          ans_buf[which].m_ans = pos_in_B;
          ans_buf[which].m_block_id = block_id;

          // Update stats.
          local_lcp_sum += lcp;
        }

        // Write LCPs to file.
        for (std::uint64_t j = 0; j < filled; ++j)
          lcp_multiwriter->write_to_ith_file(ans_buf[j].m_block_id, ans_buf[j].m_ans);

#else
        for (std::uint64_t j = 0; j < filled; ++j) {
          std::uint64_t i = idx_buf[2 * j];
          std::uint64_t phi_i = idx_buf[2 * j + 1];
          std::uint64_t left_idx = i;
          std::uint64_t right_idx = phi_i;
          if (!(left_halfsegment_beg <= left_idx && left_idx < left_halfsegment_end &&
                right_halfsegment_beg <= right_idx && right_idx < right_halfsegment_end))
            std::swap(left_idx, right_idx);

          // Compute LCP value.
          std::uint64_t lcp = 0;
          while (left_idx + lcp < left_halfsegment_ext_end && right_idx + lcp < right_halfsegment_ext_end &&
              left_halfsegment[left_idx - left_halfsegment_beg + lcp] == right_halfsegment_ptr[right_idx - right_halfsegment_beg + lcp])
            ++lcp;

          // Finish the long LCP using naive method.
          if ((left_idx + lcp == left_halfsegment_ext_end && left_halfsegment_ext_end < text_length) ||
            (right_idx + lcp == right_halfsegment_ext_end && right_halfsegment_ext_end < text_length))
            lcp = naive_lcp(left_idx, right_idx, lcp, f_text, text_length, io_vol);

          // Write LCP to file.
          std::uint64_t pos_in_B = 2 * i + lcp;
          std::uint64_t block_id = pos_in_B / max_block_size_B;
          lcp_multiwriter->write_to_ith_file(block_id, pos_in_B);
          local_lcp_sum += lcp;
        }
#endif

        pairs_processed += filled;
      }

      // Update I/O volume.
      io_vol += pair_reader->bytes_read() + extra_io + n_pairs * sizeof(ext_text_offset_type);
      total_io_volume += io_vol;

      // Clean up.
      delete pair_reader;
      utils::file_delete(pairs_filename);

      // Update statistics.
      sum_irreducible_lcps += local_lcp_sum;

      // Print summary.
      long double avg_lcp = (long double)local_lcp_sum / (long double)std::max(1UL, n_pairs);
      long double elapsed = utils::wclock() - halfsegment_process_start;
      fprintf(stderr, "time = %.1Lfs, I/O = %.1LfMiB/s, avg_lcp = %.2Lf, total I/O vol = %.2Lfn\n", elapsed,
          (1.L * io_vol / (1L << 20)) / elapsed, avg_lcp, (1.L * total_io_volume) / text_length);
    }
  }

  // Clean up.
  delete[] idx_buf;
#ifdef _OPENMP
  delete[] ans_buf;
#endif
  delete lcp_multiwriter;
  std::fclose(f_text);
  free(left_halfsegment);
  free(right_halfsegment);

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "    Total time: %.2Lfs, total I/O vol = %.2Lfn\n",
      total_time, (1.L * total_io_volume) / text_length);

  return sum_irreducible_lcps;
}

struct buf_item {
  std::uint64_t m_left_idx;
  std::uint64_t m_right_idx;
  std::uint64_t m_ans;
};

template<typename text_offset_type, typename ext_text_offset_type>
std::uint64_t process_halfsegment_pairs(std::string text_filename,
    std::uint64_t text_length, std::uint64_t max_halfsegment_size,
    std::uint64_t max_overflow_size, std::string **pairs_filenames,
    std::string output_filename, std::uint64_t &total_io_volume) {
  fprintf(stderr, "  Compute irreducible LCP values:\n");
  long double start = utils::wclock();

  // Initialize basic parameters.    
  std::uint64_t n_halfsegments = (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t sum_irreducible_lcps = 0;
  
  // Open file with text.
  std::FILE *f_text = utils::file_open(text_filename, "r");

  // Initialize writer of values 2i + PLCP[i].
  typedef async_stream_writer<ext_text_offset_type> lcp_writer_type;
  lcp_writer_type *lcp_writer = new lcp_writer_type(output_filename);

  // Allocate halfsegments.
  std::uint8_t *left_halfsegment = (std::uint8_t *)malloc(max_halfsegment_size + max_overflow_size);
  std::uint8_t *right_halfsegment = (std::uint8_t *)malloc(max_halfsegment_size + max_overflow_size);

  // Allocate buffers.
  static const std::uint64_t local_buf_size = (1UL << 20);
  text_offset_type *idx_buf = new text_offset_type[local_buf_size * 2];
#ifdef _OPENMP
  buf_item *ans_buf = new buf_item[local_buf_size];
#endif

  // Processing of halfsegment pairs follows.
  for (std::uint64_t left_halfsegment_id = 0; left_halfsegment_id < n_halfsegments; ++left_halfsegment_id) {
    std::uint64_t left_halfsegment_beg = left_halfsegment_id * max_halfsegment_size;
    std::uint64_t left_halfsegment_end = std::min(left_halfsegment_beg + max_halfsegment_size, text_length);
    std::uint64_t left_halfsegment_ext_end = std::min(left_halfsegment_end + max_overflow_size, text_length);
    std::uint64_t left_halfsegment_ext_size = left_halfsegment_ext_end - left_halfsegment_beg;
    bool left_halfsegment_loaded = false;

    // Scan all halfsegments to the right of left_halfsegment_id.
    for (std::uint64_t right_halfsegment_id = left_halfsegment_id; right_halfsegment_id < n_halfsegments; right_halfsegment_id++) {
      std::uint64_t right_halfsegment_beg = right_halfsegment_id * max_halfsegment_size;
      std::uint64_t right_halfsegment_end = std::min(right_halfsegment_beg + max_halfsegment_size, text_length);
      std::uint64_t right_halfsegment_ext_end = std::min(right_halfsegment_end + max_overflow_size, text_length);
      std::uint64_t right_halfsegment_ext_size = right_halfsegment_ext_end - right_halfsegment_beg;

      // Check if that pair of halfsegments has any associated pairs.
      std::string pairs_filename = pairs_filenames[left_halfsegment_id][right_halfsegment_id];
      if (utils::file_exists(pairs_filename) == false || utils::file_size(pairs_filename) == 0) {
        if (utils::file_exists(pairs_filename))
          utils::file_delete(pairs_filename);
        continue;
      }

      // Print initial progress message.
      fprintf(stderr, "    Process halfsegments %lu and %lu: ", left_halfsegment_id, right_halfsegment_id);
      long double halfsegment_process_start = utils::wclock();
      std::uint64_t local_lcp_sum = 0;
      std::uint64_t extra_io = 0;
      std::uint64_t io_vol = 0;

      // Initialize reading from file associated with current pair of halfsegments.
      typedef async_stream_reader<text_offset_type> pair_reader_type;
      std::uint64_t n_pairs = utils::file_size(pairs_filename) / (2 * sizeof(text_offset_type));
      pair_reader_type *pair_reader = new pair_reader_type(pairs_filename);

      // Read left halfsegment from disk (if it wasn't already)
      if (left_halfsegment_loaded == false) {
        utils::read_at_offset(left_halfsegment, left_halfsegment_beg, left_halfsegment_ext_size, text_filename);
        left_halfsegment_loaded = true;
        extra_io += left_halfsegment_ext_size;
      }

      // Read right halfsegment from disk.
      std::uint8_t *right_halfsegment_ptr = right_halfsegment;
      if (right_halfsegment_id != left_halfsegment_id) {
        utils::read_at_offset(right_halfsegment, right_halfsegment_beg, right_halfsegment_ext_size, text_filename);
        extra_io += right_halfsegment_ext_size;
      } else right_halfsegment_ptr = left_halfsegment;

      std::uint64_t pairs_processed = 0;
      while (pairs_processed < n_pairs) {
        std::uint64_t filled = std::min(n_pairs - pairs_processed, local_buf_size);
        pair_reader->read(idx_buf, filled * 2);

#ifdef _OPENMP
        std::vector<std::uint64_t> long_lcps;
        std::uint64_t max_threads = omp_get_max_threads();
        std::uint64_t max_block_size = (filled + max_threads - 1) / max_threads;
        std::uint64_t n_threads = (filled + max_block_size - 1) / max_block_size;
        #pragma omp parallel num_threads(n_threads)
        {
          std::uint64_t thread_id = omp_get_thread_num();
          std::uint64_t block_beg = thread_id * max_block_size;
          std::uint64_t block_end = std::min(block_beg + max_block_size, filled);
          std::vector<std::uint64_t> local_long_lcps;
          std::uint64_t thread_lcp_sum = 0;

          for (std::uint64_t j = block_beg; j < block_end; ++j) {
            std::uint64_t i = idx_buf[2 * j];
            std::uint64_t phi_i = idx_buf[2 * j + 1];
            std::uint64_t left_idx = i;
            std::uint64_t right_idx = phi_i;
            if (!(left_halfsegment_beg <= left_idx && left_idx < left_halfsegment_end &&
                  right_halfsegment_beg <= right_idx && right_idx < right_halfsegment_end))
              std::swap(left_idx, right_idx);

            // Compute LCP value.
            std::uint64_t lcp = 0;
            while (left_idx + lcp < left_halfsegment_ext_end &&
                right_idx + lcp < right_halfsegment_ext_end &&
                left_halfsegment[left_idx - left_halfsegment_beg + lcp] ==
                right_halfsegment_ptr[right_idx - right_halfsegment_beg + lcp])
              ++lcp;

            // If the LCP computation cannot be completed, add it to the list of unfinished LCPs.
            if ((left_idx + lcp == left_halfsegment_ext_end && left_halfsegment_ext_end < text_length) ||
                (right_idx + lcp == right_halfsegment_ext_end && right_halfsegment_ext_end < text_length)) {
              ans_buf[j].m_left_idx = left_idx;
              ans_buf[j].m_right_idx = right_idx;
              ans_buf[j].m_ans = lcp;
              local_long_lcps.push_back(j);
            } else {
              std::uint64_t pos_in_B = 2UL * i + lcp;
              ans_buf[j].m_ans = pos_in_B;
              thread_lcp_sum += lcp;
            }
          }

          #pragma omp critical
          {
            // Concatenate the list of long LCP processed by a given thread with a global list.
            long_lcps.insert(long_lcps.end(), local_long_lcps.begin(), local_long_lcps.end());
            local_lcp_sum += thread_lcp_sum;
          }
        }

        // Finish the computation of long LCPs using naive method.
        for (std::uint64_t j = 0; j < long_lcps.size(); ++j) {
          std::uint64_t which = long_lcps[j];

          // Retreive indexes from the buffer.
          std::uint64_t i = idx_buf[2 * which];
          std::uint64_t left_idx = ans_buf[which].m_left_idx;
          std::uint64_t right_idx = ans_buf[which].m_right_idx;
          std::uint64_t lcp = ans_buf[which].m_ans;

          // Compute LCP.
          lcp = naive_lcp(left_idx, right_idx, lcp, f_text, text_length, io_vol);

          // Compute answer.
          std::uint64_t pos_in_B = 2UL * i + lcp;

          // Write answer to buffer.
          ans_buf[which].m_ans = pos_in_B;

          // Update stats.
          local_lcp_sum += lcp;
        }

        // Write LCPs to file.
        for (std::uint64_t j = 0; j < filled; ++j)
          lcp_writer->write(ans_buf[j].m_ans);
#else
        for (std::uint64_t j = 0; j < filled; ++j) {
          std::uint64_t i = idx_buf[2 * j];
          std::uint64_t phi_i = idx_buf[2 * j + 1];
          std::uint64_t left_idx = i;
          std::uint64_t right_idx = phi_i;
          if (!(left_halfsegment_beg <= left_idx && left_idx < left_halfsegment_end &&
                right_halfsegment_beg <= right_idx && right_idx < right_halfsegment_end))
            std::swap(left_idx, right_idx);

          // Compute LCP value.
          std::uint64_t lcp = 0;
          while (left_idx + lcp < left_halfsegment_ext_end && right_idx + lcp < right_halfsegment_ext_end &&
              left_halfsegment[left_idx - left_halfsegment_beg + lcp] == right_halfsegment_ptr[right_idx - right_halfsegment_beg + lcp])
            ++lcp;

          // Finish the computation of long LCP using naive method.
          if ((left_idx + lcp == left_halfsegment_ext_end && left_halfsegment_ext_end < text_length) ||
            (right_idx + lcp == right_halfsegment_ext_end && right_halfsegment_ext_end < text_length))
            lcp = naive_lcp(left_idx, right_idx, lcp, f_text, text_length, io_vol);

          // Write LCP to file.
          std::uint64_t pos_in_B = 2 * i + lcp;
          lcp_writer->write(pos_in_B);
          local_lcp_sum += lcp;
        }
#endif

        pairs_processed += filled;
      }

      // Update I/O volume.
      io_vol += pair_reader->bytes_read() + extra_io + n_pairs * sizeof(ext_text_offset_type);
      total_io_volume += io_vol;

      // Clean up.
      delete pair_reader;
      utils::file_delete(pairs_filename);

      // Update statistics.
      sum_irreducible_lcps += local_lcp_sum;

      // Print summary.
      long double avg_lcp = (long double)local_lcp_sum / (long double)std::max(1UL, n_pairs);
      long double elapsed = utils::wclock() - halfsegment_process_start;
      fprintf(stderr, "time = %.1Lfs, I/O = %.1LfMiB/s, avg_lcp = %.2Lf, total I/O vol = %.2Lfn\n", elapsed,
          (1.L * io_vol / (1L << 20)) / elapsed, avg_lcp, (1.L * total_io_volume) / text_length);
    }
  }

  // Clean up.
  delete[] idx_buf;
#ifdef _OPENMP
  delete[] ans_buf;
#endif
  delete lcp_writer;
  std::fclose(f_text);
  free(left_halfsegment);
  free(right_halfsegment);

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "    Total time: %.2Lfs, total I/O vol = %.2Lfn\n",
      total_time, (1.L * total_io_volume) / text_length);

  return sum_irreducible_lcps;
}

}  // namespace em_succinct_irreducible_private

#endif  // __EM_SUCCINCT_IRREDUCIBLE_SRC_PROCESS_HALFSEGMENT_PAIRS_HPP_INCLUDED
