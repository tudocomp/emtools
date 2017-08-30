/**
 * @file    em_succinct_irreducible_src/compute_lcp_from_plcp.hpp
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

#ifndef __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_LCP_FROM_PLCP_HPP_INCLUDED
#define __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_LCP_FROM_PLCP_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <string>
#include <algorithm>
#include <omp.h>

#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"
#include "io/async_multi_stream_writer.hpp"
#include "io/async_multipart_file_writer.hpp"
#include "io/async_multipart_multifile_reader.hpp"
#include "utils.hpp"


namespace em_succinct_irreducible_private {

template<typename text_offset_type>
void compute_lcp_from_plcp(std::uint64_t text_length, std::uint64_t ram_use,
    std::string sa_filename, std::string output_filename, std::string B_filename,
    std::uint64_t &global_io_volume, std::uint64_t &max_lcp, std::uint64_t &lcp_sum,
    bool keep_plcp = false) {
  fprintf(stderr, "Convert PLCP to LCP:\n");
  long double convert_plcp_to_lcp_start = utils::wclock();

  // Initialize basic parameters.
  std::uint64_t max_block_size = ram_use / sizeof(text_offset_type);
  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;
  std::uint64_t local_lcp_sum = 0;
  std::uint64_t local_max_lcp = 0;
  std::uint64_t total_io_volume = 0;

  // Print info about blocks.
  fprintf(stderr, "  Max block size = %lu (%.2LfMiB)\n", max_block_size, (1.L * max_block_size) / (1L << 20));
  fprintf(stderr, "  Number of blocks = %lu\n", n_blocks);

  // Set the filenames of files storing SA and LCP subsequences.
  std::string *sa_subsequences_filenames = new std::string[n_blocks];
  std::string *lcp_subsequences_filenames = new std::string[n_blocks];
  for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
    sa_subsequences_filenames[block_id] = output_filename + ".sa_subseq." + utils::intToStr(block_id) + "." + utils::random_string_hash();
    lcp_subsequences_filenames[block_id] = output_filename + ".lcp_sebseq." + utils::intToStr(block_id) + "." + utils::random_string_hash();
  }

  // Compute SA subsequences.
  {
    fprintf(stderr, "  Compute SA subsequences: ");
    long double compute_sa_subseq_start = utils::wclock();
    std::uint64_t io_volume = 0;

    // Initialize streaming of suffix array.
    typedef async_stream_reader<text_offset_type> sa_reader_type;
    sa_reader_type *sa_reader = new sa_reader_type(sa_filename);

    // Initialize multifile writer of SA subsequences.
    static const std::uint64_t n_free_buffers = 4;
    std::uint64_t total_buffers_ram = ram_use;
    std::uint64_t buffer_size = std::min((16UL << 20), total_buffers_ram / (n_blocks + n_free_buffers));
    typedef async_multi_stream_writer<text_offset_type> sa_multiwriter_type;
    sa_multiwriter_type *sa_multiwriter = new sa_multiwriter_type(buffer_size, n_free_buffers);
    for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id)
      sa_multiwriter->add_file(sa_subsequences_filenames[block_id]);

    // Read SA / write SA subsequences.
    for (std::uint64_t j = 0; j < text_length; ++j) {
      std::uint64_t sa_j = sa_reader->read();
      std::uint64_t block_id = sa_j / max_block_size;
      sa_multiwriter->write_to_ith_file(block_id, sa_j);
    }

    // Update I/O volume.
    io_volume += sa_reader->bytes_read() + sa_multiwriter->bytes_written();
    total_io_volume += io_volume;

    // Clean up.
    delete sa_reader;
    delete sa_multiwriter;

    // Print summary.
    long double compute_sa_subseq_time = utils::wclock() - compute_sa_subseq_start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n", compute_sa_subseq_time,
        ((1.L * io_volume) / (1L << 20)) / compute_sa_subseq_time, (1.L * total_io_volume) / text_length);
  }

  // Compute LCP subsequences.
  {
    fprintf(stderr, "  Compute LCP subsequences: ");
    long double compute_lcp_subseq_start = utils::wclock();
    std::uint64_t io_volume = 0;

    // Allocate the array holding the block of PLCP.
    text_offset_type *plcp_block = new text_offset_type[max_block_size];

    // Initialize reading of PLCP bitvector.
    typedef async_stream_reader<std::uint64_t> plcp_bitvector_reader_type;
    plcp_bitvector_reader_type *plcp_bitvector_reader = new plcp_bitvector_reader_type(B_filename);
    std::uint64_t bitbuf = plcp_bitvector_reader->read();
    std::uint64_t bitpos = 0;
    std::uint64_t cur_plcp = 1;

    // Allocate buffer.
    static const std::uint64_t buffer_size = (1UL << 20);
    text_offset_type *buf = new text_offset_type[buffer_size];
    text_offset_type *outbuf = new text_offset_type[buffer_size];

    // Process blocks left to right.
    for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
      std::uint64_t block_beg = block_id * max_block_size;
      std::uint64_t block_end = std::min(block_beg + max_block_size, text_length);
      std::uint64_t block_size = block_end - block_beg;

      // Read a block of PLCP into RAM.
      for (std::uint64_t j = 0; j < block_size; ++j) {
        // Increment cur_plcp for every 0 in the bitvector.
        while ((bitbuf & (1UL << bitpos)) == 0) {
          ++cur_plcp;
          ++bitpos;
          if (bitpos == 64) {
            bitbuf = plcp_bitvector_reader->read();
            bitpos = 0;
          }
        }

        // We decrement last because cur_plcp is unsigned.
        --cur_plcp;
        plcp_block[j] = cur_plcp;

        // Skip the 1-bit in the bitvector.
        ++bitpos;
        if (bitpos == 64) {
          if (plcp_bitvector_reader->empty() == false)
            bitbuf = plcp_bitvector_reader->read();
          bitpos = 0;
        }
      }

      // Compute LCP subsequence and write to file.
      {
        // Initialize SA subsequence reader.
        typedef async_stream_reader<text_offset_type> sa_subseq_reader_type;
        sa_subseq_reader_type *sa_subseq_reader =
          new sa_subseq_reader_type(sa_subsequences_filenames[block_id]);

        // Initialize LCP subsequence writer.
        std::uint64_t single_file_max_bytes = text_length / (n_blocks * 2UL);  // 10UL
        typedef async_multipart_file_writer<text_offset_type> lcp_subseq_writer_type;
        lcp_subseq_writer_type *lcp_subseq_writer =
          new lcp_subseq_writer_type(lcp_subsequences_filenames[block_id], single_file_max_bytes);

        // Compute LCP subsequence.
        std::uint64_t subseq_size = utils::file_size(sa_subsequences_filenames[block_id]) / sizeof(text_offset_type);
        std::uint64_t items_processed = 0;
        while (items_processed < subseq_size) {
          std::uint64_t filled = std::min(buffer_size, subseq_size - items_processed);
          sa_subseq_reader->read(buf, filled);
#ifdef _OPENMP
          #pragma omp parallel for
          for (std::uint64_t j = 0; j < filled; ++j) {
            std::uint64_t sa_val = buf[j];
            std::uint64_t lcp_val = plcp_block[sa_val - block_beg];
            outbuf[j] = lcp_val;
          }
#else
          for (std::uint64_t j = 0; j < filled; ++j) {
            std::uint64_t sa_val = buf[j];
            std::uint64_t lcp_val = plcp_block[sa_val - block_beg];
            outbuf[j] = lcp_val;
          }
#endif
          lcp_subseq_writer->write(outbuf, filled);
          items_processed += filled;
        }

        // Update I/O volume.
        io_volume += sa_subseq_reader->bytes_read() + lcp_subseq_writer->bytes_written();

        // Clean up.
        delete sa_subseq_reader;
        delete lcp_subseq_writer;
      }

      utils::file_delete(sa_subsequences_filenames[block_id]);
    }

    // Update I/O volume.
    io_volume += plcp_bitvector_reader->bytes_read();
    total_io_volume += io_volume;

    // Clean up.
    delete[] buf;
    delete[] outbuf;
    delete[] plcp_block;
    delete plcp_bitvector_reader;
    if (keep_plcp == false)
      utils::file_delete(B_filename);

    // Print summary.
    long double compute_lcp_subseq_time = utils::wclock() - compute_lcp_subseq_start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n", compute_lcp_subseq_time,
        ((1.L * io_volume) / (1L << 20)) / compute_lcp_subseq_time, (1.L * total_io_volume) / text_length);
  }

  // Merge LCP subsequences.
  {
    fprintf(stderr, "  Merge LCP subsequences: ");
    long double merge_lcp_subseq_start = utils::wclock();
    std::uint64_t io_volume = 0;

    // Initialize the reader of LCP subsequences.
    std::uint64_t total_buffers_ram = ram_use;
    std::uint64_t buffer_size = total_buffers_ram / (2UL * n_blocks);
    typedef async_multipart_multifile_reader<text_offset_type> lcp_subseq_multireader_type;
    lcp_subseq_multireader_type *lcp_subseq_multireader =
      new lcp_subseq_multireader_type(n_blocks, buffer_size);
    for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id)
      lcp_subseq_multireader->add_file(lcp_subsequences_filenames[block_id]);

    // Initialize the writer of the final LCP array.
    typedef async_stream_writer<text_offset_type> lcp_writer_type;
    lcp_writer_type *lcp_writer = new lcp_writer_type(output_filename);

    // Initialize the reader of SA.
    typedef async_stream_reader<text_offset_type> sa_reader_type;
    sa_reader_type *sa_reader = new sa_reader_type(sa_filename);

    // Compute final LCP.
    for (std::uint64_t j = 0; j < text_length; ++j) {
      std::uint64_t sa_j = sa_reader->read();
      std::uint64_t block_id = sa_j / max_block_size;
      std::uint64_t lcp_j = lcp_subseq_multireader->read_from_ith_file(block_id);
      local_max_lcp = std::max(local_max_lcp, lcp_j);
      local_lcp_sum += lcp_j;
      lcp_writer->write(lcp_j);
    }

    // Update I/O volume.
    io_volume += sa_reader->bytes_read() + lcp_subseq_multireader->bytes_read() + lcp_writer->bytes_written();
    total_io_volume += io_volume;

    // Clean up.
    delete sa_reader;
    delete lcp_writer;
    delete lcp_subseq_multireader;

    // Print summary.
    long double merge_lcp_subseq_time = utils::wclock() - merge_lcp_subseq_start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n", merge_lcp_subseq_time,
        ((1.L * io_volume) / (1L << 20)) / merge_lcp_subseq_time, (1.L * total_io_volume) / text_length);
  }

  // Clean up.
  delete[] sa_subsequences_filenames;
  delete[] lcp_subsequences_filenames;

  // Print summary.
  long double convert_plcp_to_lcp_time = utils::wclock() - convert_plcp_to_lcp_start;
  fprintf(stderr, "Summary: time = %.2Lfs, total I/O vol = %.2Lfn\n",
      convert_plcp_to_lcp_time, (1.L * total_io_volume) / text_length);

  // Update reference variables.
  global_io_volume += total_io_volume;
  max_lcp = local_max_lcp;
  lcp_sum = local_lcp_sum;
}

template<typename text_offset_type>
void compute_lcp_from_plcp(std::uint64_t text_length, std::uint64_t ram_use, std::uint64_t *B,
    std::string sa_filename, std::string output_filename, std::uint64_t &total_io_volume,
    std::uint64_t &max_lcp, std::uint64_t &lcp_sum) {
  // Write B to disk.
  std::string B_filename = output_filename + ".plcp." + utils::random_string_hash();
  {
    // Start the timer.
    fprintf(stderr, "Write PLCP bitvector to disk: ");
    long double write_plcp_start = utils::wclock();
    std::uint64_t io_volume = 0;

    // Write the data.
    std::uint64_t length_of_B_in_words = (2UL * text_length + 63) / 64;
    utils::write_to_file(B, length_of_B_in_words, B_filename);

    // Update I/O volume.
    io_volume += length_of_B_in_words * sizeof(std::uint64_t);
    total_io_volume += io_volume;
    long double write_plcp_time = utils::wclock() - write_plcp_start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, I/O vol = %.2Lfn\n\n", write_plcp_time,
        ((1.L * io_volume) / (1L << 20)) / write_plcp_time, (1.L * io_volume) / text_length);
  }
  delete[] B;

  // Convert PLCP to LCP using EM method.
  compute_lcp_from_plcp<text_offset_type>(text_length, ram_use, sa_filename,
      output_filename, B_filename, total_io_volume, max_lcp, lcp_sum);
}

template<typename text_offset_type>
void compute_lcp_from_plcp(std::string input_filename, std::string sa_filename,
    std::string output_filename, std::uint64_t ram_use) {
  srand(time(0) + getpid());
  utils::drop_disk_pages(input_filename);
  utils::drop_disk_pages(sa_filename);
  long double start = utils::wclock();
  std::uint64_t io_volume = 0;

  // Compute basic parameters.
  std::uint64_t text_length = utils::file_size(sa_filename) / sizeof(text_offset_type);
  long double text_to_ram_ratio = (long double)text_length / (long double)ram_use;

  if (text_length == 0) {
    fprintf(stderr, "Error: the input file is empty!\n");
    std::exit(EXIT_FAILURE);
  }

  // Turn paths absolute.
  input_filename = utils::absolute_path(input_filename);
  sa_filename = utils::absolute_path(sa_filename);
  output_filename = utils::absolute_path(output_filename);

  // Print summary of basic parameters.
  fprintf(stderr, "PLCP filename = %s\n", input_filename.c_str());
  fprintf(stderr, "SA filename = %s\n", sa_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n", text_length, 1.L * text_length / (1 << 20));
  fprintf(stderr, "Text size / ram_use = %.2Lf\n", text_to_ram_ratio);
  fprintf(stderr, "RAM use = %lu (%.2LfMiB)\n", ram_use, ram_use / (1024.L * 1024));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n", sizeof(text_offset_type));
#ifdef _OPENMP
  fprintf(stderr, "Max number of threads = %d\n", omp_get_max_threads());
#endif
  fprintf(stderr, "\n");

  std::uint64_t lcp_sum = 0;
  std::uint64_t max_lcp = 0;

  // Convert the PCLP array (bitvector representation) to LCP array.
  compute_lcp_from_plcp<text_offset_type>(text_length, ram_use, sa_filename,
      output_filename, input_filename, io_volume, max_lcp, lcp_sum, true);

  // Print summary.
  long double total_time = utils::wclock() - start;
  long double avg_lcp = (long double)lcp_sum / text_length;
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  elapsed time = %.2Lfs (%.3Lfs/MiB of text)\n", total_time, total_time / (1.L * text_length / (1L << 20)));
  fprintf(stderr, "  speed = %.2LfMiB of text/s\n", (1.L * text_length / (1L << 20)) / total_time);
  fprintf(stderr, "  I/O volume = %lu (%.2Lfbytes/input symbol)\n", io_volume, (1.L * io_volume) / text_length);
  fprintf(stderr, "  sum of all LCPs = %lu\n", lcp_sum);
  fprintf(stderr, "  average LCP = %.2Lf\n", avg_lcp);
  fprintf(stderr, "  maximal LCP = %lu\n", max_lcp);
}

}  // namespace em_succinct_irreducible_private

#endif  // __EM_SUCCINCT_IRREDUCIBLE_SRC_COMPUTE_LCP_FROM_PLCP_HPP_INCLUDED
