/**
 * @file    src/lcpscan_src/sort.h
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

#ifndef __LCPSCAN_SRC_SORT_H_INCLUDED
#define __LCPSCAN_SRC_SORT_H_INCLUDED

#include <stxxl/bits/containers/sorter.h>
#include <stxxl/bits/containers/queue.h>
#include <stxxl/bits/mng/buf_istream.h>
#include <stxxl/bits/mng/buf_ostream.h>

#include <cstdio>
#include <cstdlib>
#include <list>
#include <string>
#include <numeric>
#include <algorithm>
#include <vector>

#include "./stream.h"
#include "./uint40.h"
#include "./utils.h"
#include "./tuples.h"
#include "./merge.h"
#include "./stats.h"
#include "./async_stream_reader.h"
#include "./stxxl_utils.h"


namespace lcpscan_private {

long naive_lcp(long i, long j, long lcp, std::FILE *f_text, long text_length);
void im_laca(std::string text_filename, std::string sa_filename,
    std::string out_filename, long length, bool dry_run, stats_t &stats);

//=============================================================================
// Compute a set crit_pos containing all positions i such that PLCP[i] = 0.
//=============================================================================
template<long disk_block_size>
void compute_crit_pos(
    std::string text_filename,
    std::string sa_filename,
    bool disable_cache,
    std::list<long> &crit_pos) {
  fprintf(stderr, "\rCount characters:");
  long double start = utils::wclock();

  if (disable_cache)
    utils::drop_disk_pages(text_filename);

  // Compute symbol frequencies.
  long char_count[257] = {0L}, text_length = 0L;
  typedef async_stream_reader<unsigned char> stream_reader_type;
  stream_reader_type *sr = new stream_reader_type(text_filename);
  for (; !(sr->empty()); ++text_length) {
    ++char_count[size_t(sr->read()) + 1];

    if (!(text_length & ((1L << 25) - 1))) {
      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\rCount characters: time = %.1LFs, I/O = %.2LfMiB/s",
        elapsed, (1.L * text_length / (1 << 20)) / elapsed);
    }
  }
  delete sr;

  if (disable_cache)
    utils::drop_disk_pages(text_filename);

  int sigma = 0;  // alphabet size
  for (int i = 1; i <= 256; ++i)
    if (char_count[i]) ++sigma;

  std::partial_sum(char_count, char_count + 257, char_count);

  long double elapsed = utils::wclock() - start;
  fprintf(stderr, "\rCount characters: "
      "time = %.2Lfs "
      "(+%.3Lfs/MiB), "
      "I/O = %.2LfMiB/s\n",
      elapsed,
      elapsed / (1.L * text_length / (1 << 20)),
      (1.L * text_length / (1 << 20)) / elapsed);
  fprintf(stderr, "Alphabet size = %d\n", sigma);

  if (disable_cache)
    utils::drop_disk_pages(sa_filename);

  std::FILE *fsa = utils::file_open(sa_filename, "r");
  int n_unique = std::unique(char_count, char_count + 256) - char_count;
  if (char_count[n_unique - 1] == text_length) --n_unique;

  // LCP[i] = 0 iff char_count[0..n_unique) contains i.
  for (int i = 0; i < n_unique; ++i)
    crit_pos.push_back(
        utils::read_at_offset<uint40>(char_count[i], fsa).ll() + 1);
  std::fclose(fsa);

  if (disable_cache)
    utils::drop_disk_pages(sa_filename);

  crit_pos.push_back(text_length);
  crit_pos.sort();
}

//=============================================================================
// Presort triples (SA[i], SA[i - 1], i) such that SA[i] belongs to interval
// [sblock_beg, sblock_end). The sorting key is the first element of the triple.
//
// After sorting:
//   * second components form Phi[sblock_beg, sblock_end),
//   * third components form ISA[sblock_beg, sblock_end).
//
// Output remains accessible through the triple_sorter passed as an argument.
// Returns the memory budget used to merge the runs of triples.
//=============================================================================
template<long disk_block_size>
long sort_triples(
    std::string sa_filename,
    long text_length,
    long sblock_beg,
    long sblock_end,
    stxxl::sorter<Triple, CmpTriple, disk_block_size> &triple_sorter,
    long ram_use,
    bool distribute_sa,
    bool disable_cache,
    std::vector<class stxxl::typed_block<disk_block_size, uint40>::bid_type> &sa_bids) {
  long sblock_size = sblock_end - sblock_beg;
  fprintf(stderr, "  Step I: create runs (j, ISA[j], Phi[j])\n");
  long double start = utils::wclock();

  // Scan SA and for every value in range create a triple.
  if (distribute_sa) {
    typedef stxxl::typed_block<disk_block_size, uint40> block_type;
    typedef typename block_type::bid_type bid_type;
    typedef stxxl::buf_istream<block_type, class std::vector<bid_type>::iterator> istream_type;
    istream_type SA(sa_bids.begin(), sa_bids.end(), 2);
    for (long i = 0L, prev = text_length; i < text_length; ++i) {
      uint40 SAi;
      SA >> SAi;
      long sai = SAi;
      if (sblock_beg <= sai && sai < sblock_end)
        triple_sorter.push(Triple(sai, prev, i));
      prev = sai;

      // Print progress.
      if (!(i & ((1L << 25) - 1))) {
        long double elapsed = utils::wclock() - start;
        long in_volume = sizeof(uint40) * i;
        long out_volume = sizeof(Triple) * triple_sorter.size();
        long io_volume = in_volume + out_volume;
        long double io_speed = (1.L * io_volume / (1 << 20)) / elapsed;
        fprintf(stderr, "\r    %.1Lf%%, time = %.1Lfs, "
            "I/O = %.2LfMiB/s", (100.L * i) / text_length, elapsed, io_speed);
      }
    }
  } else {
    if (disable_cache)
      utils::drop_disk_pages(sa_filename);

    async_stream_reader<uint40> SA(sa_filename);
    for (long i = 0L, prev = text_length; i < text_length; ++i) {
      long sai = SA.read();
      if (sblock_beg <= sai && sai < sblock_end)
        triple_sorter.push(Triple(sai, prev, i));
      prev = sai;

      // Print progress.
      if (!(i & ((1L << 25) - 1))) {
        long double elapsed = utils::wclock() - start;
        long in_volume = sizeof(uint40) * i;
        long out_volume = sizeof(Triple) * triple_sorter.size();
        long io_volume = in_volume + out_volume;
        long double io_speed = (1.L * io_volume / (1 << 20)) / elapsed;
        fprintf(stderr, "\r    %.1Lf%%, time = %.1Lfs, "
            "I/O = %.2LfMiB/s", (100.L * i) / text_length, elapsed, io_speed);
      }
    }

    if (disable_cache)
      utils::drop_disk_pages(sa_filename);
  }

  long double elapsed = utils::wclock() - start;
  long io_volume = sizeof(Triple) * sblock_size + 5L * text_length;
  long double io_speed = (1.L * io_volume / (1 << 20)) / elapsed;
  fprintf(stderr, "\r    Time = %.1Lfs (+%.3Lfs/MiB), I/O = %.2LfMiB/s\n",
      elapsed, elapsed / (1.L * text_length / (1 << 20)), io_speed);

  long min_triple_merge_budget = triple_sorter.min_merge_budget();
  long triple_merge_budget = std::min(ram_use / 2, min_triple_merge_budget * 2L);
  fprintf(stderr, "  Step II: merge runs (j, ISA[j], Phi[j]) / create runs (j, Phi[j])\n");
  fprintf(stderr, "    Min merge budget = %.2LfMiB. Using %.2LfMiB\n",
      1.L * min_triple_merge_budget / (1 << 20), 1.L * triple_merge_budget / (1 << 20));
  triple_sorter.sort(triple_merge_budget);
  return triple_merge_budget;
}

//=============================================================================
// Presort lexicographically pairs (Phi[i], i) for all i such that i belongs
// to superblock [sblock_beg, sblock_end) and PLCP[i] is irreducible.
//
// Pairs (Phi[i], i) and (Phi[j], j) are handled by separate sorters if i
// and j belong to different blocks.
//
// INPUT:
//   Sequence of triples (i, Phi[i], ISA[i]) passed via triple_sorter.
//
// OUTPUT:
//   As described above, returned via pair_sorters.
//   The function also stores ISA[sblock_beg..sblock_end) to disk.
//   Returns Phi[sblock_end - 1].
//=============================================================================
template<long disk_block_size>
long sort_pairs(
    long sblock_beg,
    long sblock_end,
    long max_block_size,
    long text_length,
    stats_t &stats,
    stxxl::sorter<Triple, CmpTriple, disk_block_size> *triple_sorter,
    stxxl::sorter<Pair, CmpPair, disk_block_size>** &pair_sorters,
    std::vector<class stxxl::typed_block<disk_block_size, uint40>::bid_type>* &isa_bids,
    long budget,
    std::list<long> &crit_pos,
    long prev_phi) {
  typedef stxxl::sorter<Pair, CmpPair, disk_block_size> pair_sorter_type;
  typedef stxxl::typed_block<disk_block_size, uint40> block_type;
  typedef typename block_type::bid_type bid_type;

  long sblock_size = sblock_end - sblock_beg;
  long n_block = (sblock_size + max_block_size - 1) / max_block_size;

  long double start = utils::wclock();

  // Allocate disk blocks to keep ISA.
  long disk_block_count = (sblock_size + block_type::size - 1) / block_type::size;
  isa_bids = new std::vector<bid_type>(disk_block_count);
  stxxl::block_manager *bm = stxxl::block_manager::get_instance();
  bm->new_blocks(stxxl::striping(), isa_bids->begin(), isa_bids->end());

  // Initialize ISA writer.
  typedef stxxl::buf_ostream<block_type, class std::vector<bid_type>::iterator> ostream_type;
  ostream_type isa_ostream(isa_bids->begin(), 2);

  long processed = 0;
  long irr_lcp_count = 0;

  for (long block_id = 0; block_id < n_block; ++block_id) {
    long block_beg = sblock_beg + block_id * max_block_size;
    long block_end = std::min(block_beg + max_block_size, sblock_end);

    CmpPair cmp;
    pair_sorters[block_id] = new pair_sorter_type(cmp, budget, 0L);

    for (long i = block_beg; i < block_end; ++i) {
      Triple p = *(*triple_sorter);

      // Check if PLCP[i] is irreducible.
      // Invariant:
      //   p.first = i
      //   p.second = Phi[i]
      //   p.third = ISA[i]
      if (crit_pos.front() == i || p.second.ll() != prev_phi + 1) {
        ++irr_lcp_count;
        pair_sorters[block_id]->push(Pair(p.second, p.first));
        if (crit_pos.front() == i) crit_pos.pop_front();
      }

      // Store ISA[i]. It is needed to permute PLCP into LCP.
      isa_ostream << p.third;

      // Update for next iteration.
      prev_phi = p.second;
      ++(*triple_sorter);
      ++processed;

      if (!(processed & ((1L << 25) - 1))) {
        long double elapsed = utils::wclock() - start;
        long in_volume = 15L * processed;
        long out_volume = 5L * processed + 10L * irr_lcp_count;
        long io_volume = in_volume + out_volume;
        long double io_speed = (1.L * io_volume / (1 << 20)) / elapsed;
        fprintf(stderr, "\r    %.1Lf%%, time = %.1Lfs, "
            "I/O = %.2LfMiB/s", (100.L * processed) / sblock_size,
            elapsed, io_speed);
      }
    }

    pair_sorters[block_id]->finish();
  }

  triple_sorter->finish_clear();
  delete triple_sorter;
  isa_ostream.fill(0);
  stats.irr_lcp_count += irr_lcp_count;

  long double elapsed = utils::wclock() - start;
  long io_volume = 20L * processed + 10L * irr_lcp_count;
  long double io_speed = (1.L * io_volume / (1 << 20)) / elapsed;
  fprintf(stderr, "\r    Time = %.1Lfs (+%.3Lfs/MiB), "
      "I/O = %.2LfMiB/s\n", elapsed, elapsed / (1.L * text_length / (1 << 20)),
      io_speed);

  fprintf(stderr, "    Number of irreducible LCPs = %ld\n", irr_lcp_count);

  return prev_phi;
}

//==============================================================================
// Process text block [block_beg..block_end).
//
// INPUT
//   Pairs (i, Phi[i]) where PLCP[i] is irreducible via presorted
//   phi_pair_sorter.
//
// OUTPUT
//   The collection of pairs (i, PLCP[i]) which are fed to the sorter and
//   returned.
//==============================================================================
template<long disk_block_size>
stxxl::sorter<Pair, CmpPair, disk_block_size> *process_block(
    long block_beg,
    long block_end,
    std::string text_filename,
    long text_length,
    stxxl::sorter<Pair, CmpPair, disk_block_size> *pair_sorter,
    long s_sort,
    long s_merge,
    bool distribute_text,
    bool disable_cache,
    std::vector<class stxxl::typed_block<disk_block_size, unsigned char>::bid_type> &text_bids) {
  typedef stxxl::sorter<Pair, CmpPair, disk_block_size> pair_sorter_type;

  fprintf(stderr, "\r    Process block [%ld..%ld)", block_beg, block_end);
  long double start = utils::wclock();

  // Create the output sorter.
  CmpPair cmp;
  pair_sorter_type *output = new pair_sorter_type(cmp, s_sort, 0);

  // Initialize merging.
  // NOTE: block_size was chosen so that the space required by this merger,
  // block and space needed for sorting output pairs fit in the ram limit.
  pair_sorter->sort(s_merge);

  if (pair_sorter->empty()) {
    output->finish();
    return output;
  }

  // Add one disk block to detect long LCPs.
  long overflow_end = std::min(text_length, block_end + disk_block_size);
  long ext_block_size = overflow_end - block_beg;
  unsigned char *block = new unsigned char[ext_block_size];
  if (disable_cache) utils::drop_disk_pages(text_filename);
  utils::read_at_offset(block, block_beg, ext_block_size, text_filename);
  if (disable_cache) utils::drop_disk_pages(text_filename);

  long j_prev = 0L, lcp_prev = 0L;
  long last = 0L, processed = 0L;  // to print debug

  stxxl::queue<Triple> long_lcps;  // for lcps crossing block boundary

  if (distribute_text) {
    typedef stxxl::typed_block<disk_block_size, unsigned char> block_type;
    typedef typename block_type::bid_type bid_type;
    typedef stxxl::buf_istream<block_type, class std::vector<bid_type>::iterator> istream_type;
    istream_type text(text_bids.begin(), text_bids.end(), 2);

    long text_pos = -1;
    unsigned char cur_char = 0;

    for (; !pair_sorter->empty(); ++(*pair_sorter)) {
      long i = (*pair_sorter)->second;  // s <= i < e
      long j = (*pair_sorter)->first;  // Phi[i]
      long lcp = 0L;

      lcp = std::max(lcp, lcp_prev + j_prev - j);  // PLCP property

      while (j + lcp < text_length && text_pos < j + lcp) {
        ++text_pos;
        text >> cur_char;
      }

      // Invariant: cur_char = text[text_pos]
      while (i + lcp < overflow_end && j + lcp < text_length &&
          cur_char == block[i + lcp - block_beg]) {
        ++lcp;  // extend the match
        ++text_pos;
        text >> cur_char;
      }

      if (i + lcp >= overflow_end && i + lcp < text_length && j + lcp < text_length)
        long_lcps.push(Triple(i, j, lcp));  // lcp crosses block boundary
      else output->push(Pair(i, lcp));  // PLCP[i] = lcp

      j_prev = j;
      lcp_prev = lcp;
      ++processed;

      if (j - last > (1L << 28)) {
        last = j;
        long double elapsed = utils::wclock() - start;
        long io_volume = 10L * (processed + output->size()) + j;
        long double io_speed = (1.L * io_volume / (1 << 20)) / elapsed;
        fprintf(stderr, "\r    Process block [%ld..%ld): %.1Lf%%, "
            "time = %.1Lfs, I/O = %.1LfMiB/s", block_beg, block_end,
            (100.L * j) / text_length, elapsed, io_speed);
      }
    }
  } else {
    if (disable_cache)
      utils::drop_disk_pages(text_filename);
    stream_reader<unsigned char> text(text_filename);  // text reader

    for (; !pair_sorter->empty(); ++(*pair_sorter)) {
      long i = (*pair_sorter)->second;  // s <= i < e
      long j = (*pair_sorter)->first;  // Phi[i]
      long lcp = 0L;

      lcp = std::max(lcp, lcp_prev + j_prev - j);  // PLCP property
      while (i + lcp < overflow_end && j + lcp < text_length &&
          text[j + lcp] == block[i + lcp - block_beg]) ++lcp;  // extend the match

      if (i + lcp >= overflow_end && i + lcp < text_length && j + lcp < text_length)
        long_lcps.push(Triple(i, j, lcp));  // lcp crosses block boundary
      else output->push(Pair(i, lcp));  // PLCP[i] = lcp

      j_prev = j;
      lcp_prev = lcp;
      ++processed;

      if (j - last > (1L << 28)) {
        last = j;
        long double elapsed = utils::wclock() - start;
        long io_volume = 10L * (processed + output->size()) + j;
        long double io_speed = (1.L * io_volume / (1 << 20)) / elapsed;
        fprintf(stderr, "\r    Process block [%ld..%ld): %.1Lf%%, "
            "time = %.1Lfs, I/O = %.1LfMiB/s", block_beg, block_end,
            (100.L * j) / text_length, elapsed, io_speed);
      }
    }
    if (disable_cache)
      utils::drop_disk_pages(text_filename);
  }

  delete[] block;
  pair_sorter->finish_clear();
  delete pair_sorter;

  long double elapsed = utils::wclock() - start;
  long io_volume = 10L * (processed + output->size()) + text_length;
  long double io_speed = (1.L * io_volume / (1 << 20)) / elapsed;
  fprintf(stderr, "\r    Process block [%ld..%ld): time = %.1Lfs (+%.3Lfs/MiB), "
      "I/O = %.1LfMiB/s\n", block_beg, block_end, elapsed,
      elapsed / (1.L * text_length / (1 << 20)), io_speed);

  // Compute long lcps (crossing the block boundary).
  if (!long_lcps.empty()) {
    start = utils::wclock();
    std::FILE *f_text = utils::file_open(text_filename, "r");
    for (; !long_lcps.empty(); long_lcps.pop()) {
      long i = long_lcps.front().first;
      long j = long_lcps.front().second;
      long lcp = long_lcps.front().third;

      lcp = naive_lcp(i, j, lcp, f_text, text_length);
      output->push(Pair(i, lcp));  // PLCP[i] = lcp
    }
    std::fclose(f_text);
    if (disable_cache)
      utils::drop_disk_pages(text_filename);

    elapsed = utils::wclock() - start;
    fprintf(stderr, "    Compute long lcps: %.1Lfs (+%.3Lfs/MiB)\n",
          elapsed, elapsed / (1.L * text_length / (1 << 20)));
  }

  output->finish();
  return output;
}

//==============================================================================
// Permute PLCP[sblock_beg..sblock_end) into lex order and save to disk.
// Returns PLCP[sblock_end - 1].
//==============================================================================
template<long disk_block_size>
long plcp_to_lcp(
    long sblock_beg,
    long sblock_end,
    long sblock_id,
    long max_block_size,
    stxxl::sorter<Pair, CmpPair, disk_block_size>** &lcp_sorters,
    std::vector<class stxxl::typed_block<disk_block_size, uint40>::bid_type>* &isa_bids,
    std::vector<class stxxl::typed_block<disk_block_size, uint40>::bid_type>** lcp_subsequence_bids,
    long ram_use,
    long text_length,
    stats_t &stats,
    long prev_plcp) {
  typedef stxxl::sorter<Pair, CmpPair, disk_block_size> pair_sorter_type;
  typedef stxxl::typed_block<disk_block_size, uint40> block_type;
  typedef typename block_type::bid_type bid_type;

  long sblock_size = sblock_end - sblock_beg;
  long n_block = (sblock_size + max_block_size - 1) / max_block_size;

  long min_pair_merge_budget = 0;
  for (long i = 0; i < n_block; ++i)
    min_pair_merge_budget = std::max(min_pair_merge_budget, lcp_sorters[i]->min_merge_budget());
  long pair_merge_budget = std::min(ram_use / 2, 2 * min_pair_merge_budget);

  fprintf(stderr, "  Step IV: merge runs (j, PLCP[j]) / create runs (i, LCP[i])\n");
  fprintf(stderr, "    Min merge budget = %.2LfMiB. Using %.2LfMiB\n",
      1.L * min_pair_merge_budget / (1 << 20), 1.L * pair_merge_budget / (1 << 20));

  long sorting_budget = ram_use - pair_merge_budget;

  long double start = utils::wclock();

  CmpPair cmp;
  pair_sorter_type *lex_sorter = new pair_sorter_type(cmp, sorting_budget, 0);

  typedef stxxl::buf_istream<block_type, class std::vector<bid_type>::iterator> istream_type;
  istream_type isa(isa_bids->begin(), isa_bids->end(), 2);

  long plcp_sum = 0;
  long max_lcp = 0;
  long processed = 0;
  long irr_lcp_count = 0;
  long irr_lcp_sum = 0;
  for (long block_id = 0; block_id < n_block; ++block_id) {
    long block_beg = sblock_beg + block_id * max_block_size;
    long block_end = std::min(block_beg + max_block_size, sblock_end);

    pair_sorter_type *lcp_sorter = lcp_sorters[block_id];
    lcp_sorter->sort(pair_merge_budget);

    long plcp = 0;
    for (long i = block_beg; i < block_end; ++i) {
      if (!(lcp_sorter->empty()) && (*lcp_sorter)->first.ll() == i) {
        // PLCP[i] is irreducible
        plcp = (*lcp_sorter)->second;
        ++(*lcp_sorter);
        irr_lcp_count += 1;
        irr_lcp_sum += plcp;
      } else {
        // PLCP[i] is reducible
        plcp = std::max(0L, prev_plcp - 1);
      }

      // Invariant:
      //   PLCP[i] = plcp
      //   ISA[i] = *isa
      //   LCP[*isa] = plcp
      lex_sorter->push(Pair(*isa, plcp));
      prev_plcp = plcp;
      plcp_sum += plcp;
      max_lcp = std::max(max_lcp, plcp);
      ++processed;
      ++isa;

      if (!(i & ((1L << 25) - 1))) {
        long double elapsed = utils::wclock() - start;
        long double io_volume = 15L * processed + 10L * irr_lcp_count;
        long double io_speed = (1.L * io_volume / (1 << 20)) / elapsed;
        fprintf(stderr, "\r    %.1Lf%%, time = %.1Lfs, "
            "I/O = %.2LfMiB/s", (100.L * processed) / sblock_size, elapsed, io_speed);
      }
    }

    lcp_sorter->finish_clear();
    delete lcp_sorter;
  }

  stats.irr_lcp_sum += irr_lcp_sum;
  stats.lcp_sum += plcp_sum;
  stats.max_lcp = std::max(stats.max_lcp, max_lcp);

  stxxl::block_manager *bm = stxxl::block_manager::get_instance();
  bm->delete_blocks(isa_bids->begin(), isa_bids->end());
  delete isa_bids;
  delete[] lcp_sorters;

  long double elapsed = utils::wclock() - start;
  long double io_speed = ((15.L * processed + 10.L * irr_lcp_count) / (1 << 20)) / elapsed;
  fprintf(stderr, "\r    Time = %.1Lfs (+%.3Lfs/MiB), "
      "I/O = %.2LfMiB/s\n", elapsed, elapsed / (1.L * text_length / (1 << 20)), io_speed);

  fprintf(stderr, "  Step V: merge runs (i, LCP[i])\n");
  fprintf(stderr, "    Min merge budget = %.2LfMiB. Using %.2LfMiB\n",
      1.L * lex_sorter->min_merge_budget() / (1 << 20), 1.L * ram_use / (1 << 20));

  start = utils::wclock();

  long lcp_disk_block_count = (sblock_size + block_type::size - 1) / block_type::size;
  lcp_subsequence_bids[sblock_id] = new std::vector<bid_type>(lcp_disk_block_count);

  bm->new_blocks(stxxl::striping(), lcp_subsequence_bids[sblock_id]->begin(),
      lcp_subsequence_bids[sblock_id]->end());
  typedef stxxl::buf_ostream<block_type, class std::vector<bid_type>::iterator> ostream_type;
  ostream_type lcp_ostream(lcp_subsequence_bids[sblock_id]->begin(), 2);

  // Write PLCP[sblock_beg..sblock_end) to disk permuted to lex order.
  lex_sorter->sort(ram_use);
  for (long i = 0L; i < sblock_size; ++i) {
    //   (*lex_sorter)->first = j for some j s.t. sa_beg <= SA[j] < sb_end
    //   (*lex_sorter)->second = LCP[j]
    lcp_ostream << (*lex_sorter)->second;
    ++(*lex_sorter);

    if (!(i & ((1L << 25) - 1))) {
      long double elapsed = utils::wclock() - start;
      long double io_speed = ((15.L * i) / (1 << 20)) / elapsed;
      fprintf(stderr, "\r    %.1Lf%%, time = %.1Lfs, "
          "I/O = %.2LfMiB/s", (100.L * i) / sblock_size, elapsed, io_speed);
    }
  }
  lcp_ostream.fill(0L);

  lex_sorter->finish_clear();
  delete lex_sorter;

  elapsed = utils::wclock() - start;
  io_speed = ((15.L * sblock_size) / (1 << 20)) / elapsed;
  fprintf(stderr, "\r    Time = %.1Lfs (+%.3Lfs/MiB), I/O = %.2LfMiB/s\n",
      elapsed, elapsed / (1.L * text_length / (1 << 20)), io_speed);

  return prev_plcp;
}

//==============================================================================
// LCPscan main function.
//==============================================================================
template<long disk_block_size>
void LCPscan(
    std::string text_filename,
    std::string sa_filename,
    long text_length,
    long ram_use,
    long n_sblock,
    stats_t &stats,
    std::vector<class stxxl::typed_block<disk_block_size, uint40>::bid_type> &output_bids,
    bool distribute_text,
    bool distribute_sa,
    bool disable_cache,
    std::vector<class stxxl::typed_block<disk_block_size, unsigned char>::bid_type> &text_bids,
    std::vector<class stxxl::typed_block<disk_block_size, uint40>::bid_type> &sa_bids) {
  typedef stxxl::sorter<Triple, CmpTriple, disk_block_size> triple_sorter_type;
  typedef stxxl::sorter<Pair, CmpPair, disk_block_size> pair_sorter_type;
  typedef stxxl::typed_block<disk_block_size, uint40> block_type;
  typedef typename block_type::bid_type bid_type;

  // Compute superblock sizes.
  long sblock_size[n_sblock];
  sblock_size[0] = text_length * (powl(6, n_sblock - 1) / (powl(6, n_sblock) - powl(5, n_sblock)));
  for (long i = 1; i + 1 < n_sblock; ++i)
    sblock_size[i] = (5.L * sblock_size[i - 1]) / 6.L;
  sblock_size[n_sblock - 1] = text_length;
  for (long i = 0; i + 1 < n_sblock; ++i)
    sblock_size[n_sblock - 1] -= sblock_size[i];

  fprintf(stderr, "\nNumber of superblocks = %ld\n", n_sblock);
  fprintf(stderr, "Superblock sizes:\n");
  for (long i = 0; i < n_sblock; ++i)
    fprintf(stderr, "  #%ld: %ld (%.2LfMiB, %.1Lf%%)\n", i + 1, sblock_size[i],
        1.L * sblock_size[i] / (1 << 20), (100.L * sblock_size[i]) / text_length);
  fprintf(stderr, "\n");

  long ndisks = (long)stxxl::config::get_instance()->disks_number();
  long min_merge_blocks = 4 * ndisks + 3;

  if (ram_use < min_merge_blocks * disk_block_size) {
    fprintf(stderr, "Error: not enough RAM. Need at least %ldMiB. "
        "Increase the memory limit with an option -m and try again.\n",
        (long)((1.L * min_merge_blocks * disk_block_size) / (1024L * 1024)));
    std::exit(EXIT_FAILURE);
  }

  // Compute the budget for merging pairs in 'process_block'.
  // The function below gives smooth transition when ram_use decreases.
  long s_merge = 0;
  if (ram_use >= 160L * disk_block_size)
    s_merge = (long)((80.L / (1 + (160.L * disk_block_size) / ram_use)) * disk_block_size);
  else if (ram_use >= 92L * disk_block_size)
    s_merge = 40L * disk_block_size;
  else s_merge = ram_use / 2 - 6L * disk_block_size;

  // Compute the block size for 'process_block'. The formula
  // below always gives max_block_size >= ram_use / 2.
  long max_block_size = 0;
  if (ram_use >= 160L * disk_block_size)
    max_block_size = (long)((1.L * ram_use) / (1 + (160.L * disk_block_size) / ram_use));
  else max_block_size = ram_use / 2;

  // Compute the budget for sorting pairs in 'process_block'.
  // The function below gives smooth transition when ram_use decreases.
  long s_sort = 0;
  if (ram_use >= 160L * disk_block_size)
    s_sort = (long)((80.L / (1 + (160.L * disk_block_size) / ram_use)) * disk_block_size);
  else if (ram_use >= 92L * disk_block_size)
    s_sort = ram_use / 2 - 40 * disk_block_size;
  else s_sort = 6L * disk_block_size;

  // Note: s_merge + max_block_size + s_sort = ram_use.
  fprintf(stderr, "Budget for merging runs (j, Phi[j]) = %.2LfMiB\n", 1.L * s_merge / (1 << 20));
  fprintf(stderr, "Budget for creating runs (j, PLCP[j]) = %.2LfMiB\n", 1.L * s_sort / (1 << 20));
  fprintf(stderr, "Max block size = %.2LfMiB\n\n", 1.L * max_block_size / (1 << 20));

  // Collect critical positions.
  std::list<long> crit_pos;
  compute_crit_pos<disk_block_size>(text_filename, sa_filename, disable_cache, crit_pos);

  // Process superblocks.
  std::vector<bid_type> *lcp_subsequence_bids[n_sblock];
  long sblock_beg = 0;
  long prev_phi = text_length;
  long prev_plcp = 0;
  for (long sblock_id = 0; sblock_id < n_sblock; ++sblock_id) {
    long sblock_end = sblock_beg + sblock_size[sblock_id];
    long this_sblock_size = sblock_size[sblock_id];

    fprintf(stderr, "\nProcess superblock %ld/%ld [%ld..%ld)\n",
        sblock_id + 1, n_sblock, sblock_beg, sblock_end);
    long double sblock_start = utils::wclock();

    // STEP I: sort triples (SA[i], SA[i - 1], i) by the first component.
    //         Only triples with SA[i] inside the superblock are considered.
    CmpTriple cmp_triple;
    triple_sorter_type *triple_sorter = new triple_sorter_type(cmp_triple, ram_use, 0L);
    long triple_merge_budget = sort_triples<disk_block_size>(sa_filename, text_length,
        sblock_beg, sblock_end, *triple_sorter, ram_use, distribute_sa, disable_cache, sa_bids);
    long pair_sort_budget = ram_use - triple_merge_budget;

    // STEP II: sort pairs (i, Phi[i]) by second component.
    long n_block = (this_sblock_size + max_block_size - 1) / max_block_size;
    pair_sorter_type **pair_sorters = new pair_sorter_type*[n_block];
    std::vector<bid_type> *isa_bids;
    prev_phi = sort_pairs<disk_block_size>(sblock_beg, sblock_end, max_block_size,
        text_length, stats, triple_sorter, pair_sorters, isa_bids, pair_sort_budget,
        crit_pos, prev_phi);

    // STEP III: compute irreducible PLCP values inside the superblock.
    fprintf(stderr, "  Step III: merge runs (j, Phi[j]) / create runs (j, PLCP[j])\n");
    pair_sorter_type **lcp_sorters = new pair_sorter_type*[n_block];
    for (long block_id = 0; block_id < n_block; ++block_id) {
      long block_beg = sblock_beg + block_id * max_block_size;
      long block_end = std::min(block_beg + max_block_size, sblock_end);
      lcp_sorters[block_id] = process_block<disk_block_size>(block_beg,
          block_end, text_filename, text_length, pair_sorters[block_id],
          s_sort, s_merge, distribute_text, disable_cache, text_bids);
    }
    delete[] pair_sorters;

    // STEP IV, V: add reducible PLCPs and permute (with the
    //             help of ISA) to lex-order.
    prev_plcp = plcp_to_lcp<disk_block_size>(sblock_beg, sblock_end, sblock_id,
        max_block_size, lcp_sorters, isa_bids, lcp_subsequence_bids, ram_use,
        text_length, stats, prev_plcp);

    // Print summary.
    long double sblock_total_time = utils::wclock() - sblock_start;
    fprintf(stderr, "Summary for superblock: time = %.1Lfs (+%.3Lfs/MiB), "
        "max-alloc = %.2LfB/B\n", sblock_total_time, sblock_total_time /
        (1.L * text_length / (1 << 20)), 1.L * stxxl_utils::max_alloc() / text_length);

    sblock_beg += this_sblock_size;
  }

  // STEP VI: merge LCP subsequences into final LCP array.
  if (n_sblock > 1) {
    merge<disk_block_size>(sa_filename, text_length, ram_use, n_sblock, sblock_size,
        lcp_subsequence_bids, output_bids, stats, distribute_sa, disable_cache, sa_bids);
  } else output_bids = *lcp_subsequence_bids[0];

  for (long i = 0; i < n_sblock; ++i)
    delete lcp_subsequence_bids[i];
}

//==============================================================================
// Main dispatching function. Compute the LCP array of text_filename and write
// into out_filename using at most ram_use bytes of RAM.
//==============================================================================
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
  typedef stxxl::typed_block<disk_block_size, uint40> block_type;
  typedef typename block_type::bid_type bid_type;
  long double start = utils::wclock();

  static const long kSmallInstance = (1L << 10);

  long text_length = utils::file_size(text_filename);

  // Force STXXL to print info about disks.
  (void)stxxl_utils::max_alloc();

  // Turn paths absolute.
  text_filename = utils::absolute_path(text_filename);
  sa_filename = utils::absolute_path(sa_filename);
  out_filename = utils::absolute_path(out_filename);

  // Print input parameters.
  fprintf(stderr, "\nText filename = %s\n", text_filename.c_str());
  fprintf(stderr, "SA filename = %s\n", sa_filename.c_str());
  if (!dry_run)
    fprintf(stderr, "Output filename = %s\n", out_filename.c_str());
  fprintf(stderr, "Text length = %ld (%.2LfMiB)\n", text_length, 1.L * text_length / (1 << 20));
  fprintf(stderr, "RAM use = %ld (%.2LfMiB)\n", ram_use,
      ram_use / (1024.L * 1024));
  fprintf(stderr, "Test run = %s\n", dry_run ? "TRUE" : "FALSE");
  fprintf(stderr, "Distribute text = %s\n", distribute_text ? "TRUE" : "FALSE");
  fprintf(stderr, "Distribute SA = %s\n", distribute_sa ? "TRUE" : "FALSE");
  fprintf(stderr, "Disable caching = %s\n", disable_cache ? "TRUE" : "FALSE");

  typedef stxxl::typed_block<disk_block_size, uint40> sa_block_type;
  typedef typename sa_block_type::bid_type sa_bid_type;
  long sa_block_count = (text_length + sa_block_type::size - 1) / sa_block_type::size;
  std::vector<sa_bid_type> sa_bids;

  typedef stxxl::typed_block<disk_block_size, unsigned char> text_block_type;
  typedef typename text_block_type::bid_type text_bid_type;
  long text_block_count = (text_length + text_block_type::size - 1) / text_block_type::size;
  std::vector<text_bid_type> text_bids;

  if ((distribute_text || distribute_sa) && text_length > kSmallInstance) {
    fprintf(stderr, "\nPreprocess:\n");

    if (distribute_text) {
      // Read suffix array into STXXL space.
      fprintf(stderr, "  Copy text into working space: ");
      text_bids.resize(text_block_count);
      stxxl::block_manager *bm = stxxl::block_manager::get_instance();
      bm->new_blocks(stxxl::striping(), text_bids.begin(), text_bids.end());
      typedef stxxl::buf_ostream<text_block_type, class std::vector<text_bid_type>::iterator> ostream_type;
      ostream_type text_ostream(text_bids.begin(), 2);

      if (disable_cache)
        utils::drop_disk_pages(text_filename);
      typedef stream_reader<unsigned char> stream_reader_type;
      stream_reader_type *text_istream = new stream_reader_type(text_filename);

      long double start = utils::wclock();
      for (long i = 0; i < text_length; ++i) {
        unsigned char sym = text_istream->read();
        text_ostream << sym;

        if (!(i & ((1L << 28) - 1))) {
          long double elapsed = utils::wclock() - start;
          long double io_speed = ((2.L * i) / (1 << 20)) / elapsed;
          fprintf(stderr, "\r  Copy text into working space: "
              "%.1Lf%%, I/O = %.2LfMiB/s", (100.L * i) / text_length, io_speed);
        }
      }

      text_ostream.fill(0);
      delete text_istream;

      if (disable_cache)
        utils::drop_disk_pages(text_filename);

      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\r  Copy text into working space: 100.0%%, I/O = %.2LfMiB/s\n",
          ((2.L * text_length) / (1 << 20)) / elapsed);

    }

    if (distribute_sa) {
      // Read suffix array into STXXL space.
      fprintf(stderr, "  Copy SA into working space: ");
      sa_bids.resize(sa_block_count);
      stxxl::block_manager *bm = stxxl::block_manager::get_instance();
      bm->new_blocks(stxxl::striping(), sa_bids.begin(), sa_bids.end());
      typedef stxxl::buf_ostream<sa_block_type, class std::vector<sa_bid_type>::iterator> ostream_type;
      ostream_type sa_ostream(sa_bids.begin(), 2);

      if (disable_cache)
        utils::drop_disk_pages(sa_filename);
      typedef stream_reader<uint40> stream_reader_type;
      stream_reader_type *sa_istream = new stream_reader_type(sa_filename);

      long double start = utils::wclock();
      for (long i = 0; i < text_length; ++i) {
        uint40 sai = sa_istream->read();
        sa_ostream << sai;

        if (!(i & ((1L << 25) - 1))) {
          long double elapsed = utils::wclock() - start;
          long double io_speed = ((2.L * sizeof(uint40) * i) / (1 << 20)) / elapsed;
          fprintf(stderr, "\r  Copy SA into working space: "
              "%.1Lf%%, I/O = %.2LfMiB/s", (100.L * i) / text_length, io_speed);
        }
      }

      sa_ostream.fill(0);
      delete sa_istream;

      if (disable_cache)
        utils::drop_disk_pages(sa_filename);

      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\r  Copy SA into working space: 100.0%%, I/O = %.2LfMiB/s\n",
          ((2.L * sizeof(uint40) * text_length) / (1 << 20)) / elapsed);
    }
  }

  std::vector<bid_type> lcp_bids;
  stats_t stats;
  if (text_length > kSmallInstance) {
    fprintf(stderr, "\nRunning LCPscan...\n");
    LCPscan<disk_block_size>(text_filename, sa_filename, text_length, ram_use,
        n_sblock, stats, lcp_bids, distribute_text, distribute_sa, disable_cache,
        text_bids, sa_bids);
  } else {  // for small instances run in-memory algorithm
    fprintf(stderr, "\nRunning internal memory algorithm...\n");
    im_laca(text_filename, sa_filename, out_filename, text_length, dry_run, stats);
  }

  // Delete text and SA from STXXL space.
  if ((distribute_text || distribute_sa) && text_length > kSmallInstance) {
    stxxl::block_manager *bm = stxxl::block_manager::get_instance();
    if (distribute_text) bm->delete_blocks(text_bids.begin(), text_bids.end());
    if (distribute_sa) bm->delete_blocks(sa_bids.begin(), sa_bids.end());
  }

  // Stop the timer.
  long double total_computation_time = utils::wclock() - start;

  // Write the LCP array as a proper file.
  if (!dry_run && text_length > kSmallInstance) {
    fprintf(stderr, "\nPostprocess:\n");
    fprintf(stderr, "  Copy LCP array to filesystem: ");
    start = utils::wclock();

    typedef stxxl::buf_istream<block_type, class std::vector<bid_type>::iterator> istream_type;
    long prefetch_blocks = (ram_use / 2) / disk_block_size;
    istream_type istream(lcp_bids.begin(), lcp_bids.end(), std::max(prefetch_blocks, 4L));
    stream_writer<uint40> writer(out_filename, ram_use / 2);

    for (long i = 0; i < text_length; ++i) {
      uint40 lcp_value;
      istream >> lcp_value;
      writer.write(lcp_value);

      if (!(i & ((1L << 25) - 1))) {
        long double elapsed = utils::wclock() - start;
        long double io_speed = ((2.L * sizeof(uint40) * i) / (1 << 20)) / elapsed;
        fprintf(stderr, "\r  Copy LCP array to filesystem: "
            "%.1Lf%%, I/O = %.2LfMiB/s", (100.L * i) / text_length, io_speed);
      }
    }
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "\r  Copy LCP array to filesystem: 100.0%%, I/O = %.2LfMiB/s\n",
        ((2.L * sizeof(uint40) * text_length) / (1 << 20)) / elapsed);
  }

  // Delete the output bids.
  if (text_length > kSmallInstance) {
    stxxl::block_manager *bm = stxxl::block_manager::get_instance();
    bm->delete_blocks(lcp_bids.begin(), lcp_bids.end());
  }

  // Print final summary and stats.
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  elapsed time = %.2Lfs (%.3Lfs/MiB)\n",
      total_computation_time, total_computation_time /
      (1.L * text_length / (1 << 20)));
  fprintf(stderr, "  maximal LCP = %ld\n", stats.max_lcp);
  fprintf(stderr, "  sum of LCPs = %ld\n", stats.lcp_sum);
  fprintf(stderr, "  number of irreducible LCPs = %ld\n", stats.irr_lcp_count);
  fprintf(stderr, "  sum of irreducible LCPs = %ld\n", stats.irr_lcp_sum);
}

}  // namespace lcpscan_private

#endif  // __LCPSCAN_SRC_SORT_H_INCLUDED
