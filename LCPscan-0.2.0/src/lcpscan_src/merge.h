/**
 * @file    src/lcpscan_src/merge.h
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

#ifndef __LCPSCAN_SRC_MERGE_H_INCLUDED
#define __LCPSCAN_SRC_MERGE_H_INCLUDED

#include <stxxl/bits/mng/typed_block.h>
#include <stxxl/bits/mng/mng.h>
#include <stxxl/bits/mng/buf_writer.h>
#include <stxxl/bits/mng/config.h>
#include <stxxl/bits/mng/block_prefetcher.h>
#include <stxxl/bits/algo/async_schedule.h>

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>

#include "./stream.h"
#include "./uint40.h"
#include "./utils.h"
#include "./stats.h"
#include "./stxxl_utils.h"


namespace lcpscan_private {

// Custom buf_ostream that allocates new blocks on the fly.
template <typename BlockType, class AllocStr = STXXL_DEFAULT_ALLOC_STRATEGY>
class buf_ostream_creat {
private:
    typedef BlockType block_type;
    typedef AllocStr alloc_strategy_type;
    typedef typename block_type::bid_type bid_type;

    alloc_strategy_type m_alloc_strategy;
    stxxl::buffered_writer<block_type> m_writer;
    stxxl::int_type m_current_elem;
    stxxl::unsigned_type m_block_count;
    std::vector<bid_type> *m_bids;
    stxxl::block_manager *m_blockmanager;
    block_type *m_current_block;

public:
    typedef typename block_type::const_reference const_reference;
    typedef typename block_type::reference reference;
    typedef buf_ostream_creat<block_type, alloc_strategy_type> self_type;

    buf_ostream_creat(std::vector<bid_type> &bids, stxxl::int_type nbuffers = 2)
      : m_writer(nbuffers, nbuffers / 2),
        m_current_elem(0),
        m_block_count(0),
        m_bids(&bids),
        m_blockmanager(stxxl::block_manager::get_instance()) {
      m_current_block = m_writer.get_free_block();
    }

    self_type& operator<< (const_reference record) {
      m_current_block->elem[m_current_elem++] = record;
      if (UNLIKELY(m_current_elem >= block_type::size)) {
        bid_type newbid;
        m_blockmanager->new_block(m_alloc_strategy, newbid, m_block_count++);
        m_bids->push_back(newbid);
        m_current_block = m_writer.write(m_current_block, newbid);
        m_current_elem = 0;
      }
      return (*this);
    }

    self_type& finish(const_reference padding) {
      while (m_current_elem != 0)
        operator<< (padding);
      return (*this);
    }
};

// Custom buf_istream that immediatelly deallocates disk blocks after reading.
template <typename BlockType, typename BIDIteratorType>
class buf_istream_and_destroy {
  public:
    typedef BlockType block_type;
    typedef BIDIteratorType bid_iterator_type;
    typedef stxxl::block_prefetcher<block_type, bid_iterator_type> prefetcher_type;
    typedef typename block_type::reference reference;
    typedef buf_istream_and_destroy<block_type, bid_iterator_type> self_type;

    prefetcher_type * prefetcher;
    bid_iterator_type begin_bid, end_bid;
    stxxl::int_type current_elem;
    stxxl::int_type cur_block;
    block_type *current_blk;
    stxxl::int_type *prefetch_seq;
    bool not_finished;
    long n_blocks;

    buf_istream_and_destroy(bid_iterator_type _begin, bid_iterator_type _end, stxxl::int_type nbuffers)
      : current_elem(0),
        not_finished(true) {
      begin_bid = _begin;
      cur_block = 0;
      n_blocks = _end - _begin;
      const stxxl::unsigned_type ndisks = stxxl::config::get_instance()->disks_number();
      const stxxl::int_type seq_length = _end - _begin;
      prefetch_seq = new stxxl::int_type[seq_length];
      nbuffers = stxxl::STXXL_MAX(2 * ndisks, stxxl::unsigned_type(nbuffers - 1));
      stxxl::compute_prefetch_schedule(_begin, _end, prefetch_seq, nbuffers, ndisks);
      prefetcher = new prefetcher_type(_begin, _end, prefetch_seq, nbuffers);
      current_blk = prefetcher->pull_block();
    }

    self_type& operator>> (reference record) {
      assert(not_finished);
      record = current_blk->elem[current_elem++];
      if (UNLIKELY(current_elem >= block_type::size)) {
        stxxl::block_manager *bm = stxxl::block_manager::get_instance();
        bm->delete_block(*(begin_bid + cur_block++));
        current_elem = 0;
        not_finished = prefetcher->block_consumed(current_blk);
      }
      return (*this);
    }

    ~buf_istream_and_destroy() {
      if (cur_block < n_blocks) {
        stxxl::block_manager *bm = stxxl::block_manager::get_instance();
        while (cur_block < n_blocks)
          bm->delete_block(*(begin_bid + cur_block++));
      }

      delete prefetcher;
      delete[] prefetch_seq;
    }
};


//=============================================================================
// Merge (with the help of SA) LCP subsequences stored on disk into the final
// LCP array and write the result to disk.
//=============================================================================
template<long disk_block_size>
void merge(
    std::string sa_fname,
    long text_length,
    long ram_use,
    long n_sblock,
    long *sblock_size,
    std::vector<class stxxl::typed_block<disk_block_size, uint40>::bid_type> **lcp_subsequence_bids,
    std::vector<class stxxl::typed_block<disk_block_size, uint40>::bid_type> &lcp_complete_bids,
    stats_t &stats,
    bool distribute_sa,
    bool disable_cache,
    std::vector<class stxxl::typed_block<disk_block_size, uint40>::bid_type> &sa_bids) {
  typedef stxxl::typed_block<disk_block_size, uint40> block_type;
  typedef typename block_type::bid_type bid_type;
    
  fprintf(stderr, "\nSTEP VI: merge LCP subsequences\n");
  long double merge_start = utils::wclock();

  long blocks_per_stream = std::max(2L, (ram_use / (2L * disk_block_size * (n_sblock + 1))));
  typedef buf_ostream_creat<block_type, stxxl::striping> ostream_type;
  ostream_type ostream(lcp_complete_bids, blocks_per_stream);

  typedef buf_istream_and_destroy<block_type, class std::vector<bid_type>::iterator> istream_type;
  istream_type *readers[n_sblock];
  for (long i = 0; i < n_sblock; ++i)
    readers[i] = new istream_type(lcp_subsequence_bids[i]->begin(), lcp_subsequence_bids[i]->end(), blocks_per_stream);

  long sblock_beg[n_sblock];
  sblock_beg[0] = 0L;
  for (long i = 0L; i + 1 < n_sblock; ++i)
    sblock_beg[i + 1] = sblock_beg[i] + sblock_size[i];

  uint40 lcp;
  long max_lcp = 0;

  if (distribute_sa) {
    typedef stxxl::buf_istream<block_type, class std::vector<bid_type>::iterator> istream_type;
    long istream_blocks = std::max(2L, ram_use / (2L * disk_block_size));
    istream_type SA(sa_bids.begin(), sa_bids.end(), istream_blocks);
    for (long i = 0; i < text_length; ++i) {
      uint40 SAi;
      SA >> SAi;
      long sblock_id = 0L, sai = SAi;
      while (sblock_id + 1 < n_sblock && sblock_beg[sblock_id + 1] <= sai) ++sblock_id;

      *readers[sblock_id] >> lcp;
      ostream << lcp;
      max_lcp = std::max(max_lcp, lcp.ll());

      if (!(i & ((1L << 26) - 1))) {
        long double elapsed = utils::wclock() - merge_start;
        long double io_speed = ((15.L * i) / ((1L << 20) * elapsed));
        fprintf(stderr, "\r  %.1Lf%%, time = %.1Lfs, I/O = %.2LfMiB/s",
            (100.L * i) / text_length, elapsed, io_speed);
      }
    }
  } else {
    if(disable_cache) utils::drop_disk_pages(sa_fname);
    stream_reader<uint40> *sa_reader = new stream_reader<uint40>(sa_fname, ram_use / 2);
    for (long i = 0; i < text_length; ++i) {
      long sai = sa_reader->read(), sblock_id = 0L;
      while (sblock_id + 1 < n_sblock && sblock_beg[sblock_id + 1] <= sai) ++sblock_id;

      *readers[sblock_id] >> lcp;
      ostream << lcp;
      max_lcp = std::max(max_lcp, lcp.ll());

      if (!(i & ((1L << 26) - 1))) {
        long double elapsed = utils::wclock() - merge_start;
        long double io_speed = ((15.L * i) / ((1L << 20) * elapsed));
        fprintf(stderr, "\r  %.1Lf%%, time = %.1Lfs, I/O = %.2LfMiB/s",
            (100.L * i) / text_length, elapsed, io_speed);
      }
    }
    delete sa_reader;
    if (disable_cache) utils::drop_disk_pages(sa_fname);
  }

  ostream.finish(0);
  for (long i = 0; i < n_sblock; ++i)
    delete readers[i];
  stats.max_lcp = max_lcp;

  long double text_length_mib = 1.L * text_length / (1 << 20);
  long double merge_time = utils::wclock() - merge_start;
  long double merge_time_relative = merge_time / text_length_mib;
  long double io_speed = 15.L / merge_time_relative;
  long double max_alloc = 1.L * stxxl_utils::max_alloc() / text_length;
  fprintf(stderr, "\r  Time = %.1Lfs (+%.3Lfs/MiB), I/O = %.2LfMiB/s, "
      "max-alloc = %.2LfB/B\n", merge_time, merge_time_relative,
      io_speed, max_alloc);
}

}  // namespace lcpscan_private

#endif  // __LCPSCAN_SRC_MERGE_H_INCLUDED
