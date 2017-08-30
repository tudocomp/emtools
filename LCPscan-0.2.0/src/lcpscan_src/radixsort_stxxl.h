/**
 * @file    src/lcpscan_src/radixsort_stxxl.h
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

#ifndef __LCPSCAN_SRC_RADIXSORT_STXXL_H_INCLUDED
#define __LCPSCAN_SRC_RADIXSORT_STXXL_H_INCLUDED

#include <stxxl/bits/mng/typed_block.h>
#include <stxxl/bits/mng/adaptor.h>
#include <stxxl/bits/common/types.h>

#include <cstdlib>
#include <climits>
#include <vector>
#include <algorithm>


namespace lcpscan_private {

// NOTE: A page is valid if it's aligned with the page grid and lies
// entirely in the input sequence. We do not allow partial pages.
// None of the iteratora that are input for functions below are guaranteed
// to be aligned with the page boundaries.

// Tell if `it' is the beginning of a valid page.
template<typename Iterator, unsigned pagesize_log, unsigned modulo>
inline bool is_valid_page_beg(const Iterator &end, const Iterator &it) {
  static const unsigned pagesize = (1U << pagesize_log);
  static const unsigned pagesize_mask = pagesize - 1;

  return (!(it.offset & pagesize_mask)) &&
    (pagesize <= (modulo - it.offset)) && (pagesize <= (end - it));
}


// Tell whether `it' is inside a valid page.
template<typename Iterator, unsigned pagesize_log, unsigned modulo>
inline bool is_valid_page(const Iterator &begin, const Iterator &end, const Iterator &it) {
  static const unsigned pagesize = (1U << pagesize_log);

  long pbeg = ((it.offset >> pagesize_log) << pagesize_log);
  long pend = pbeg + pagesize;
  return (pend <= modulo) && 
    ((begin.base_element != it.base_element) || ((long)begin.offset <= pbeg)) &&
    ((end.base_element != it.base_element) || (pend <= (long)end.offset));
}

// Tell the page id of the page `it' belongs to (we assume the page is valid).
template<typename Iterator, unsigned pagesize_log, unsigned modulo>
inline long get_page_id(const Iterator &begin, const Iterator &it) {
  static const unsigned pagesize = (1U << pagesize_log);

  long valid_pages_per_block = (modulo >> pagesize_log);
  long tmp = (it - begin) + begin.offset - it.offset;
  long full_blocks = tmp / modulo;
  long ret = full_blocks * valid_pages_per_block;
  ret += (it.offset >> pagesize_log); 
  ret -= std::min(valid_pages_per_block, (long)((begin.offset + pagesize - 1) >> pagesize_log));
  return std::max(0L, ret);
}

template<typename value_type>
struct page_stxxl {
  typedef page_stxxl<value_type> page_type;
  value_type *m_ptr;
  uint32_t m_id;

  page_stxxl() {}
  page_stxxl(const page_type &p) { *this = p; }

  inline page_type& operator = (const page_type &p) {
    m_ptr = p.m_ptr;
    m_id = p.m_id;
    return *this;
  }
} __attribute__((packed));


template<typename array_type, typename value_type, stxxl::unsigned_type modulo, typename Compare>
static inline void radixsort8msb_copy(stxxl::ArrayOfSequencesIterator<array_type, value_type, modulo> begin,
    stxxl::ArrayOfSequencesIterator<array_type, value_type, modulo> end,
    Compare &cmp, size_t depth) {
  typedef stxxl::ArrayOfSequencesIterator<array_type, value_type, modulo> Iterator;
  long length = end - begin;

  if (length < 128) {
    std::sort(begin, end, cmp);
    return;
  }

  const long n_buckets = 256;
  uint64_t bucketsize[n_buckets];
  std::fill(bucketsize, bucketsize + n_buckets, (uint64_t)0);

  // Fill oracle and count character occurences
  for (Iterator i = begin; i != end; ++i) {
    uint8_t v = ((uint8_t*)&(*i))[depth];
    ++bucketsize[v];
  }

  // Prefix sum
  uint64_t bucketindex[n_buckets];
  bucketindex[0] = 0;
  for (long i = 1; i < n_buckets; ++i)
    bucketindex[i] = bucketindex[i - 1] + bucketsize[i - 1];

  // Out-of-place permutation
  value_type* sorted = new value_type[length];
  for (Iterator i = begin; i != end; ++i) {
    uint8_t v = ((uint8_t*)&(*i))[depth];
    sorted[bucketindex[v]++] = *i;
  }

#if 1
  std::copy(sorted, sorted + length, begin);
#else
  Iterator dest = begin;
  value_type *src = sorted;

  long left = length;

  // Copy first disk block.
  long firstblock = std::min(left, (long)(modulo - dest.offset));
  std::copy(src, src + firstblock, dest.base_element + dest.offset);

  left -= firstblock;
  src += firstblock;
  if (left > 0) {
    // Handle all full blocks in the middle.
    long full_blocks = left / modulo;
    for (long j = 0; j < full_blocks; ++j) {
      dest.base++;
      dest.base_element = dest.base->elem;
      std::copy(src, src + modulo, dest.base_element);
      src += modulo;
    }
    left -= full_blocks * modulo;

    // Handle last disk block.
    if (left > 0) {
      dest.base++;
      dest.base_element = dest.base->elem;
      std::copy(src, src + left, dest.base_element);
    }
  }
#endif
  delete[] sorted;
  
  // Recursion into bucket
  if (depth > 0) {
    Iterator i = begin;
    for (long j = 0; j < n_buckets; i += bucketsize[j++])
      if (bucketsize[j] > 1)
        radixsort8msb_copy(i, i + bucketsize[j], cmp, depth - 1);
  }
}

template<typename array_type, typename value_type, stxxl::unsigned_type modulo,
typename Compare, unsigned pagesize_log = 10U>
void radixsort8msb_page(stxxl::ArrayOfSequencesIterator<array_type, value_type, modulo> begin,
    stxxl::ArrayOfSequencesIterator<array_type, value_type, modulo> end,
    Compare &cmp, size_t depth) {
  static const unsigned pagesize = (1U << pagesize_log);
  static const unsigned pagesize_mask = (1U << pagesize_log) - 1;
  typedef stxxl::ArrayOfSequencesIterator<array_type, value_type, modulo> Iterator;

  static_assert((pagesize <= modulo),
      "Pagesize is larger than disk block size (in elements)");

  long length = end - begin;
  if (length < 0x10000) {
    radixsort8msb_copy(begin, end, cmp, depth);
    return;
  }

  const long n_buckets = 256;
  typedef page_stxxl<value_type> page_type;
  static std::list<page_type> freepages;
  static std::list<page_type> buckets[n_buckets];

  // Actual number of pages allocated in the input can be smaller.
  long n_pages_upper_bound = ((length + pagesize - 1) >> pagesize_log);
  std::vector<page_type*> backlinks(n_pages_upper_bound, (page_type *)NULL);
  long bucketsize[n_buckets];

  for (long i = 0; i < n_buckets; ++i) {
    bucketsize[i] = 0;
    buckets[i].clear();
  }

  // Distribute the items into buckets.
  long processed = 0;
  long n_input_pages = 0;
  for (Iterator it = begin; it != end; ++it, ++processed) {
    uint8_t bucket_id = ((uint8_t*)&(*it))[depth];

    // Get new page for the bucket, if necessary.
    if (!(bucketsize[bucket_id] & pagesize_mask)) {
      page_type p;
      p.m_ptr = NULL;
      p.m_id = UINT_MAX;
      if (freepages.empty()) p.m_ptr = new value_type[pagesize];
      else {
        p = freepages.front();
        freepages.pop_front();
      }

      buckets[bucket_id].push_back(p);
      if (p.m_id != UINT_MAX)
        backlinks[p.m_id] = &(buckets[bucket_id].back());
    }

    // Place the element in the bucket.
    buckets[bucket_id].back().m_ptr[bucketsize[bucket_id] & pagesize_mask] = *it;
    ++bucketsize[bucket_id];

    // We processed the last element of the page. Add it to empty pages.
    if (processed >= pagesize_mask && ((it.offset & pagesize_mask) == pagesize_mask)) {
      page_type p;
      Iterator page_beg_it = it;
      page_beg_it.offset -= pagesize_mask;
      page_beg_it.pos -= pagesize_mask;
      p.m_ptr = &(*page_beg_it);
      p.m_id = n_input_pages++;
      freepages.push_back(p);
    }
  }


  // Process each bucket, and copy all items in that bucket to proper
  // place in the original array. This means that those positions that
  // are occupied by other pages must be moved to free space etc.
  Iterator dest = begin;
  for (long bucket_id = 0; bucket_id < n_buckets; ++bucket_id) {
    if (bucketsize[bucket_id] == 0) continue;

    typename std::list<page_type>::iterator it = buckets[bucket_id].begin();
    for (long pagebeg = 0; pagebeg < bucketsize[bucket_id]; pagebeg += pagesize, ++it) {
      long items = std::min(bucketsize[bucket_id] - pagebeg, (long)pagesize);

      // Special case, page is already in the right place.
      if (it->m_ptr == &(*dest)) {
        backlinks[it->m_id] = NULL;
        dest.offset += items;
        dest.pos += items;
        if (dest.offset == modulo) {
          // This is always allowed since any Iterator
          // must have a sentinel element.
          dest.base++;
          dest.base_element = dest.base->elem;
          dest.offset = 0;
        }
        continue;
      }

      if (dest.offset + items <= modulo) {
        Iterator it2 = dest;
        it2.offset += items - 1;
        it2.pos += items - 1;
        if (is_valid_page<Iterator, pagesize_log, modulo>(begin, end, it2)) {
          long page_evict = get_page_id<Iterator, pagesize_log, modulo>(begin, it2);
          if (backlinks[page_evict]) {
            // Find non-stale page.
            page_type tmp;
            tmp.m_ptr = NULL;
            while (!freepages.empty()) {
              tmp = freepages.front();
              freepages.pop_front();
              if (tmp.m_id < page_evict)
                tmp.m_ptr = NULL;
              else break;
            }
            if (tmp.m_ptr == NULL) {
              tmp.m_ptr = new value_type[pagesize];
              tmp.m_id = UINT_MAX;
            }

            // Copy page_evict to tmp.
            long page_evict_in_block = ((dest.offset + items - 1) >> pagesize_log);
            value_type *e =  dest.base_element + (page_evict_in_block << pagesize_log);
            std::copy(e, e + pagesize, tmp.m_ptr);

            // Update backlinks.
            *(backlinks[page_evict]) = tmp;
            if (tmp.m_id != UINT_MAX)
              backlinks[tmp.m_id] = backlinks[page_evict];
            backlinks[page_evict] = NULL;
          }
        }

        std::copy(it->m_ptr, it->m_ptr + items, dest.base_element + dest.offset);
        dest.offset += items;
        dest.pos += items;
        if (dest.offset == modulo) {
          dest.base++;
          dest.base_element = dest.base->elem;
          dest.offset = 0;
        }

        if (it->m_id != UINT_MAX)
          backlinks[it->m_id] = NULL;

        // Put page in the set of free pages.
        if (it->m_id != UINT_MAX) freepages.push_back(*it);
        else freepages.push_front(*it);
      } else {
        // dest.offset + items > modulo, i.e., the current page from the
        // bucket spans two blocks. First, copy modulo - dest.offset elements
        // from the page into the end of the current disk block. This for sure
        // does not require evicting any page from the current disk block.
        long tail_length = modulo - dest.offset;
        std::copy(it->m_ptr, it->m_ptr + tail_length, dest.base_element + dest.offset);
        dest.base++;
        dest.base_element = dest.base->elem;
        dest.offset = 0;
        dest.pos += tail_length;

        long head_length = items - tail_length;
        // Now, copy head_length element from the end of the current bucket
        // page into the new disk block. First, check if the first page of the
        // new disk block requires evicting.
        if (is_valid_page_beg<Iterator, pagesize_log, modulo>(end, dest)) {
          long evict_page = get_page_id<Iterator, pagesize_log, modulo>(begin, dest);
          if (backlinks[evict_page]) {
            // Find non-stale page.
            page_type tmp;
            tmp.m_ptr = NULL;
            while (!freepages.empty()) {
              tmp = freepages.front();
              freepages.pop_front();
              if (tmp.m_id < evict_page)
                tmp.m_ptr = NULL;
              else break;
            }
            if (tmp.m_ptr == NULL) {
              tmp.m_ptr = new value_type[pagesize];
              tmp.m_id = UINT_MAX;
            }

            // Copy page_evict to tmp.
            value_type *e = dest.base_element;
            std::copy(e, e + pagesize, tmp.m_ptr);

            // Update backlinks.
            *(backlinks[evict_page]) = tmp;
            if (tmp.m_id != UINT_MAX)
              backlinks[tmp.m_id] = backlinks[evict_page];
            backlinks[evict_page] = NULL;
          }
        }

        std::copy(it->m_ptr + tail_length, it->m_ptr + items, dest.base_element);
        dest.offset += head_length;
        dest.pos += head_length;
        if (dest.offset == modulo) {
          dest.base++;
          dest.base_element = dest.base->elem;
          dest.offset = 0;
        }

        if (it->m_id != UINT_MAX)
          backlinks[it->m_id] = NULL;

        // Put page in the set of free pages.
        if (it->m_id != UINT_MAX) freepages.push_back(*it);
        else freepages.push_front(*it);
      }
    }
  }

  // Release temporary pages.
  while (!freepages.empty()) {
    page_type p = freepages.front();
    freepages.pop_front();

    // Skip stale pages.
    if (p.m_id == UINT_MAX)
      delete[] p.m_ptr;
  }
  freepages.clear();

  // Recurse into each bucket.
  if (depth > 0) {
    Iterator i = begin;
    for (long j = 0; j < n_buckets; i += bucketsize[j++]) {
      if (bucketsize[j] > 1) {
        radixsort8msb_page<array_type, value_type, modulo, Compare,
          pagesize_log>(i, i + bucketsize[j], cmp, depth - 1);
      }
    }
  }
}

}  // namespace lcpscan_private

#endif  // __LCPSCAN_SRC_RADIXSORT_STXXL_H_INCLUDED
