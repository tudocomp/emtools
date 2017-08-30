/**
 * @file    src/lcpscan_src/radixsort.h
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

#ifndef __LCPSCAN_SRC_RADIXSORT_H_INCLUDED
#define __LCPSCAN_SRC_RADIXSORT_H_INCLUDED

#include <stdint.h>
#include <vector>
#include <stack>
#include <list>
#include <algorithm>


namespace lcpscan_private {

template <typename Iterator, typename Compare>
static inline void radixsort8msb(Iterator begin, Iterator end, Compare &cmp, size_t depth) {
  if (end - begin < 128) {
    std::sort(begin, end, cmp);
    return;
  }

  const size_t K = 256;

  size_t bucketsize[K];
  memset(bucketsize, 0, K * sizeof(bucketsize[0]));

  // fill oracle and count character occurances
  typedef uint8_t oracle_type;
  oracle_type* oracle = (oracle_type*)malloc((end - begin) * sizeof(oracle_type));
  size_t ic = 0, jc;
  for (Iterator i = begin; i != end; ++i, ++ic) {
    uint8_t v = ((uint8_t*)&(*i))[depth];
    ++bucketsize[ oracle[ic] = v ];
  }

  // prefix sum
  ssize_t bucketindex[K];
  bucketindex[0] = bucketsize[0];
  size_t last_bucket_size = bucketsize[0];
  for (size_t i = 1; i < K; ++i) {
    bucketindex[i] = bucketindex[i-1] + bucketsize[i];
    if (bucketsize[i]) last_bucket_size = bucketsize[i];
  }

  // in-place permutation
  ic = 0;
  for (Iterator i = begin, j; i < end - last_bucket_size; ) {
    while ( (jc = --bucketindex[ oracle[ic] ]) > ic ) {
      j = begin + jc;
      std::swap(*i, *j);
      std::swap(oracle[ic], oracle[jc]);
    }
    i  += bucketsize[ oracle[ic] ];
    ic += bucketsize[ oracle[ic] ];
  }
  free(oracle);

  if (depth == 0) return;

  // recursion into bucket
  Iterator i = begin;
  for (size_t j = 0; j < K; i += bucketsize[j++]) {
    if (bucketsize[j] <= 1) continue;
    radixsort8msb(i, i + bucketsize[j], cmp, depth-1);
  }
}


template<typename value_type, typename Compare>
static inline void radixsort8msb_copy(value_type* tab, size_t length,
    Compare &cmp, size_t depth) {
  if (length < 128) {
    std::sort(tab, tab + length, cmp);
    return;
  }

  const size_t n_buckets = 256;
  uint64_t bucketsize[n_buckets];
  std::fill(bucketsize, bucketsize + n_buckets, (uint64_t)0);

  // fill oracle and count character occurences
  for (size_t i = 0; i < length; ++i) {
    uint8_t v = ((uint8_t*)&(tab[i]))[depth];
    ++bucketsize[v];
  }

  // Prefix sum.
  uint64_t bucketindex[n_buckets];
  bucketindex[0] = 0;
  for (size_t i = 1; i < n_buckets; ++i)
    bucketindex[i] = bucketindex[i - 1] + bucketsize[i - 1];

  // Out-of-place permutation
  value_type *sorted = new value_type[length];
  for (size_t i = 0; i < length; ++i) {
    uint8_t v = ((uint8_t*)&(tab[i]))[depth];
    sorted[bucketindex[v]++] = tab[i];
  }

  std::copy(sorted, sorted + length, tab);
  delete[] sorted;
  
  // recursion into bucket
  if (depth > 0) {
    value_type *ptr = tab;
    for (size_t j = 0; j < n_buckets; ptr += bucketsize[j++])
      if (bucketsize[j] > 1)
        radixsort8msb_copy(ptr, bucketsize[j], cmp, depth - 1);
  }
}


template<typename value_type, typename Compare, unsigned pagesize_log = 10U>
static inline void radixsort8msb_page(value_type* tab, size_t length,
    Compare &cmp, size_t depth) {
  static const unsigned pagesize = (1U << pagesize_log);
  static const unsigned pagesize_mask = (1U << pagesize_log) - 1;
  size_t n_pages = (length + pagesize - 1) / pagesize;

  if (length < 0x10000) {
    radixsort8msb_copy(tab, length, cmp, depth);
    return;
  }

  const size_t n_buckets = 256;
  static std::list<value_type *> freepages;
  static std::list<value_type *> buckets[n_buckets];

  std::vector<value_type **> backlinks(n_pages, (value_type**)NULL);
  size_t bucketsize[n_buckets];

  for (uint64_t i = 0; i < n_buckets; ++i) {
    bucketsize[i] = 0;
    buckets[i].clear();
  }

  // Distribute the items into buckets.
  value_type *ptr = tab;
  for (uint64_t i = 0; i < length; ++i, ++ptr) {
    uint8_t bucket_id = ((uint8_t*)&(*ptr))[depth];

    // Get new page for the bucket, if necessary.
    if (!(bucketsize[bucket_id] & pagesize_mask)) {
      value_type *p = NULL;
      if (freepages.empty()) p = new value_type[pagesize];
      else {
        p = freepages.front();
        freepages.pop_front();
      }

      buckets[bucket_id].push_back(p);
      if (tab <= p && p < tab + length)
        backlinks[(p - tab) >> pagesize_log] = &(buckets[bucket_id].back());
    }

    // Place the element in the bucket.
    buckets[bucket_id].back()[bucketsize[bucket_id] & pagesize_mask] = *ptr;
    ++bucketsize[bucket_id];

    // We processed the last element of the page. Add it to empty pages.
    if ((i & pagesize_mask) == pagesize_mask)
      freepages.push_back(ptr - pagesize_mask);
  }

  // Process each bucket, and copy all items in that bucket to proper
  // place in the original array. This means that those positions that
  // are occupied by other pages must be moved to free space etc.
  uint64_t filled = 0;
  for (uint64_t bucket_id = 0; bucket_id < n_buckets; ++bucket_id) {
    if (bucketsize[bucket_id] == 0) continue;

    typename std::list<value_type*>::iterator it = buckets[bucket_id].begin();
    for (uint64_t pagebeg = 0; pagebeg < bucketsize[bucket_id]; pagebeg += pagesize, ++it) {
      uint64_t items = std::min(bucketsize[bucket_id] - pagebeg, (size_t)pagesize);
      uint64_t page_evict = ((filled + items - 1) >> pagesize_log);

      // Special case, page is already in the right place.
      if (*it == tab + filled) {
        backlinks[filled >> pagesize_log] = NULL;
        filled += items;
        continue;
      }

      if (backlinks[page_evict]) {
        // Find non-stale free page.
        value_type *tmp = NULL;
        while (!freepages.empty()) {
          tmp = freepages.front();
          freepages.pop_front();
          if (tab <= tmp && tmp < tab + filled)
            tmp = NULL;
          else break;
        }

        // If cannot find non-stale page, allocate new.
        if (tmp == NULL)
          tmp = new value_type[pagesize];

        // Copy page_evict to tmp.
        value_type *e = tab + (page_evict << pagesize_log);
        std::copy(e, e + pagesize, tmp);

        // Update backlinks.
        *(backlinks[page_evict]) = tmp;
        if (tab <= tmp && tmp < tab + length)
          backlinks[(tmp - tab) >> pagesize_log] = backlinks[page_evict];
        backlinks[page_evict] = NULL;
      }

      if (tab <= *it && *it < tab + length)
        backlinks[(*it - tab) >> pagesize_log] = NULL;

      std::copy(*it, (*it) + items, tab + filled);
      filled += items;

      // Put page in the set of free pages.
      if (tab <= (*it) && (*it) < tab + length) freepages.push_back(*it);
      else freepages.push_front(*it);
    }
  }

  // Release free pages.
  while (!freepages.empty()) {
    value_type *p = freepages.front();
    freepages.pop_front();

    // Skip stale pages.
    if (p < tab || p >= tab + length)
      delete[] p;
  }
  freepages.clear();

  // Recurse into each bucket.
  if (depth > 0) {
    ptr = tab;
    for (uint64_t j = 0; j < n_buckets; ptr += bucketsize[j++])
      if (bucketsize[j] > 1)
        radixsort8msb_page(ptr, bucketsize[j], cmp, depth - 1);
  }
}


template<typename Iterator>
struct page {
  typedef typename std::iterator_traits<Iterator>::value_type value_type;
  typedef page<Iterator> page_type;

  Iterator m_it;
  value_type *m_ptr;

  page() {}
  page(const page_type &p) { *this = p; }

  inline page_type& operator = (const page_type &p) {
    m_it = p.m_it;
    m_ptr = p.m_ptr;
    return *this;
  }

  inline void set_and_advance(const value_type &v) {
    if (!m_ptr) *m_it++ = v;
    else *m_ptr++ = v;
  }

  inline value_type& operator[] (size_t i) const {
    if (m_ptr == NULL) return *(m_it + i);
    else return m_ptr[i];
  }
};


template<typename Iterator, typename Compare>
static inline void radixsort8msb_copy(Iterator begin, Iterator end,
    Compare &cmp, size_t depth) {
  typedef typename std::iterator_traits<Iterator>::value_type value_type;

  long length = end - begin;
  if (length < 128) {
    std::sort(begin, end, cmp);
    return;
  }

  const size_t n_buckets = 256;
  uint64_t bucketsize[n_buckets];
  std::fill(bucketsize, bucketsize + n_buckets, (uint64_t)0);

  // Compute bucket sizes.
  for (Iterator it = begin; it != end; ++it) {
    uint8_t v = ((uint8_t*)&(*it))[depth];
    ++bucketsize[v];
  }

  // Prefix sum.
  uint64_t bucketindex[n_buckets];
  bucketindex[0] = 0;
  for (size_t i = 1; i < n_buckets; ++i)
    bucketindex[i] = bucketindex[i - 1] + bucketsize[i - 1];

  // Out-of-place permutation.
  value_type *sorted = new value_type[length];
  for (Iterator it = begin; it != end; ++it) {
    uint8_t v = ((uint8_t*)&(*it))[depth];
    sorted[bucketindex[v]++] = *it;
  }

  std::copy(sorted, sorted + length, begin);
  delete[] sorted;
  
  // recursion into bucket
  if (depth > 0) {
    Iterator i = begin;
    for (size_t j = 0; j < n_buckets; i += bucketsize[j++])
      if (bucketsize[j] > 1)
        radixsort8msb_copy(i, i + bucketsize[j], cmp, depth - 1);
  }
}

template<typename Iterator, typename Compare, unsigned pagesize_log = 10U>
static inline void radixsort8msb_page(Iterator begin, Iterator end, Compare &cmp, size_t depth) {
  static const unsigned pagesize = (1U << pagesize_log);
  static const unsigned pagesize_mask = (1U << pagesize_log) - 1;
  typedef typename std::iterator_traits<Iterator>::value_type value_type;

  // The upper bound on the number of pages from in input.
  size_t length = end - begin;
  size_t n_pages_upper_bound = (length + pagesize - 1) / pagesize;

  if (length < 0x10000) {
    radixsort8msb_copy(begin, end, cmp, depth);
    return;
  }

  const size_t n_buckets = 256;
  typedef page<Iterator> page_type;
  static std::list<page_type> freepages;
  static std::list<page_type> buckets[n_buckets];
  static page_type bucket_ptr[n_buckets];
  std::vector<page_type*> backlinks(n_pages_upper_bound, (page_type *)NULL);
  size_t bucketsize[n_buckets];

  for (uint64_t i = 0; i < n_buckets; ++i) {
    bucketsize[i] = 0;
    buckets[i].clear();
  }

  // Distribute the items into buckets.
  long processed = 0;
  Iterator prev_page_it = begin;
  for (Iterator it = begin; it != end; ++it, ++processed) {
    uint8_t bucket_id = ((uint8_t*)&(*it))[depth];

    // Get new page for the bucket, if necessary.
    if (!(bucketsize[bucket_id] & pagesize_mask)) {
      page_type p;
      p.m_it = Iterator();
      p.m_ptr = NULL;
      if (freepages.empty()) p.m_ptr = new value_type[pagesize];
      else {
        p = freepages.front();
        freepages.pop_front();
      }

      buckets[bucket_id].push_back(p);
      bucket_ptr[bucket_id] = buckets[bucket_id].back();
      if (p.m_ptr == NULL)
        backlinks[(p.m_it - begin) >> pagesize_log] = &(buckets[bucket_id].back());
    }

    // Place the element in the bucket.
    bucket_ptr[bucket_id].set_and_advance(*it);
    ++bucketsize[bucket_id];

    // We processed the last element of the page. Add it to empty pages.
    if ((processed & pagesize_mask) == pagesize_mask) {
      page_type p;
      p.m_it = prev_page_it;
      p.m_ptr = NULL;
      freepages.push_back(p);
      prev_page_it = it;
      ++prev_page_it;
    }
  }

  // Process each bucket, and copy all items in that bucket to proper
  // place in the original array. This means that those positions that
  // are occupied by other pages must be moved to free space etc.
  Iterator dest = begin;
  for (uint64_t bucket_id = 0; bucket_id < n_buckets; ++bucket_id) {
    if (bucketsize[bucket_id] == 0) continue;

    typename std::list<page_type>::iterator it = buckets[bucket_id].begin();
    for (uint64_t pagebeg = 0; pagebeg < bucketsize[bucket_id]; pagebeg += pagesize, ++it) {
      uint64_t items = std::min(bucketsize[bucket_id] - pagebeg, (size_t)pagesize);
      long page_evict = (((dest - begin) + items - 1) >> pagesize_log);

      // Special case, page is already in the right place.
      if (it->m_ptr == NULL && it->m_it == dest) {
        backlinks[(it->m_it - begin) >> pagesize_log] = NULL;
        dest += items;
        continue;
      }

      if (backlinks[page_evict]) {
        // Find non-stale free page.
        page_type tmp;
        tmp.m_it = Iterator();
        tmp.m_ptr = NULL;

        bool found = false;
        while (!freepages.empty()) {
          tmp = freepages.front();
          freepages.pop_front();

          // Check if the page is stale.
          if (!tmp.m_ptr && (long)((tmp.m_it - begin) >> pagesize_log) < page_evict) {}
          else {
            found = true;
            break;
          }
        }

        // If cannot find non-stale page, allocate new.
        if (!found)
          tmp.m_ptr = new value_type[pagesize];

        // Copy page_evict to tmp.
        Iterator e = begin + (page_evict << pagesize_log);
        if (tmp.m_ptr != NULL) std::copy(e, e + pagesize, tmp.m_ptr);
        else std::copy(e, e + pagesize, tmp.m_it);

        // Update backlinks.
        *(backlinks[page_evict]) = tmp;
        if (tmp.m_ptr == NULL)
          backlinks[(tmp.m_it - begin) >> pagesize_log] = backlinks[page_evict];
        backlinks[page_evict] = NULL;
      }

      if (it->m_ptr == NULL) {
        backlinks[(it->m_it - begin) >> pagesize_log] = NULL;
        std::copy(it->m_it, it->m_it + items, dest);
      } else std::copy(it->m_ptr, it->m_ptr + items, dest);
      dest += items;

      // Put page in the set of free pages.
      if (it->m_ptr == NULL) freepages.push_back(*it);
      else freepages.push_front(*it);
    }
  }

  // Release temporary pages.
  while (!freepages.empty()) {
    page_type p = freepages.front();
    freepages.pop_front();

    // Skip stale pages.
    if (p.m_ptr != NULL)
      delete[] p.m_ptr;
  }
  freepages.clear();

  // Recurse into each bucket.
  if (depth > 0) {
    Iterator i = begin;
    for (uint64_t j = 0; j < n_buckets; i += bucketsize[j++])
      if (bucketsize[j] > 1)
        radixsort8msb_page(i, i + bucketsize[j], cmp, depth - 1);
  }
}

}  // namespace lcpscan_private

#endif  // __LCPSCAN_SRC_RADIXSORT_H_INCLUDED
