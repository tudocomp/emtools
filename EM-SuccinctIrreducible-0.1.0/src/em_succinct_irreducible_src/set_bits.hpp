/**
 * @file    em_succinct_irreducible_src/set_bits.hpp
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

#ifndef __EM_SUCCINCT_IRREDUCIBLE_SRC_SET_BITS_HPP_INCLUDED
#define __EM_SUCCINCT_IRREDUCIBLE_SRC_SET_BITS_HPP_INCLUDED

#include <cstdint>
#include <vector>
#include <algorithm>
#include <omp.h>


namespace em_succinct_irreducible_private {

void set_bits(std::uint64_t *bv, std::uint64_t *tab, std::uint64_t tab_size) {
  for (std::uint64_t i = 0; i < tab_size; ++i) {
    std::uint64_t idx = tab[i];
    bv[idx >> 6] |= (1UL << (idx & 63));
  }
}

#ifdef _OPENMP
template<typename int_type>
void permute_into_small_buckets(int_type *tab,
    int_type *temp, std::uint64_t length,
    std::uint64_t lower_bound, std::uint64_t upper_bound,
    std::uint64_t max_bucket_size,
    std::vector<std::uint64_t> &output_bucket_sizes) {
  // Move all items into temp array.
  #pragma omp parallel for
  for (std::uint64_t j = 0; j < length; ++j)
    temp[j] = tab[j];

  // Compute bucket range. Note that bucket range is understood
  // as the length of ranges of keys assigned to a bucket;
  // bucket size is the number of items inside the bucket.
  static const std::uint64_t max_buckets = 1024;
  std::uint64_t value_range = upper_bound - lower_bound;
  std::uint64_t bucket_range_log = 6;
  std::uint64_t bucket_range = 64;
  while ((value_range + bucket_range - 1) / bucket_range > max_buckets) {
    ++bucket_range_log;
    bucket_range <<= 1;
  }
  std::uint64_t n_buckets = (value_range + bucket_range - 1) / bucket_range;

  // Allocate bucket counts.
  std::uint64_t max_threads = omp_get_max_threads();
  std::uint64_t max_range_size = (length + max_threads - 1) / max_threads;
  std::uint64_t n_threads = (length + max_range_size - 1) / max_range_size;
  std::uint64_t **bucket_ptr = new std::uint64_t*[n_threads];
  for (std::uint64_t thread_id = 0; thread_id < n_threads; ++thread_id)
    bucket_ptr[thread_id] = new std::uint64_t[n_buckets];

  // Permute items into buckets.
  #pragma omp parallel num_threads(n_threads)
  {
    std::uint64_t thread_id = omp_get_thread_num();
    std::uint64_t range_beg = thread_id * max_range_size;
    std::uint64_t range_end = std::min(range_beg + max_range_size, length);
    std::uint64_t *local_bucket_ptr = bucket_ptr[thread_id];
    std::fill(local_bucket_ptr, local_bucket_ptr + n_buckets, 0UL);

    // Compute bucket counts.
    for (std::uint64_t j = range_beg; j < range_end; ++j) {
      std::uint64_t bucket_id = (((std::uint64_t)tab[j] - lower_bound) >> bucket_range_log);
      ++local_bucket_ptr[bucket_id];
    }

    // Compute destination pointers.
    #pragma omp barrier
    #pragma omp single
    {
      std::uint64_t total_buckets_size = 0;
      for (std::uint64_t bucket_id = 0; bucket_id < n_buckets; ++bucket_id) {
        std::uint64_t this_bucket_size = 0;
        for (std::uint64_t i = 0; i < n_threads; ++i) {
          std::uint64_t local_bucket_size = bucket_ptr[i][bucket_id];
          bucket_ptr[i][bucket_id] = total_buckets_size + this_bucket_size;
          this_bucket_size += local_bucket_size;
        }
        total_buckets_size += this_bucket_size;
      }
    }

    // Move items into buckets.
    for (std::uint64_t j = range_beg; j < range_end; ++j) {
      std::uint64_t bucket_id = ((temp[j] - lower_bound) >> bucket_range_log);
      std::uint64_t dest_pos = local_bucket_ptr[bucket_id]++;
      tab[dest_pos] = temp[j];
    }
  }

  // Free the memory for bucket_ptr. Keep only bucket sizes.
  std::vector<std::uint64_t> unrefined_bucket_sizes(n_buckets);
  for (std::uint64_t bucket_id = 0; bucket_id < n_buckets; ++bucket_id) {
    unrefined_bucket_sizes[bucket_id] = bucket_ptr[n_threads - 1][bucket_id];
    if (bucket_id > 0)
      unrefined_bucket_sizes[bucket_id] -= bucket_ptr[n_threads - 1][bucket_id - 1];
  }
  for (std::uint64_t thread_id = 0; thread_id < n_threads; ++thread_id)
    delete[] bucket_ptr[thread_id];
  delete[] bucket_ptr;

  // Compute the output. If necessary, refine large buckets recursively.
  std::uint64_t cur_bucket_beg = 0;
  for (std::uint64_t bucket_id = 0; bucket_id < n_buckets; ++bucket_id) {
    if (unrefined_bucket_sizes[bucket_id] > max_bucket_size) {
      std::uint64_t lower_bound_rec = lower_bound + bucket_id * bucket_range;
      std::uint64_t upper_bound_rec = std::min(lower_bound_rec + bucket_range, upper_bound);
      permute_into_small_buckets(tab + cur_bucket_beg, temp, unrefined_bucket_sizes[bucket_id],
          lower_bound_rec, upper_bound_rec, max_bucket_size, output_bucket_sizes);
    } else output_bucket_sizes.push_back(unrefined_bucket_sizes[bucket_id]);
    cur_bucket_beg += unrefined_bucket_sizes[bucket_id];
  }
}

template<typename int_type>
void set_bits(std::uint64_t *bv,
    std::uint64_t bv_size, int_type *tab,
    std::uint64_t tab_size, int_type *temp) {
  std::uint64_t max_threads = omp_get_max_threads();

  // Partition the input array into buckets.
  std::vector<std::uint64_t> bucket_sizes;
  {
    // First, partition the input array into small buckets. There may be
    // a lot of them, so they need to be merged into larger buckets.
    std::vector<std::uint64_t> small_bucket_sizes;
    std::uint64_t ideal_bucket_size = std::max(512UL, (tab_size + max_threads - 1) / max_threads);
    std::uint64_t max_bucket_size = 2UL * ideal_bucket_size;
    permute_into_small_buckets(tab, temp, tab_size, 0,
        bv_size, max_bucket_size, small_bucket_sizes);

    // Merge small buckets into at most max_threads final buckets.
    std::uint64_t n_small_buckets = small_bucket_sizes.size();
    std::uint64_t small_bucket_ptr = 0;
    for (std::uint64_t bucket_id = 0; bucket_id < max_threads; ++bucket_id) {
      if (small_bucket_ptr < n_small_buckets) {
        std::uint64_t cur_bucket_total_size = small_bucket_sizes[small_bucket_ptr++];
        std::uint64_t cur_bucket_range_end = small_bucket_ptr;

        // Keep adding buckets as long as we are
        // getting closer to the ideal bucket size.
        while (cur_bucket_range_end < n_small_buckets && (std::abs((std::int64_t)(cur_bucket_total_size +
                small_bucket_sizes[cur_bucket_range_end]) - (std::int64_t)ideal_bucket_size) <=
                std::abs((std::int64_t)cur_bucket_total_size - (std::int64_t)ideal_bucket_size) ||
                (bucket_id + 1 == max_threads)))
          cur_bucket_total_size += small_bucket_sizes[cur_bucket_range_end++];

        // Add the final bucket to the list.
        bucket_sizes.push_back(cur_bucket_total_size);
        small_bucket_ptr = cur_bucket_range_end;
      }
    }
  }

  // Update the bits in bv. The above partitioning guarantees that
  // no thread will attempt to update the same word in bv and
  // that all threads will update roughly the same amount of bits.
  // Lastly, the above guarantees bucket_sizes.size() <= max_threads.
  {
    // Partial (exclusive) sum over bucket_sizes.
    std::uint64_t total_bucket_size = 0;
    std::uint64_t n_buckets = bucket_sizes.size();
    for (std::uint64_t bucket_id = 0; bucket_id < n_buckets; ++bucket_id) {
      std::uint64_t this_bucket_size = bucket_sizes[bucket_id];
      bucket_sizes[bucket_id] = total_bucket_size;
      total_bucket_size += this_bucket_size;
    }

    // Set the bits in the bitvector.
    #pragma omp parallel num_threads(n_buckets)
    {
      std::uint64_t bucket_id = omp_get_thread_num();
      std::uint64_t bucket_beg = bucket_sizes[bucket_id];
      std::uint64_t bucket_end = (bucket_id + 1 == n_buckets) ? total_bucket_size : bucket_sizes[bucket_id + 1];
      for (std::uint64_t j = bucket_beg; j < bucket_end; ++j) {
        std::uint64_t bv_idx = tab[j];
        bv[bv_idx >> 6] |= (1UL << (bv_idx & 63));
      }
    }
  }
}
#endif

}  // namespace em_succinct_irreducible_private

#endif  // __EM_SUCCINCT_IRREDUCIBLE_SRC_SET_BITS_HPP_INCLUDED
