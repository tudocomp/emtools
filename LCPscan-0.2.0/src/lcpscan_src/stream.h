/**
 * @file    src/lcpscan_src/stream.h
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

#ifndef __LCPSCAN_SRC_STREAM_H_INCLUDED
#define __LCPSCAN_SRC_STREAM_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <string>

#include "./utils.h"


namespace lcpscan_private {

template<typename T>
struct stream_reader {
  stream_reader(std::string fname, int buf_bytes = (4L << 20))
      : m_bufelems(buf_bytes / sizeof(T)),
        m_filled(0),
        m_offset(0L) {
    m_file = utils::file_open(fname, "r");
    m_buffer = new T[m_bufelems];
    refill();
  }

  ~stream_reader() {
    delete[] m_buffer;
    std::fclose(m_file);
  }

  inline T peek() {
    return m_buffer[m_pos];
  }

  inline T& operator[] (long i) {
    assert(i >= m_offset);
    if (i >= m_offset + m_filled) {  // elem not in a buffer
      if (i < m_offset + m_filled + m_bufelems) {
        // elem in the next buffer, no need to fseek
        refill();
      } else {
        // elem potentially far away, fseek
        refill(i);
      }
    }

    return m_buffer[i - m_offset];
  }

  stream_reader& operator++ () {
    ++m_pos;
    if (m_pos == m_filled) refill();

    return *this;
  }

  inline T read() {
    T ret = m_buffer[m_pos++];
    if (m_pos == m_filled) refill();

    return ret;
  }

  inline bool empty() {
    return (!m_filled && !refill());
  }

 private:
  int refill(long new_offset = -1) {
    if (new_offset != -1) {
      m_offset = new_offset;
      std::fseek(m_file, m_offset, SEEK_SET);
    } else {
      m_offset += m_filled;
    }

    m_filled = (int)std::fread(m_buffer, sizeof(T), m_bufelems, m_file);
    m_pos = 0;

    return m_filled;
  }

  std::FILE *m_file;

  int m_bufelems, m_filled, m_pos;
  long m_offset;
  T *m_buffer;
};

template<typename T>
struct stream_writer {
  stream_writer(std::string fname, int bufsize = (4 << 20))
      : m_bufelems(bufsize / sizeof(T)),
        m_filled(0) {
    m_file = utils::file_open(fname.c_str(), "w");
    m_buffer = new T[m_bufelems];
  }

  ~stream_writer() {
    if (m_filled) flush();
    delete[] m_buffer;
    std::fclose(m_file);
  }

  inline void write(T x) {
    m_buffer[m_filled++] = x;
    if (m_filled == m_bufelems) flush();
  }

 private:
  void flush() {
    utils::write_to_file(m_buffer, m_filled, m_file);
    m_filled = 0;
  }

  std::FILE *m_file;

  int m_bufelems, m_filled;
  T *m_buffer;
};

}  // namespace lcpscan_private

#endif  // __LCPSCAN_SRC_STREAM_H_INCLUDED
