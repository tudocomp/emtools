/* example.cpp */

#include <cstdio>
#include <iomanip>
#include <vector>
#include <iostream>
#include <climits>

#include <stxxl/io>
#include <stxxl/aligned_alloc>

#include <stxxl/mng>
#include <stxxl/sort>
#include <stxxl/vector>

#include <stdlib.h>
#include <string.h>

#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <limits.h>
#include <inttypes.h>
#include <math.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <dirent.h>
#include <getopt.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <stack>
#include <set>
#include <map>
// #include <ext/algorithm>
// #include <parallel/algorithm>
// #include <ext/pb_ds/assoc_container.hpp>

#include <stxxl/bits/algo/sort.h>
#include <stxxl/bits/containers/priority_queue.h>
#include <stxxl/bits/containers/vector.h>
#include <stxxl/bits/containers/queue.h>
#include <stxxl/bits/containers/stack.h>
#include <stxxl/bits/containers/sorter.h>
#include <stxxl/bits/containers/deque2.h>
#include <stxxl/bits/io/iostats.h>
#include <stxxl/bits/mng/buf_istream.h>
#include <stxxl/bits/mng/buf_ostream.h>
#include <stxxl/bits/stream/sort_stream.h>
#include <stxxl/bits/stream/stream.h>


struct counter_object
{
    // This stream produces a sequence of integers.
    typedef int         value_type;

private:
    // A class attribute to save the current value.
    int                 m_current_value;

public:
    // A constructor to set the initial value to 1.
    counter_object()
        : m_current_value(1)
    {
    }

    // The retrieve operator returning the current value.
    const value_type& operator* () const
    {
        return m_current_value;
    }

    // Increment operator advancing to the next integer.
    counter_object& operator++ ()
    {
        ++m_current_value;
        return *this;
    }

    // Empty indicator, which in this case can check the current value.
    bool empty() const
    {
        return (m_current_value > 1000);
    }
};


// define comparator class: compare right-most decimal and then absolute value
struct CompareMod10
{
    // comparison operator() returning true if (a < b)
    inline bool operator() (int a, int b) const
    {
        if ((a % 10) == (b % 10))
            return a < b;
        else
            return (a % 10) < (b % 10);
    }

    // smallest possible integer value
    int min_value() const { return INT_MIN; }
    // largest possible integer value
    int max_value() const { return INT_MAX; }
};

int main()
{
  static const int ram_use = 10*1024*1024;   // amount of memory to use in runs creation

  // define a runs sorter which accepts imperative push()s and orders by CompareMod10 object.
  typedef stxxl::sorter<int, CompareMod10> sr_counter_type;

  // instance of CompareMod10 comparator class.
  CompareMod10    comparemod10;

  // instance of sorter which waits for input.
  sr_counter_type sr_counter (comparemod10, ram_use);

  // write sequence of integers into sorter, which creates sorted runs
  for (int i = 1; i <= 1000; ++i)
    sr_counter.push(i);

  // signal sorter that the input stream is finished and switch to output mode.
  sr_counter.sort();

  // read sorted stream: sorter also conforms to the stream interface.
  while (!sr_counter.empty())
  {
      std::cout << *sr_counter << " ";
      ++sr_counter;
  }
  std::cout << std::endl;

  return 0;
}

