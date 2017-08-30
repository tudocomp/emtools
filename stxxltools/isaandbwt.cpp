/***************************************************************************
 *  test1.cpp
 *
 *  Part of a simple STXXL example. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iostream>
#include <limits>

#include <stxxl/vector>
#include <stxxl/random>
#include <stxxl/sort>
#include <stxxl/bits/algo/ksort.h>
//#include "/scripts/code/dcheck.hpp"

#include <string>
#include <sstream>

#ifndef DCHECK
#define DCHECK_(x,y,z) \
  if (!(x)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x) + ", we got " + std::to_string(y) + " vs " + std::to_string(z))
#define DCHECK(x) \
  if (!(x)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x))
#define DCHECK_EQ(x, y) DCHECK_((x) == (y), x,y)
#define DCHECK_NE(x, y) DCHECK_((x) != (y), x,y)
#define DCHECK_LE(x, y) DCHECK_((x) <= (y), x,y)
#define DCHECK_LT(x, y) DCHECK_((x) < (y) ,x,y)
#define DCHECK_GE(x, y) DCHECK_((x) >= (y),x,y )
#define DCHECK_GT(x, y) DCHECK_((x) > (y) ,x,y)
#endif //DCHECK





// struct my_less_int : std::less<int>
// {
//     int min_value() const { return std::numeric_limits<int>::min(); };
//     int max_value() const { return std::numeric_limits<int>::max(); };
// };
//
// int main(int argv,)
// {
//     // create vector
//     stxxl::VECTOR_GENERATOR<int>::result vector;
//
//     // fill vector with random integers
//     stxxl::random_number32 random;
//
//     for (size_t i = 0; i < 100*1024*1024; ++i) {
//         vector.push_back(random());
//     }
//
//     // sort vector using 16 MiB RAM
//     stxxl::sort(vector.begin(), vector.end(), my_less_int(), 16*1024*1024);
//
//     // output first and last items:
//     std::cout << vector.size() << " items sorted ranging from "
//               << vector.front() << " to " << vector.back() << std::endl;
//
//     return 0;
// }



#include <iostream>
#include <fstream>
#include <stxxl/bits/common/uint_types.h>

size_t filesize( const char*const filepath ){
	std::ifstream file(filepath, std::ios::binary | std::ios::ate | std::ios::in);
	if(!file.good()) return 0;
	return file.tellg();
}
bool file_exists(const char *const filepath) {
	std::ifstream infile(filepath);
		return infile.good();
}


template<class int_t>
class IntegerFileForwardIterator {
	const size_t m_size;
	std::ifstream m_is;
	size_t m_index;
	char m_buf[sizeof(int_t)];
	public:

	IntegerFileForwardIterator(const char*const filename) 
	: m_size { filesize(filename) }
	, m_is {filename, std::ios::binary | std::ios::in }
	, m_index {0}
	{}

	size_t size() const { return m_size/sizeof(int_t); }
	size_t index() const { return m_index; }
	int_t operator*() { return *reinterpret_cast<int_t*>(m_buf); }
	IntegerFileForwardIterator& operator++(int) { 
		m_is.read(m_buf, sizeof(int_t));
		++m_index;
		return *this;
	}
};

template<class int_t>
class IntegerFileArray {
	const size_t m_size;
	std::ifstream m_is;
	public:
	IntegerFileArray(const char*const filename) 
	: m_size { filesize(filename) }
	, m_is {filename, std::ios::binary | std::ios::in }
	{}
	int_t operator[](size_t i) {
//		DCHECK_LT(i, size());
		m_is.seekg(i*sizeof(int_t), std::ios_base::beg);
		char buf[sizeof(int_t)];
		m_is.read(buf, sizeof(int_t));
		return *reinterpret_cast<int_t*>(buf);
	}
	size_t size() const { return m_size/sizeof(int_t); }
};


using namespace stxxl;
void bwt() {
	IntegerFileForwardIterator<uint40> sa   { "/bighome/workspace/eSAIS/build/src/a.sa5" };
	std::ifstream is("/bighome/workspace/eSAIS/build/src/a", std::ios::binary | std::ios::in);
	while(is) {

	}
	
}

		// struct KeyExtractor {
		// 	typedef uint40 key_type;
		// 	typedef std::pair<uint40,uint40> value_type;
		// 	key_type m_key;
        //
		// 	KeyExtractor() {}
		// 	KeyExtractor(const key_type& k) : m_key(k) {}
		// 	key_type operator()(const value_type& v) const { return v.first; }
		// 	value_type min_value() const { return value_type(0,0); }
		// 	//value_type max_value() const { return std::make_pair(std::numeric_limits<key_type>::max(),0); }
		// 	value_type max_value() const { return value_type(m_key,0); }
		// };

template<class pair_t>
		struct KeyExtractor {
			typedef pair_t value_type;
			typedef typename pair_t::first_type key_type;
			key_type m_key;

			KeyExtractor() {}
			KeyExtractor(const key_type& k) : m_key(k) {}
			key_type operator()(const value_type& v) const { return v.first; }
			value_type min_value() const { return pair_t(0, (typename value_type::second_type)0); }
			//value_type max_value() const { return std::make_pair(std::numeric_limits<key_type>::max(),0); }
			value_type max_value() const { return pair_t(m_key, (typename value_type::second_type)0); }
		};
		
using namespace std;
int main(int argc, char** argv) {
	if(argc != 2) {
		cout << "Usage: " << argv[0] << " text-file" << std::endl;
		return 1;
	}
	const std::string textfilename = argv[1];
	const std::string safilename = textfilename  + ".sa5";
	const std::string isafilename = textfilename  + ".isa5";
	const std::string bwtfilename = textfilename  + ".bwt";
	if(!file_exists(textfilename.c_str())) {
		cout << "Could not open text file " << textfilename << std::endl;
		return 1;
	}
	if(!file_exists(safilename.c_str())) {
		cout << "Could not open SA file " << safilename << std::endl;
		return 1;
	}
    stxxl::VECTOR_GENERATOR<std::pair<uint40,uint40>>::result isa; // (text_position, factor_length)
	IntegerFileForwardIterator<uint40> safile { safilename.c_str() };
	while(safile.index() < safile.size()) {
		uint40 index = static_cast<uint64>(safile.index());
		isa.push_back(std::make_pair(*safile++,index));
	}
	stxxl::ksort(isa.begin(), isa.end(), KeyExtractor<std::pair<uint40,uint40>>(isa.size()),512*1024*1024); //, STXXL_DEFAULT_ALLOC_STRATEGY());
	std::ofstream isa_out(isafilename, std::ios::binary);
	for(auto it = isa.begin(); it != isa.end(); ++it) {
		isa_out.write((char*)(&it->second), sizeof(uint40));
	}
	isa_out.close();
	
    stxxl::VECTOR_GENERATOR<std::pair<uint40,char>>::result bwt; 
	ifstream textfile(textfilename, ios::in | ios::binary);
	const uint40 isa_zero = isa.begin()->second;
	{
		auto it = isa.begin();
		++it;
		for(; it != isa.end(); ++it) {
			bwt.push_back(std::make_pair(it->second+1, textfile.get()));
		DCHECK(textfile.good());
		}
		DCHECK(textfile.good());
		bwt.push_back(std::make_pair(isa_zero,textfile.get()));
		DCHECK(textfile.good());
	}
	stxxl::ksort(bwt.begin(), bwt.end(), KeyExtractor<std::pair<uint40,char>>(isa.size()),512*1024*1024); //, STXXL_DEFAULT_ALLOC_STRATEGY());
	std::ofstream bwt_out(bwtfilename, std::ios::binary);
	for(auto it = bwt.begin(); it != bwt.end(); ++it) {
		// if(it->second == 0) bwt_out.put(1); // TODO BUG: prevent writing the 0-byte by writing 1
		bwt_out.put(it->second);
	}


    // stxxl::VECTOR_GENERATOR<std::pair<uint40,uint40>>::result bwt; // (text_position, factor_length)
	// ifstream textfile(textfilename, ios::in | ios::binary);
	
	//
	// {
	// 	IntegerFileForwardIterator<uint40> sa   { "/bighome/workspace/eSAIS/build/src/a.sa5" };
	// 	IntegerFileForwardIterator<uint40> isa  { "/bighome/workspace/eSAIS/build/src/a.isa5" };
	// 	IntegerFileForwardIterator<uint40> plcp { "/bighome/workspace/eSAIS/build/src/a.plcp5" };
	// 	while(sa.index() < sa.size()) {
	// 		std::cout << *sa++ << "," << *isa++ << "," << *plcp++ << endl; 
	// 	}
	// }
	// std::cout << endl;
	// {
	// 	IntegerFileArray<uint40> sa   { "/bighome/workspace/eSAIS/build/src/a.sa5" };
	// 	IntegerFileArray<uint40> isa  { "/bighome/workspace/eSAIS/build/src/a.isa5" };
	// 	IntegerFileArray<uint40> plcp { "/bighome/workspace/eSAIS/build/src/a.plcp5" };
	// 	for(size_t i = 0; i < sa.size(); ++i) {
	// 		std::cout << sa[i] << "," << isa[i] << "," << plcp[i] << endl; 
	// 	}
	// }
}
