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

#include <iostream>
#include <fstream>
#include <stxxl/bits/common/uint_types.h>
//#include "/home/niki/opt/include/sdsl/bits.hpp"
#include <sdsl/bits.hpp>

size_t filesize( const char*const filepath ){
	std::ifstream file(filepath, std::ios::binary | std::ios::ate | std::ios::in);
	if(!file.good()) return 0;
	return file.tellg();
}
bool file_exists(const char *const filepath) {
	std::ifstream infile(filepath);
		return infile.good();
}
typedef size_t len_t;




class PLCPFileForwardIterator {
	std::ifstream m_is;
	
	uint64_t m_chunk = 0; // current data chunk
	len_t m_idx = 0; // current select parameter
	len_t m_block = 0; // block index
	len_t m_blockrank = 0; //number of ones up to previous block
	uint_fast8_t m_ones; // number of ones in the current block `m_block`

	void read_chunk() {
		m_is.read(reinterpret_cast<char*>(&m_chunk), sizeof(decltype(m_chunk)));
		m_ones = sdsl::bits::cnt(m_chunk);
	}

	public:
	static constexpr const len_t eof = -1;
	PLCPFileForwardIterator(const char* filepath) 
		: m_is(filepath) 
	{
		read_chunk();
	}

	len_t index() const { return m_idx; }
	bool has_next() const {
		return m_is.good();
	}

	len_t next_select() {
		while(m_blockrank+m_ones < m_idx+1) {
			if(!m_is) {break;}
			++m_block;
			m_blockrank += m_ones;
			read_chunk();
		}
		return 64*m_block + sdsl::bits::sel(m_chunk, m_idx+1-m_blockrank);
	}
	len_t operator()() {
		const len_t ret = next_select() - 2*m_idx;
		return ret;
	}
	void advance() {
		++m_idx;
	}
};

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


using namespace std;
int main(int argc, char** argv) {
	if(argc != 2) {
		cout << "Usage: " << argv[0] << " text-file" << std::endl;
		return 1;
	}
	const std::string textfilename = argv[1];
	const std::string plcpfilename = textfilename  + ".plcp";
	if(!file_exists(plcpfilename.c_str())) {
		cout << "Could not open text file " << textfilename << std::endl;
		return 1;
	}
	PLCPFileForwardIterator p(plcpfilename.c_str());
	size_t i=0;
//	while(i < 100) {
	while(p.has_next()) {
		size_t entry = p();
//		if(!p.has_next()) break;
		std::cout << i++ << "->" << entry << "->" << p.has_next() << endl;
		p.advance();
	}


}
