#include <iostream>
#include <limits>

#include <stxxl/vector>
#include <stxxl/sorter>
#include <stxxl/bits/common/uint_types.h>
#include <stxxl/io>

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

// Select IO System to be used to read/write input/output
// LinuxAIO is fasted but -suprise- only available on Linux machines
using FileDriver = stxxl::linuxaio_file;
constexpr auto FileModeRead = stxxl::file::RDONLY | stxxl::file::DIRECT | stxxl::file::NO_LOCK;
constexpr auto FileModeWrite = stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT | stxxl::file::TRUNC;


// Configure Integer type used for read/writing input/output
using FileInt = stxxl::uint40;
using IntPair = std::pair<FileInt, FileInt>;
using IntCharPair = std::pair<FileInt, char>;

using IntVector = typename stxxl::VECTOR_GENERATOR<FileInt, 4, 8, sizeof(FileInt)<<19>::result;
using CharVector = typename stxxl::VECTOR_GENERATOR<char>::result;


static bool file_exists(const char *const filepath) {
    std::ifstream infile(filepath);
    return infile.good();
}

template<typename pair_t, bool OnlyFirst = false>
struct PairCompare {
    using first_type = typename pair_t::first_type;
    using second_type = typename pair_t::second_type;

    bool operator() (const pair_t& a, const pair_t& b) const {
        return OnlyFirst ? (a.first < b.first) : (a < b);
    }

    pair_t min_value() const {
      return {std::numeric_limits<first_type >::min(),
              std::numeric_limits<second_type>::min()};
    }

    pair_t max_value() const {
      return {std::numeric_limits<first_type >::max(),
            std::numeric_limits<second_type>::max()};
    }
};

int main(int argc, char* argv[]) {
    const size_t sorter_mem = 512llu << 20;

    if(argc != 2) {
        std::cout << "Usage: " << argv[0] << " text-file" << std::endl;
        return 1;
    }
    const std::string textfilename = argv[1];
    const std::string safilename = textfilename  + ".sa5";
    const std::string isafilename = textfilename  + ".isa5";
    const std::string bwtfilename = textfilename  + ".bwt";

    if(!file_exists(textfilename.c_str())) {
        std::cout << "Could not open text file " << textfilename << std::endl;
        return 1;
    }
    if(!file_exists(safilename.c_str())) {
        std::cout << "Could not open SA file " << safilename << std::endl;
        return 1;
    }


    auto& stats = *stxxl::stats::get_instance();
    stxxl::stats_data stats_begin(stats);

    // Compute ISA
    stxxl::sorter<IntPair, PairCompare<IntPair, true> > isa_sorter(
      PairCompare<IntPair, true>(), sorter_mem);

    {
        // Open SA File for input
        FileDriver safile(safilename, FileModeRead);
        IntVector savector(&safile);

        IntVector::bufreader_type reader(savector);
        for(stxxl::unsigned_type index = 0; !reader.empty(); ++reader, ++index) {
            isa_sorter.push( std::make_pair(*reader, index) );
        }

        isa_sorter.sort();
    }


    // Write ISA file and fill BWT file
    stxxl::sorter<IntCharPair, PairCompare<IntCharPair, true> >
        bwt_sorter(PairCompare<IntCharPair, true>(), sorter_mem);
    {
        // open isa output file
        FileDriver isa_file(isafilename, FileModeWrite);
        IntVector isa_vector(&isa_file);
        isa_vector.resize(isa_sorter.size());
        typename IntVector::bufwriter_type isa_writer(isa_vector);

        // open text file for reading
        FileDriver text_file(textfilename, FileModeRead);
        CharVector text_vector(&text_file);
        typename CharVector::bufreader_type text_reader(text_vector);
        DCHECK_EQ(text_vector.size(), isa_sorter.size());

        // the first entry in isa needs special treatment ...
        const auto isa_zero = isa_sorter->second;
        isa_writer << isa_zero;


        for(++isa_sorter; !isa_sorter.empty(); ++isa_sorter, ++text_reader) {
            // write isa output file
            isa_writer << isa_sorter->second;

            // fill sorter for next stage
            bwt_sorter.push(IntCharPair{isa_sorter->second, *text_reader});
        }

        bwt_sorter.push(IntCharPair{isa_zero, *text_reader});

        isa_sorter.finish_clear();
        isa_writer.finish();
    }

    // Sort and write out BWT
    {
        // open bwt output file
        FileDriver bwt_file(bwtfilename, FileModeWrite);
        CharVector bwt_vector(&bwt_file);
        bwt_vector.resize(bwt_sorter.size());
        typename CharVector::bufwriter_type bwt_writer(bwt_vector);

        bwt_sorter.sort();

        for(; !bwt_sorter.empty(); ++bwt_sorter)
            bwt_writer << bwt_sorter->second;

        bwt_writer.finish();
    }

    stxxl::stats_data stats_end(stats);
    std::cout << (stats_end - stats_begin) << std::endl;

    return 0;
}
