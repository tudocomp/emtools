/**
 * @file    main.cpp
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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <string>
#include <getopt.h>
#include <unistd.h>

#include "uint40.hpp"
#include "uint48.hpp"
#include "em_succinct_irreducible_src/compute_plcp_bitvector.hpp"

char *program_name;

void usage(int status) {
  printf(

"Usage: %s [OPTION]... FILE\n"
"Construct the PLCP array (bitvector representation) for text stored in FILE.\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -b, --bwt=BWTFILE       specify the location of the Burrows-Wheeler\n"
"                          transform of FILE (default: FILE.bwt)\n"
"  -h, --help              display this help and exit\n"
"  -i, --intsize=SIZE      use integers of SIZE bytes (default: 5). Currently\n"
"                          supported values are 4, 5, 6, and 8\n"
"  -m, --mem=MEM           use MEM MiB of RAM for computation (default: 3584)\n"
"  -o, --output=OUTFILE    specify output filename (default: FILE.plcp)\n"
"  -s, --sa=SUFARRAY       specify the location of the suffix array of FILE\n"
"                          (default: FILE.saX, X = integer size, see -i flag)\n",
    program_name);

  std::exit(status);
}

bool file_exists(std::string filename) {
  std::FILE *f = std::fopen(filename.c_str(), "r");
  bool ret = (f != NULL);
  if (f != NULL) std::fclose(f);

  return ret;
}

std::FILE *file_open(std::string filename, std::string mode) {
  std::FILE *f = std::fopen(filename.c_str(), mode.c_str());
  if (f == NULL) {
    std::perror(filename.c_str());
    std::exit(EXIT_FAILURE);
  }
  return f;
}

std::uint64_t file_size(std::string filename) {
  std::FILE *f = file_open(filename, "r");
  std::fseek(f, 0, SEEK_END);
  long size = std::ftell(f);
  if (size < 0) {
    std::perror(filename.c_str());
    std::exit(EXIT_FAILURE);
  }
  std::fclose(f);
  return (std::uint64_t)size;
}

template<typename int_type>
std::string intToStr(int_type x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

template<typename text_offset_type>
void compute_plcp_bitvector(std::string text_filename, std::string sa_filename,
    std::string bwt_filename, std::string output_filename, std::uint64_t ram_use) {
  std::uint64_t text_length = file_size(text_filename);
  if (2UL * text_length <= std::numeric_limits<text_offset_type>::max()) {
    em_succinct_irreducible_private::compute_plcp_bitvector<text_offset_type, text_offset_type>(text_filename, sa_filename, bwt_filename, output_filename, ram_use);
  } else {
    if (sizeof(text_offset_type) < 4) em_succinct_irreducible_private::compute_plcp_bitvector<text_offset_type, std::uint32_t>(text_filename, sa_filename, bwt_filename, output_filename, ram_use);
    else if (sizeof(text_offset_type) == 4) em_succinct_irreducible_private::compute_plcp_bitvector<text_offset_type, uint40>(text_filename, sa_filename, bwt_filename, output_filename, ram_use);
    else if (sizeof(text_offset_type) == 5) em_succinct_irreducible_private::compute_plcp_bitvector<text_offset_type, uint48>(text_filename, sa_filename, bwt_filename, output_filename, ram_use);
    else em_succinct_irreducible_private::compute_plcp_bitvector<text_offset_type, std::uint64_t>(text_filename, sa_filename, bwt_filename, output_filename, ram_use);
  }
}

int main(int argc, char **argv) {
  srand(time(0) + getpid());
  program_name = argv[0];

  static struct option long_options[] = {
    {"bwt",     required_argument, NULL, 'b'},
    {"help",    no_argument,       NULL, 'h'},
    {"intsize", required_argument, NULL, 'i'},
    {"mem",     required_argument, NULL, 'm'},
    {"output",  required_argument, NULL, 'o'},
    {"sa",      required_argument, NULL, 's'},
    {NULL, 0, NULL, 0}
  };

  std::uint64_t int_size = 5;
  std::uint64_t ram_use = 3584UL << 20;
  std::string out_filename("");
  std::string sa_filename("");
  std::string bwt_filename("");

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "b:hi:m:o:s:", long_options, NULL)) != -1) {
    switch(c) {
      case 'b':
        bwt_filename = std::string(optarg);
        break;
      case 'h':
        usage(EXIT_FAILURE);
      case 'i':
        int_size = std::atol(optarg);
        if (!(int_size == 4 || int_size == 5 || int_size == 6 || int_size == 8)) {
          fprintf(stderr, "Error: invalid int size (%lu)\n\n", int_size);
          usage(EXIT_FAILURE);
        }
        break;
      case 'm':
        ram_use = std::atol(optarg) << 20;
        if (ram_use == 0) {
          fprintf(stderr, "Error: invalid RAM limit (%lu)\n\n", ram_use);
          usage(EXIT_FAILURE);
        }
        break;
      case 'o':
        out_filename = std::string(optarg);
        break;
      case 's':
        sa_filename = std::string(optarg);
        break;
      default:
        usage(EXIT_FAILURE);
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Error: FILE not provided\n\n");
    usage(EXIT_FAILURE);
  }

  // Parse the text filename.
  std::string text_filename = std::string(argv[optind++]);
  if (optind < argc) {
    fprintf(stderr, "Warning: multiple input files provided. "
    "Only the first will be processed.\n");
  }

  // Set default filenames (if not provided).
  if (sa_filename.empty()) sa_filename = text_filename + ".sa" + intToStr(int_size);
  if (out_filename.empty()) out_filename = text_filename + ".plcp";
  if (bwt_filename.empty()) bwt_filename = text_filename + ".bwt";

  // Check if input text, suffix array, and BWT exist.
  if (!file_exists(text_filename)) {
    fprintf(stderr, "Error: input file (%s) does not exist\n\n",
        text_filename.c_str());
    usage(EXIT_FAILURE);
  }
  if (!file_exists(sa_filename)) {
    fprintf(stderr, "Error: suffix array (%s) does not exist\n\n",
        sa_filename.c_str());
    usage(EXIT_FAILURE);
  }
  if (!file_exists(bwt_filename)) {
    fprintf(stderr, "Error: BWT of input text (%s) does not exist\n\n",
        bwt_filename.c_str());
    usage(EXIT_FAILURE);
  }

  if (file_exists(out_filename)) {
    // Output file exists, should we proceed?
    char *line = NULL;
    std::uint64_t buflen = 0;
    std::int64_t len = 0L;

    do {
      printf("Output file (%s) exists. Overwrite? [y/n]: ",
          out_filename.c_str());
      if ((len = getline(&line, &buflen, stdin)) == -1) {
        printf("\nError: failed to read answer\n\n");
        std::fflush(stdout);
        usage(EXIT_FAILURE);
      }
    } while (len != 2 || (line[0] != 'y' && line[0] != 'n'));

    if (line[0] == 'n') {
      free(line);
      std::exit(EXIT_FAILURE);
    }
    free(line);
  }

  // Run the algorithm.
  if (int_size == 4) compute_plcp_bitvector<std::uint32_t>(text_filename, sa_filename, bwt_filename, out_filename, ram_use);
  else if (int_size == 5) compute_plcp_bitvector<uint40>(text_filename, sa_filename, bwt_filename, out_filename, ram_use);
  else if (int_size == 6) compute_plcp_bitvector<uint48>(text_filename, sa_filename, bwt_filename, out_filename, ram_use);
  else compute_plcp_bitvector<std::uint64_t>(text_filename, sa_filename, bwt_filename, out_filename, ram_use);
}
