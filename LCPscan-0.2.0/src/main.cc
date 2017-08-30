/**
 * @file    src/main.cc
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

#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <string>
#include <fstream>

#include "lcpscan_src/lcpscan.h"


// Main settings.
const long kDiskBlockSize = (1L << 20);  // disk block size
const long kSBlockCnt = 10;              // number of superblocks

char *program_name;

void usage(int status) {
  printf(
"Usage: %s [OPTION]... FILE\n"
"Construct the LCP array for text stored in FILE.\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -d, --distribute=DATA   copy the input into STXXL space prior to computation;\n"
"                          valid choices for DATA are: text, sa, all\n"
"  -h, --help              display this help and exit\n"
"  -m, --mem=MEM           use MEM MiB of RAM for computation (default: 3584)\n"
"                          NOTE: using less than 3584MiB of RAM may result in\n"
"                                multiple passes of external-memory sorting for\n"
"                                very large inputs\n"
"  -n, --nocache           disable disk page caching\n"
"  -o, --output=OUTFILE    specify output file (default: FILE.lcp5)\n"
"  -s, --sa=SUFARRAY       specify the location of the suffix array of FILE\n"
"                          (default: FILE.sa5)\n"
"  -t, --test              test run (do not write output to file)\n",
  program_name);

  std::exit(status);
}

bool file_exists(std::string fname) {
  std::FILE *f = std::fopen(fname.c_str(), "r");
  bool ret = (f != NULL);
  if (f != NULL) std::fclose(f);

  return ret;
}

std::FILE *file_open(std::string fname, std::string mode) {
  std::FILE *f = std::fopen(fname.c_str(), mode.c_str());
  if (!f) {
    std::perror(fname.c_str());
    std::exit(EXIT_FAILURE);
  }

  return f;
}

void file_copy(std::string src_fname, std::string dest_fname) {
  std::ifstream src(src_fname.c_str(), std::ios::binary);
  std::ofstream dst(dest_fname.c_str(), std::ios::binary);
  dst << src.rdbuf();
}

void file_delete(std::string fname) {
  int res = std::remove(fname.c_str());
  if (res) {
    fprintf(stderr, "Failed to delete %s: %s\n",
        fname.c_str(), strerror(errno));
    std::exit(EXIT_FAILURE);
  }
}

std::string absolute_path(std::string fname) {
  char path[1 << 18];
  bool created = false;

  if (!file_exists(fname)) {
    // We need to create the file, since realpath fails on non-existing files.
    std::fclose(file_open(fname, "w"));
    created = true;
  }
  if (!realpath(fname.c_str(), path)) {
    fprintf(stderr, "Error: realpath failed for %s\n", fname.c_str());
    std::exit(EXIT_FAILURE);
  }

  if (created)
    file_delete(fname);

  return std::string(path);
}

void find_stxxl_config() {
  if (file_exists("./.stxxl")) {
    fprintf(stderr, "STXXL config file detected.\n");
    return;
  } else if (file_exists(std::string(std::getenv("HOME")) + "/.stxxl")) {
    fprintf(stderr, "Cannot find STXXL config file. Using $HOME/.stxxl\n");
    std::string src = std::string(std::getenv("HOME")) + "/.stxxl";
    std::string dest = absolute_path("./.stxxl");
    file_copy(src, dest);
    return;
  } else {
    fprintf(stderr, "Error: failed to find/copy STXXL config file!\n");
    std::exit(EXIT_FAILURE);
  }
}

int main(int argc, char **argv) {
  program_name = argv[0];

  static struct option long_options[] = {
    {"help",       no_argument,       NULL, 'h'},
    {"nocache",    no_argument,       NULL, 'n'},
    {"test",       no_argument,       NULL, 't'},
    {"distribute", required_argument, NULL, 'd'},
    {"mem",        required_argument, NULL, 'm'},
    {"output",     required_argument, NULL, 'o'},
    {"sa",         required_argument, NULL, 's'},
    {NULL, 0, NULL, 0}
  };

  bool dry_run = false;
  bool distr_text = false;
  bool distr_sa = false;
  bool distr_flag_set = false;
  bool disable_cache = false;
  long ram_use = 3584L << 20;
  std::string sa_fname("");
  std::string out_fname("");
  std::string distr_what("");

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "hntd:m:o:s:", long_options, NULL)) != -1) {
    switch(c) {
      case 'm':
        ram_use = std::atol(optarg) << 20;
        if (ram_use <= 0L) {
          fprintf(stderr, "Error: invalid RAM limit (%ld)\n\n", ram_use);
          usage(EXIT_FAILURE);
        }
        break;
      case 'o':
        out_fname = std::string(optarg);
        break;
      case 's':
        sa_fname = std::string(optarg);
        break;
      case 'h':
        usage(EXIT_FAILURE); 
      case 't':
        dry_run = true;
        break;
      case 'd':
        distr_flag_set = true;
        distr_what = std::string(optarg);
        break;
      case 'n':
        disable_cache = true;
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
  std::string text_fname = std::string(argv[optind++]);
  if (optind < argc) {
    fprintf(stderr, "Warning: multiple input files provided. "
        "Only the first will be processed.\n\n");
  }

  // Set default SA and output filenames (if not provided).
  if (sa_fname.empty()) sa_fname = text_fname + ".sa5";
  if (out_fname.empty()) out_fname = text_fname + ".lcp5";

  // Check for their existence.
  if (!file_exists(text_fname)) {
    fprintf(stderr, "Error: input file (%s) does not exist\n\n",
        text_fname.c_str());
    usage(EXIT_FAILURE);
  }
  if (!file_exists(sa_fname)) {
    fprintf(stderr, "Error: suffix array (%s) does not exist\n\n",
        sa_fname.c_str());
    usage(EXIT_FAILURE);
  }

  if (distr_flag_set) {
    std::transform(distr_what.begin(), distr_what.end(),
        distr_what.begin(), ::tolower);
    if (distr_what == "text") distr_text = true;
    else if (distr_what == "sa") distr_sa = true;
    else if (distr_what == "all") {
      distr_text = true;
      distr_sa = true;
    } else {
      fprintf(stderr, "Error: invalid argument for -d (%s). Valid "
          "choices are: `text', `sa' or `all'.\n", distr_what.c_str());
      usage(EXIT_FAILURE);
    }
  }

  if (!dry_run && file_exists(out_fname)) {
    // Output file exists -- should we proceed?
    char *line = NULL;
    size_t buflen = 0;
    long len = 0L;

    do {
      printf("Output file (%s) exists. Overwrite? [y/n]: ",
          out_fname.c_str());
      if ((len = getline(&line, &buflen, stdin)) == -1) {
        fprintf(stderr, "\nError: failed to read answer\n\n");
        usage(EXIT_FAILURE);
      }
    } while (len != 2 || (line[0] != 'y' && line[0] != 'n'));

    if (line[0] == 'n') {
      free(line);
      std::exit(EXIT_FAILURE);
    }
    free(line);
  }
  
  find_stxxl_config();

  construct_lcp<kDiskBlockSize>(text_fname, sa_fname, out_fname,
      ram_use, kSBlockCnt, dry_run, distr_text, distr_sa, disable_cache);
}
