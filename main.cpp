/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <iostream>
#include <fstream>
#include <getopt.h>

#include "vg/vg.hpp"
#include "Variant.h"

#include "snpmerge.h"

using namespace vcflib;
using namespace vg;
using namespace std;

static const int DefaultWindowSize = 100;

void help_main(char** argv)
{
  cerr << "usage: " << argv[0] << " [options] VGFILE VCFFILE" << endl
       << "Pull apart adjacent snps when genotype information permits in"
       << " order to reduce number of paths that do not reflect haplotypes."
       << "\nThe input vg file must have been created from the input vcf file."
       << endl
       << "options:" << endl
       << "    -h, --help          print this help message" << endl
       << "    -w, --window-size N maximum distance between adjacent snps to be"
       << " merged (default=" << DefaultWindowSize << ")" << endl
       << "    -o, --offset N      vcf-coordinate of first position in vg path"
       << " (default=1)" << endl;
}

int main(int argc, char** argv) {
    
  if(argc == 1) {
    // Print the help
    help_main(argv);
    return 1;
  }

  int windowSize = DefaultWindowSize;
  int offset = 1;
    
  optind = 1; // Start at first real argument
  bool optionsRemaining = true;
  while(optionsRemaining) {
    static struct option longOptions[] = {
      {"window-size", required_argument, 0, 'w'},
      {"offset", required_argument, 0, 'o'},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}
    };

    int optionIndex = 0;

    switch(getopt_long(argc, argv, "w:o:h", longOptions, &optionIndex)) {
      // Option value is in global optarg
    case -1:
      optionsRemaining = false;
      break;
    case 'w': 
      windowSize = atol(optarg);
      break;
    case 'o':
      offset = atol(optarg);
      break;
    case 'h': // When the user asks for help
      help_main(argv);
      exit(1);
      break;
    default:
      cerr << "Illegal option" << endl;
      exit(1);
    }
  }

  if(argc - optind < 2) {
    // We don't have two positional arguments
    // Print the help
    help_main(argv);
    return 1;
  }

  string vgFile = argv[optind++];
  string vcfFile = argv[optind++]; 
    
  // Open the vg file
  ifstream vgStream(vgFile);
  if(!vgStream.good())
  {
    cerr << "Could not read " << vgFile << endl;
    exit(1);
  }
  VG vg(vgStream);
  cerr << "vg open" << endl;

  // Open the vcf file
  VariantCallFile vcf;
  vcf.open(vcfFile);
  cerr << "vcf open" << endl;

  SNPMerge snpMerge;

  // Process all adjacant variants my merging them in the graph
  // when possible
  snpMerge.processGraph(&vg, &vcf, offset);

  // Above inserts new nodes between existing nodes.  So we revise ids
  // to be sorted
  //vg.sort();
  //vg.compact_ids();

  // output modified graph to cout
  vg.serialize_to_ostream(cout);
    
  return 0;
}


