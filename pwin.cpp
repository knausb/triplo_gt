// Compile with:
// g++ -std=c++0x pwin.cpp -lgmpxx -lgmp -o pwin
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
//#include <math.h> /* log10 */
//#include <gmp.h> // GNU multiple precision arithmetic library, not standard: libgmp3-dev.
#include <boost/algorithm/string.hpp> // Not standard on Macs or Ubuntu: libboost1.46-dev.

using namespace std;
using namespace boost;

/* ----- ----- ***** ----- ----- */
/*           Functions           */
/* ----- ----- ***** ----- ----- */



/* ----- ----- ***** ----- ----- */
/*        Print functions        */
/* ----- ----- ***** ----- ----- */

void print_usage(){
  cerr << "  -h print this help message.\n";
  cerr << "  -i infile.\n";
  cerr << "  -w window size [default = 1e3].\n";




  cerr << "  -c print allele counts in genotpye section.\n";
  cerr << "  -e allowable genotyping error [default = 1e-9]; must not be zero.\n";
  cerr << "  -m print vcf header (meta) information.\n";
  cerr << "  -p print phred scaled likelihoods in genotype section.\n";
  cerr << "  -s file with sample names in same order as in\n     the s/bam file, one name per line.\n";
  cerr << "  -t minimum threshold for calling an allele [default = 0].\n";

  cerr << "\n";

}

/* ----- ----- ***** ----- ----- */
/*             Main              */
/* ----- ----- ***** ----- ----- */


int main(int argc, char **argv) {
  string line;
  vector <string> fields;
  vector <string> snames;
//  int opt, header = 0, counts = 0, phred = 0;
//  float error = 0.000000001; // Can not be zero!!!
//  int min_cnt = 0; // Minimum threshold.
  string infile = "NA";
  int opt, win = 999;
  int nsamp;

  /* Parse command line options. */
  while ((opt = getopt(argc, argv, "hi:w:")) != -1) {
    switch (opt) {
      case 'h':
        print_usage();
        exit(EXIT_FAILURE);
      case 'i':
        infile = optarg;
        break;
      case 'w':
        win = atoi(optarg) - 1;
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s [-ce:hmpst:]\n", argv[0]);
        print_usage();
        exit(EXIT_FAILURE);
    }
  }

  if(infile == "NA"){
    cerr << "Error: no input file specified.\n";
    print_usage();
    exit(EXIT_FAILURE);
  }

//  cout << "Got here\n";



  /* Open file and parse line by line. */
  ifstream myfile (infile);
  if (myfile.is_open()){
    getline (myfile,line);

    /* Case no header */
    if(line[0] != '#'){
      split( fields, line, is_any_of( "\t" ) );
      nsamp = fields.size() - 9;  // Determine the number of samples.
//      string snames[nsamp];
//      for(int i=10; i<fields.size(); i++){
//        snames[i-10] = fields[i];
//        snames[i-10].append("Sample_", i-10);
//      }
    }

    /* Omit meta lines. */
    while(line[0] == '#' & line[1] == '#') getline(myfile,line);

    /* Manage header line. */
    if(line[0] == '#' & line[1] == 'C'){
      split( fields, line, is_any_of( "\t" ) );
      nsamp = fields.size() - 9;  // Determine the number of samples.
//      for(int i=10; i<fields.size(); i++){
//        snames[i-10] = fields[i];
//      }
    }

//    for(int i=0; i < nsamp; i++){cout << snames[i] << "\n";}

    cout << line << '\n';


    while ( getline (myfile,line) ){

//      cout << line << '\n';
      split( fields, line, is_any_of( "\t" ) );
      cout << fields[0] << "\t" << fields[1] << "\n";

    }
    myfile.close();
  }

  else cout << "Unable to open file"; 


  /*  */


  return 0;
}

