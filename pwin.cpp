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

void print_header(vector <string> snames){
  cout << "##fileformat=PVCFv0.0\n";
  cout << "##source=pwinv0.0\n";
  cout << "#CHROM\tSTART\tEND\tFORMAT";
  for(int i=0; i<snames.size(); i++){
    cout << "\t" << snames[i];
  }
  cout << "\n";
}


/* ----- ----- ***** ----- ----- */
/*             Main              */
/* ----- ----- ***** ----- ----- */


int main(int argc, char **argv) {
  string line;
  vector <string> fields;
  vector <string> snames;
  int opt, win = 999;
  int nsamp;

  /* Parse command line options. */
  while ((opt = getopt(argc, argv, "h:w:")) != -1) {
    switch (opt) {
      case 'h':
        print_usage();
        exit(EXIT_FAILURE);
      case 'w':
        win = atoi(optarg) - 1;
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s [-h:w:]\n", argv[0]);
        print_usage();
        exit(EXIT_FAILURE);
    }
  }

  /* Parse line by line. */
    getline (cin,line);

    /* Case no header */
    if(line[0] != '#'){
      split( fields, line, is_any_of( "\t" ) );
      nsamp = fields.size() - 9;  // Determine the number of samples.
      for(int i=9; i<fields.size(); i++){
        cout << fields[i] << "\n";
        snames.push_back(fields[i]);
      }
    }

    /* Omit meta lines. */
    while(line[0] == '#' & line[1] == '#') getline(cin,line);

    /* Manage header line. */
    if(line[0] == '#' & line[1] == 'C'){
      split( fields, line, is_any_of( "\t" ) );
      nsamp = fields.size() - 9;  // Determine the number of samples.
      for(int i=9; i<fields.size(); i++){
//        cout << fields[i] << "\n";
        snames.push_back(fields[i]);
      }
    }

//    cout << "\n\n";
//    for(int i=0; i<snames.size(); i++){cout << snames[i] << "\n";}



//    cout << line << '\n';

    /* 
       Probably have a line reading issue here.
       If there was no header we're on a data line.
       If there was a header, we're on that.
       Perhaps just require that there is a meta and header line.

   */

    print_header(snames);

    int start = 1;
    int stop = start + win;
    while ( getline (cin,line) ){
      vector < vector<int> > rds(nsamp);


//      cout << line << '\n';
      split( fields, line, is_any_of( "\t" ) );
//      cout << fields[0] << "\t" << fields[1] << "\n";




    }

  /*  */

  return 0;
}

