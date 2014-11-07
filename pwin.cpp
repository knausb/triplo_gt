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
//  cerr << "  -i infile.\n";
  cerr << "  -w window size [default = 1e3].\n";

/*
  cerr << "  -c print allele counts in genotpye section.\n";
  cerr << "  -e allowable genotyping error [default = 1e-9]; must not be zero.\n";
  cerr << "  -m print vcf header (meta) information.\n";
  cerr << "  -p print phred scaled likelihoods in genotype section.\n";
  cerr << "  -s file with sample names in same order as in\n     the s/bam file, one name per line.\n";
  cerr << "  -t minimum threshold for calling an allele [default = 0].\n";
*/

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

void print_null(string chromo, int start, int stop, vector <string> snames){
  cout << chromo << "\t" << start << "\t" << stop << "\t";
  cout << "GT:RD:P1,P2,P3,P4";
  for(int i=0; i<snames.size(); i++){
    cout << "\t" << "0:0:0,0,0,0";
  }
  cout << "\n";
}

/* ----- ----- ***** ----- ----- */
/*             Main              */
/* ----- ----- ***** ----- ----- */


int main(int argc, char **argv) {
  string line;
  vector <string> snames;
  vector <string> fields;
//  vector <string> lines;
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
//        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: %s [-h:w:]\n", argv[0]);
        print_usage();
        exit(EXIT_FAILURE);
    }
  }

  /* Parse line by line. */
  getline (cin,line);

  /* Case no header */
  if(line[0] != '#'){
    cout << "Error, this does not appear to be a properly formatted file:\n";
    cout << line << "\n\n";
    exit(EXIT_FAILURE);
  }

  /* Omit meta lines. */
  while(line[0] == '#' & line[1] == '#') getline(cin,line);

  /* Manage header line. */
  if(line[0] == '#' & line[1] == 'C'){
    split( fields, line, is_any_of( "\t" ) );
    nsamp = fields.size() - 9;  // Determine the number of samples.
    for(int i=9; i<fields.size(); i++){
      snames.push_back(fields[i]);
    }
  }

  /* Process records. */
  print_header(snames);

  int start = 1;
  int stop = start + win;

  /* Read in the first line. */
  getline (cin,line);
  split( fields, line, is_any_of( "\t" ) );


  /* If there are no records in first windows. */
  if(atoi(fields[1].c_str()) > stop){
    int i = atoi(fields[1].c_str())/(win+1);
    for(int j=0; j<i; j++){
      print_null(fields[0], start, stop, snames);
      start = stop + 1;
      stop = start + win;
    }
  }

  /* Windows with records. */
//  while ( getline (cin,line) ){
  std::vector <string> lines;
  lines.push_back(line);
  while(getline (cin,line)){
    split( fields, line, is_any_of( "\t" ) );

    if(atoi(fields[1].c_str()) <= stop){
      lines.push_back(line);
    } else {
      cout << "Window " << start << " to " << stop << "\n";
      cout << "lines contains " << lines.size() << " elements\n";
      cout << lines[0] << "\n";
      cout << lines[lines.size()-1] << "\n";
      cout << "\n";

//      proc_win();
      lines.clear();

      /* Begin the next window. */
      lines.push_back(line);
      // last???
      start = stop + 1;
      stop = start + win;
      split( fields, line, is_any_of( "\t" ) );
      if(atoi(fields[1].c_str()) > stop){
        int i = atoi(fields[1].c_str())/(win+1);
        for(int j=0; j<i; j++){
          print_null(fields[0], start, stop, snames);
          start = stop + 1;
          stop = start + win;
        }
      }
    }
  }
//    vector < vector<int> > rds(nsamp);
//      cout << line << '\n';
//    split( fields, line, is_any_of( "\t" ) );
//      cout << fields[0] << "\t" << fields[1] << "\n";


//  }




  /*  */

  return 0;
}

