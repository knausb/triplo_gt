// Compile with:
// g++-4.9 -std=c++0x pwin2.cpp -o pwin2
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <regex>

#include <boost/algorithm/string.hpp> // Not standard on Macs or Ubuntu: libboost1.46-dev.

using namespace std;



/* ----- ----- ***** ----- ----- */
/*           Functions           */
/* ----- ----- ***** ----- ----- */

void init_win_sums(int nsamp, int nGTs[], int nNAs[],
                   int RDs[], int nPDs[][7], int MPs[]){
  for(int i=0; i<nsamp; i++){
    nGTs[i] = 0;
    nNAs[i] = 0;
    RDs[i]  = 0;
    MPs[i]  = 0;
    for(int j=0; j<7; j++){
      nPDs[i][j] = 0;
    }
  }
}


void proc_sample(int RD, int GT, int &nGTs, int &nNAs, int &RDs, int nPDs[7], string data){
  vector <string> fields;
  boost::algorithm::split( fields, data, boost::algorithm::is_any_of( ":" ) );

  RDs = RDs + stoi(fields[RD]);

}



void proc_win(int nsamp, int nGTs[], int nNAs[], int RDs[],
              int nPDs[][7], int MPs[],  
              std::vector <string> lines){

  vector <string> data;
  vector <string> fields;
//  cout << lines[0] << "\n";

  for(int i=0; i<lines.size(); i++){
    boost::algorithm::split( data, lines[i], boost::algorithm::is_any_of( "\t" ) );
    boost::algorithm::split( fields, data[8], boost::algorithm::is_any_of( ":" ) );
//    cout << data[8] << "\n";

    /* Determine format positions */
    int RD = 0, GT = 0;
//    int CT = 0;
    for(int i=0; i<fields.size(); i++){
      if(fields[i] == "RD"){RD = i;}
//      if(fields[i] == "CT"){CT = i;}
      if(fields[i] == "GT"){GT = i;}
    }


    for(int j=9; j<data.size(); j++){
//      cout << fields[j] << "\t";
      proc_sample(RD, GT, nGTs[j-9], nNAs[j-9], RDs[j-9], nPDs[j-9], data[j]);
    }
//  cout << "\n";
  }
}



/* ----- ----- ***** ----- ----- */
/*        Print functions        */
/* ----- ----- ***** ----- ----- */

void print_usage(){
  cerr << "  -h print this help message.\n";
  cerr << "  -w window size [default = 1e3].\n";

  cerr << "\n";
}

void print_header(vector <string> snames){
  cout << "##fileformat=PVCFv0.0\n";
  cout << "##source=pwinv0.0\n";
  cout << "##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Number of genotypes called\">\n";
  cout << "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Read depth over all called genotypes\">\n";
  cout << "##FORMAT=<ID=PD,Number=7,Type=Integer,Description=\"Number for each class of genotypes (1,2,3,4,5,6,7)\">\n";
  cout << "##FORMAT=<ID=MP,Number=1,Type=Integer,Description=\"Majority rule ploidy\">\n";
  cout << "#CHROM\tSTART\tEND\tVARIANTS\t";
  cout << "POS" << "\t" << "QUAL" << "\t" << "." << "\t" << "." << "\t";
  cout << "FORMAT";
  for(int i=0; i<snames.size(); i++){
    cout << "\t" << snames[i];
  }
  cout << "\n";
}

void print_null(string chromo, int start, int stop, vector <string> snames){
  cout << chromo << "\t" << start << "\t" << stop << "\t";
  cout << 0 << "\t";
  cout << "." << "\t" << "." << "\t" << "." << "\t" << "." << "\t";
//  cout << "GT:RD:P1,P2,P3,P4";
  cout << "GT:RD:PD:MP";
  for(int i=0; i<snames.size(); i++){
    cout << "\t" << "0:0:0,0,0,0,0,0,0";
  }
  cout << "\n";
}


/* ----- ----- ***** ----- ----- */
/*             Main              */
/* ----- ----- ***** ----- ----- */


int main(int argc, char **argv) {
  string line;
  std::vector <string> lines;
  vector <string> snames;
  vector <string> fields;
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


  /* Read in the first line. */
  getline (cin,line);

  /* In case there is no header */
  if(line[0] != '#'){
    cerr << "Error, this does not appear to be a properly formatted file:\n";
    cerr << "Expecting a header beginning with '#'.\n";
    cerr << line << "\n\n";
    exit(EXIT_FAILURE);
  }


  while(line[0] == '#' & line[1] == '#'){
    /* Omit meta lines which begin with '##'. */
    getline (cin,line);
//    cout << line << "\n\n";
  }
//  cout << line << "\n\n";

  if (line[0] == '#' & line[1] == 'C'){
    /* Manage header line which begins with '#CHROM'. */
    boost::algorithm::split( fields, line, boost::algorithm::is_any_of( "\t" ) );
    nsamp = fields.size() - 9;  // Determine the number of samples.
    for(int i=9; i<fields.size(); i++){
      snames.push_back(fields[i]);
    }
  } else {
    cerr << "Poorly formed file\n";
    cerr << "Expecting a header line begining with '#CHROM'\n";
    cerr << line;
    cerr << "\n";
    exit(EXIT_FAILURE);
  }


  /* Process records. */
  int start = 1;
  int stop = start + win;
  int last = start;

  int nGTs[nsamp];
  int nNAs[nsamp];
  int  RDs[nsamp];
  int nPDs[nsamp][7];
  int  MPs[nsamp];

  init_win_sums(nsamp, nGTs, nNAs, RDs, nPDs, MPs);

  print_header(snames);


  while(getline (cin,line) ){
    boost::algorithm::split( fields, line, boost::algorithm::is_any_of( "\t" ) );
    if(atoi(fields[1].c_str()) > stop){
      /* Begin a new window. */
      proc_win(nsamp, nGTs, nNAs, RDs, nPDs, MPs, lines);
//      print_win();
      init_win_sums(nsamp, nGTs, nNAs, RDs, nPDs, MPs);
      start = stop + 1;
      stop = start + win;
      lines.clear();
      lines.push_back(line);
    } else {
      lines.push_back(line);
    }
  }

  return 0;
}

// EOF.
