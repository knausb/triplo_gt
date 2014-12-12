// Compile with:
// g++-4.9 -std=c++0x pwin2.cpp -o pwin2
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
//#include <regex>

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



void maj_rule_ploid(int &MPs, int nPDs[]){

//  if(nPDs[1] > nPDs[MPs]){MPs = 1;}
  MPs = 1; // Don't count homozygotes
  if(nPDs[2] > nPDs[MPs]){MPs = 2;}
  if(nPDs[3] > nPDs[MPs]){MPs = 3;}
  if(nPDs[4] > nPDs[MPs]){MPs = 4;}
  if(nPDs[5] > nPDs[MPs]){MPs = 5;}
  if(nPDs[6] > nPDs[MPs]){MPs = 6;}

/*
  cout << "\n";
  cout << "maj_rule_ploid nPDs[0]: " << nPDs[0] << "," << nPDs[1] << ",";
  cout << nPDs[2] << "," << nPDs[3] << "," << nPDs[4] << ",";
  cout << nPDs[5] << "," << nPDs[6];
  cout << ":" << MPs;
  cout << "\n";
*/

}

void proc_ploid(int nsamp, int MPs[], int nPDs[][7]){
  /* Determine most abundant ploidy. */
  for(int i=0; i<nsamp; i++){
    maj_rule_ploid(MPs[i], nPDs[i]);
  }
}




void proc_sample_len(int RD, int GT, int &nGTs, int &nNAs, int &RDs, int nPDs[7], string data){
  vector <string> fields;
  boost::algorithm::split( fields, data, boost::algorithm::is_any_of( ":" ) );
  string qGT = fields[GT];
  RDs = RDs + stoi(fields[RD]);

  if (qGT == "./." || fields[GT] == ".|."){
    nNAs = nNAs + 1;
  } else if (qGT=="A/A"||qGT=="C/C"||qGT=="G/G"||qGT=="T/T"){
    nGTs = nGTs + 1;
    nPDs[0] = nPDs[0] + 1;
  } else if (qGT=="A/C"||qGT=="A/G"||qGT=="A/T"||qGT=="C/G"||qGT=="C/T"||qGT=="G/T"){
    nGTs = nGTs + 1;
    nPDs[1] = nPDs[1] + 1;
  } else if (qGT.length() == 5){
    // Triploid
    nGTs = nGTs + 1;
    nPDs[2] = nPDs[2] + 1;
  } else if (qGT.length() == 7){
    // Tetraploid
    nGTs = nGTs + 1;
    nPDs[3] = nPDs[3] + 1;
  } else if (qGT.length() == 9){
    // Pentaploid
    nGTs = nGTs + 1;
    nPDs[4] = nPDs[4] + 1;
  } else if (qGT.length() == 11){
    // Hexaploid
    nGTs = nGTs + 1;
    nPDs[5] = nPDs[5] + 1;
  } else if (qGT.length() == 13){
    // Septaploid
    nGTs = nGTs + 1;
    nPDs[6] = nPDs[6] + 1;
  }
}






void proc_win(int nsamp, int nGTs[], int nNAs[], int RDs[],
              int nPDs[][7], std::vector <string> lines){

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
//      proc_sample(RD, GT, nGTs[j-9], nNAs[j-9], RDs[j-9], nPDs[j-9], data[j]);
      proc_sample_len(RD, GT, nGTs[j-9], nNAs[j-9], RDs[j-9], nPDs[j-9], data[j]);
//      cout << "nPDs[j][0]: " << nPDs[j] << "\n";
    }
  }



/*
  cout << "\n";
  cout << "proc_win nPDs[0]: " << nPDs[0][0] << "," << nPDs[0][1] << ",";
  cout << nPDs[0][2] << "," << nPDs[0][3] << "," << nPDs[0][4] << ",";
  cout << nPDs[0][5] << "," << nPDs[0][6];
  cout << ":" << MPs[0];
  cout << "\n";
*/

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



void print_win(vector <string> fields, int start, int stop,
               int nsamp, int nGTs[], int nNAs[], int RDs[],
               int nPDs[][7], int MPs[]){
  cout <<  fields[0];
  cout << "\t";
  cout <<  start;
  cout << "\t";
  cout <<  stop;
  cout << "\t" << "." << "\t";
  cout << "." << "\t" << "." << "\t" << "." << "\t" << "." << "\t";
  cout << "GT:RD:PD:MP";
  cout << "\t";
  for(int i=0; i<nsamp; i++){
    cout << nGTs[i];
    cout << ":";
    cout << RDs[i];
    cout << ":";

    cout << nPDs[i][0];
    for(int j=1; j<7; j++){
      cout << "," << nPDs[i][j];
    }

    cout << ":";
    cout << MPs[i] + 1; // Convert zero based ploidy to one based.

    cout << "\t";
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

//  cout << "Initial nPDs[0]: " << nPDs[0][0] << "\n";

  print_header(snames);


  while(getline (cin,line) ){
    boost::algorithm::split( fields, line, boost::algorithm::is_any_of( "\t" ) );
    if(atoi(fields[1].c_str()) > stop){
      /* Process and print the last window. */
      proc_win(nsamp, nGTs, nNAs, RDs, nPDs, lines);
      proc_ploid(nsamp, MPs, nPDs);

/*
      cout << "Post proc_win nPDs[0]: " << nPDs[0][0] << "," << nPDs[0][1] << ",";
      cout << nPDs[0][2] << "," << nPDs[0][3] << "," << nPDs[0][4] << ",";
      cout << nPDs[0][5] << "," << nPDs[0][6];
      cout << ":" << MPs[0];

      cout << "\n";
*/

      print_win(fields, start, stop, nsamp, nGTs, nNAs, RDs, nPDs, MPs);

      /* Begin a new window. */
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
