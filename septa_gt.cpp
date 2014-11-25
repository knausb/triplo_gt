// Compile with:
// g++ -std=c++0x septa_gt.cpp -lgmpxx -lgmp -o septa_gt
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <math.h> /* log10 */
#include <gmp.h> // GNU multiple precision arithmetic library, not standard: libgmp3-dev.
#include <boost/algorithm/string.hpp> // Not standard on Macs or Ubuntu: libboost1.46-dev.

using namespace std;
//using namespace boost;

/* ----- ----- ***** ----- ----- */
/*           Functions           */
/* ----- ----- ***** ----- ----- */

void refA(std::string nucs, int sampn, int nuc_cnt[][8]){
  for (int i = 0; i < nucs.size(); i++){
    if (nucs[i] == '.') nuc_cnt[sampn][0]++;
    if (nucs[i] == ',') nuc_cnt[sampn][1]++;
    if (nucs[i] == 'C') nuc_cnt[sampn][2]++;
    if (nucs[i] == 'c') nuc_cnt[sampn][3]++;
    if (nucs[i] == 'G') nuc_cnt[sampn][4]++;
    if (nucs[i] == 'g') nuc_cnt[sampn][5]++;
    if (nucs[i] == 'T') nuc_cnt[sampn][6]++;
    if (nucs[i] == 't') nuc_cnt[sampn][7]++;
  }
}

void refC(std::string nucs, int sampn, int nuc_cnt[][8]){
  for (int i = 0; i < nucs.size(); i++){
    if (nucs[i] == 'A') nuc_cnt[sampn][0]++;
    if (nucs[i] == 'a') nuc_cnt[sampn][1]++;
    if (nucs[i] == '.') nuc_cnt[sampn][2]++;
    if (nucs[i] == ',') nuc_cnt[sampn][3]++;
    if (nucs[i] == 'G') nuc_cnt[sampn][4]++;
    if (nucs[i] == 'g') nuc_cnt[sampn][5]++;
    if (nucs[i] == 'T') nuc_cnt[sampn][6]++;
    if (nucs[i] == 't') nuc_cnt[sampn][7]++;
  }
}

void refG(std::string nucs, int sampn, int nuc_cnt[][8]){
  for (int i = 0; i < nucs.size(); i++){
    if (nucs[i] == 'A') nuc_cnt[sampn][0]++;
    if (nucs[i] == 'a') nuc_cnt[sampn][1]++;
    if (nucs[i] == 'C') nuc_cnt[sampn][2]++;
    if (nucs[i] == 'c') nuc_cnt[sampn][3]++;
    if (nucs[i] == '.') nuc_cnt[sampn][4]++;
    if (nucs[i] == ',') nuc_cnt[sampn][5]++;
    if (nucs[i] == 'T') nuc_cnt[sampn][6]++;
    if (nucs[i] == 't') nuc_cnt[sampn][7]++;
  }
}

void refT(std::string nucs, int sampn, int nuc_cnt[][8]){
  for (int i = 0; i < nucs.size(); i++){
    if (nucs[i] == 'A') nuc_cnt[sampn][0]++;
    if (nucs[i] == 'a') nuc_cnt[sampn][1]++;
    if (nucs[i] == 'C') nuc_cnt[sampn][2]++;
    if (nucs[i] == 'c') nuc_cnt[sampn][3]++;
    if (nucs[i] == 'G') nuc_cnt[sampn][4]++;
    if (nucs[i] == 'g') nuc_cnt[sampn][5]++;
    if (nucs[i] == '.') nuc_cnt[sampn][6]++;
    if (nucs[i] == ',') nuc_cnt[sampn][7]++;
  }
}

void get_rd(int rds[], int sampn, int nuc_cnts[][8]){
//  cout << "\nget_rd sampn: " << sampn << "\n";
  for(int i=0; i<sampn; i++){ // Samples
    rds[i] = nuc_cnts[i][0];
    for(int j=1; j<8; j++){
      rds[i] = rds[i] + nuc_cnts[i][j];
    }
  }
}





/* ----- ----- ***** ----- ----- */
/*        Print functions        */
/* ----- ----- ***** ----- ----- */

void print_usage(){
  cerr << "  -c print allele counts in genotpye section.\n";
  cerr << "  -e allowable genotyping error [default = 1e-9]; must not be zero.\n";
  cerr << "  -h print this help message.\n";
  cerr << "  -m print vcf header (meta) information.\n";
  cerr << "  -p print phred scaled likelihoods in genotype section.\n";
  cerr << "  -s file with sample names in same order as in\n     the s/bam file, one name per line.\n";
  cerr << "  -t minimum threshold for calling an allele [default = 0].\n";
  cerr << "\n";
}


void print_header(float error, int min_cnt, string sfile){
  cout << "##fileformat=VCFv4.2\n";
  cout << "##source=gtV0.0.0\n";
  cout << "##FILTER=<ID=NA,Description=\"Per nucleotide minimum threshold of " << min_cnt << "\">\n";
  cout << "##FILTER=<ID=NA,Description=\"Error rate of " << error << " used for likelihood calculation\">\n";
  cout << "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Read depth\">\n";
  cout << "##FORMAT=<ID=CT,Number=4,Type=Integer,Description=\"Count of each nucleotide (A,C,G,T)\">\n";
  cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
  cout << "##FORMAT=<ID=PL,Number=26,Type=Integer,Description=\"Phred scaled likelihood for 26 genotpyes\">\n";
  cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
//
  if(sfile != "NA"){
    string line;
    ifstream myfile (sfile);
    if (myfile.is_open())
    {
      while ( getline (myfile,line) )
      {
      cout << "\t" << line;
      }
    myfile.close();
    }
    else cout << "Unable to open file"; 
  }
//
  cout << "\n";
}




/* ----- ----- ***** ----- ----- */
/*             Main              */
/* ----- ----- ***** ----- ----- */


int main(int argc, char **argv) {
  string lineInput;
  vector <string> fields;
  int opt; // options
  int header = 0; // print header/meta data
  int counts = 0; // print counts 
  int phred = 0; // Scale likelihoods as phred values
  float error = 0.000000001; // Can not be zero!!!
  int min_cnt = 0; // Minimum threshold.
  string sfile = "NA";

  /* Parse command line options. */
  while ((opt = getopt(argc, argv, "ce:hmps:t:")) != -1) {
    switch (opt) {
      case 'c': // print counts
        counts = 1;
        break;
      case 'e': // permisable error
        error = atof(optarg);
        break;
      case 'h': // print usage and exit
        print_usage();
        exit(EXIT_FAILURE);
      case 'm': // print header/meta data
        header = 1;
        break;
      case 'p': // phred scale likelihoods
        phred = 1;
        break;
      case 's': // Figure this out!!!
        sfile = optarg;
        break;
      case 't': // Minimum threshold
        min_cnt = atoi(optarg);
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s [-ce:hmpst:]\n", argv[0]);
        print_usage();
        exit(EXIT_FAILURE);
    }
  }

  /* Print header */
  if(header == 1){print_header(error, min_cnt, sfile);}


  /* Parse line by line or site by site. */
  while (getline(cin,lineInput)) {
    split( fields, lineInput, boost::algorithm::is_any_of( "\t " ) );
    /* Declare variables */
    int nsamp = (fields.size()-3)/3;  // Determine the number of samples.
//    int nsamp = (fields.size()-2)/2;  // Determine the number of samples.
    int nuc_cnts [nsamp][8]; // A,a,C,c,G,g,T,t.
    int rds [nsamp]; // Read depth.
    string gts [nsamp]; // Genotypes.
//    int pls [nsamp][51];  // Phred scaled likelihoods. 
    int pls [nsamp][26];  // Phred scaled likelihoods. 

    /* Initialize variables. */
    for(int i=0; i<nsamp; i++){
      rds[i] = 0;
      gts[i] = "./.";
    }

    /* Process each sample in the line. */
    int sampn = -1;  // sample counter.
    for(int i = 4; i <= fields.size(); i++){
      if((i-1) % 3 == 0){
        /* New sample */
        sampn++;
        for(int j=0; j<8; j++){nuc_cnts[sampn][j] = 0;}
        if(fields[i] == "*"){
          // No data.
          //for(int j=0; j<8; j++){nuc_cnts[sampn][j] = 0;}
        } else {
          // Count each nucleotide.
          if(fields[2] == "A"){refA(fields[i], sampn, nuc_cnts);}
          if(fields[2] == "C"){refC(fields[i], sampn, nuc_cnts);}
          if(fields[2] == "G"){refG(fields[i], sampn, nuc_cnts);}
          if(fields[2] == "T"){refT(fields[i], sampn, nuc_cnts);}
        }
      }
    }

    /* Get read depths */
    get_rd(rds, nsamp, nuc_cnts);

    /* Calculate Phred-scaled likelihoods */
//    mult_pl(pls, nsamp, error, nuc_cnts, min_cnt);

    /* Determine a genotype */
//    det_gt(gts, nsamp, rds, nuc_cnts, pls);

    /* Minimum threshold. */
//    min_tresh(gts, nsamp, min_cnt, nuc_cnts);

    /* Print locus. */
//    int unique_gt = cnt_gts(nsamp, gts);
//    if(unique_gt > 1){
//      print_locus(fields, counts, phred, nsamp, rds, nuc_cnts, pls, gts);
      /* Debug */
//      debug1(fields, counts, phred, nsamp, rds, nuc_cnts, pls, gts, error, min_cnt);
//    }
  }
  return 0;
}



/* ----- ----- ***** ----- ----- */
// EOF.
