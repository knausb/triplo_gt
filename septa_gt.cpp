// Compile with:
// g++ -std=c++0x septa_gt.cpp -lgmpxx -lgmp -o septa_gt
#include <iostream>
#include <fstream>
#include <string>
#include <array>
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


/* Calculate possibe ways of getting counts */
double possible_counts(int nuc_cnt[4]){
  int rd = nuc_cnt[0] + nuc_cnt[1] + nuc_cnt[2] + nuc_cnt[3];

  /* Multiple precision variables. */
  mpz_t fac [5]; // Factorials n, A, C, G, T.
  mpf_t pos [3]; // Floats for num, den, pos;

  // Initialize mp ints.
  for(int j=0; j<5; j++){
    mpz_init(fac[j]);
    mpz_set_ui(fac[j],0);
  }

  // Initialize mp floats.
  for(int j=0; j<3; j++){
    mpf_init(pos[j]);
    mpf_set_ui(pos[j],0); 
  }

  // Factorials.
  mpz_fac_ui(fac[0], rd);
  mpz_fac_ui(fac[1], nuc_cnt[0]);
  mpz_fac_ui(fac[2], nuc_cnt[1]);
  mpz_fac_ui(fac[3], nuc_cnt[2]);
  mpz_fac_ui(fac[4], nuc_cnt[3]);

  // Multiply denominators.
  mpz_mul(fac[1],fac[1],fac[2]);
  mpz_mul(fac[1],fac[1],fac[3]);
  mpz_mul(fac[1],fac[1],fac[4]);

  // Recast mp ints to floats for division.
  mpf_set_z(pos[1], fac[0]);
  mpf_set_z(pos[2], fac[1]);

  // Divide.
  mpf_div(pos[0], pos[1], pos[2]);
//    cout << "pos0=" << pos[0] << "\n";
  double posd = mpf_get_d(pos[0]);

  // Clean up the mpz_t handles or else we will leak memory
  for (int j=0; j<5; j++){mpz_clear(fac[j]);}
  for (int j=0; j<3; j++){mpf_clear(pos[j]);}

  return(posd);
}

void int_to_nucs(vector <string> gts[26], vector<pair<int,int>> moves){



}


/* Create type and function to help sort nucleotides. */
typedef std::pair<int,int> mypair;
bool comparator ( const mypair& l, const mypair& r)
    { return l.first > r.first; }


//void counts_2_plh(int mlhs[51], int nuc_cnts[8], float error, int min_cnt, int debug=0){
void counts_2_plh(int mlhs[26], int nuc_cnts[8], float error, int debug=0){

  /* Add forward and reverse counts. */
  int nuc_cnt[4];
  nuc_cnt[0] = nuc_cnts[0] + nuc_cnts[1];
  nuc_cnt[1] = nuc_cnts[2] + nuc_cnts[3];
  nuc_cnt[2] = nuc_cnts[4] + nuc_cnts[5];
  nuc_cnt[3] = nuc_cnts[6] + nuc_cnts[7];


//  cout << "counts_2_plh\t";
//  cout << nuc_cnt[0] << "," << nuc_cnt[1] << "," << nuc_cnt[2] << "," << nuc_cnt[3] << ":"; 


  /* Sort nucleotide counts. */
  vector<pair<int,int>> moves = {
    {nuc_cnt[0], 0},
    {nuc_cnt[1], 1},
    {nuc_cnt[2], 2},
    {nuc_cnt[3], 3}
  };

  std::sort(moves.begin(), moves.end(), comparator);

//  cout << moves[0].first << "," << moves[1].first;
//  cout << "," << moves[2].first << "," << moves[3].first;
//  cout << ":" << moves[0].second << "," << moves[1].second;
//  cout << "," << moves[2].second << "," << moves[3].second;
//  cout << "\n";

  /* Possible ways to obtain observed counts */
  double poss_counts = possible_counts(nuc_cnt);
//  cout << " poss_counts=" << poss_counts;
//  cout << "\n";


  /* --*-- --*-- --*-- */
  /*       Models.     */
  /* --*-- --*-- --*-- */

  /* Likelihoods. */
  float mls [26];
  for(int j=0; j<26; j++){mls[j]=0;}

  /* Homozygote. */
  mls[0] = poss_counts * 
             pow(1-(3*error)/4, moves[0].first) * 
             pow(error/4, moves[1].first+moves[2].first+moves[3].first);

  /* Biallelic heterozygote. */
  mls[1] = poss_counts *
             pow(0.5-error/4, moves[0].first+moves[1].first) *
             pow(error/4, moves[2].first+moves[3].first);

  /* Biallelic triploid. */
  float prop3 = 0.3333333;
  float prop6 = 0.6666667;
  mls[2] = poss_counts * 
              pow(prop6-error/4, moves[0].first) * 
              pow(prop3-error/4, moves[1].first) *
              pow(error/4, moves[2].first+moves[3].first);

  /* Triallelic triploid loci. */
  mls[3] = poss_counts *
             pow(prop3-error/4, moves[0].first+moves[1].first+moves[2].first) *
             pow(error/4, moves[3].first);

  /* Biallelic tetraploid loci */
  mls[4] = poss_counts *
             pow(0.75-error/4, moves[0].first) *
             pow(0.25-error/4, moves[1].first) *
             pow(error/4, moves[2].first+moves[3].first);

  /* Triallelic tetraploid loci */
  mls[5] = poss_counts *
             pow(0.5-error/4, moves[0].first) *
             pow(0.25-error/4, moves[1].first) *
             pow(0.25-error/4, moves[2].first) *
             pow(error/4, moves[3].first);

  /* Tetra-allelic tetraploid loci */
  mls[6] = poss_counts *
             pow(0.25-error/4, moves[0].first) *
             pow(0.25-error/4, moves[1].first) *
             pow(0.25-error/4, moves[2].first) *
             pow(0.25-error/4, moves[3].first);

  /* Bi-allelic pentaploid */
  mls[7] = poss_counts *
             pow(0.8-error/4, moves[0].first) *
             pow(0.2-error/4, moves[1].first) *
             pow(error/4, moves[2].first+moves[3].first);
  mls[8] = poss_counts *
             pow(0.6-error/4, moves[0].first) *
             pow(0.4-error/4, moves[1].first) *
             pow(error/4, moves[2].first+moves[3].first);

  /* Tri-allelic pentaploid */
  mls[9] = poss_counts *
             pow(0.6-error/4, moves[0].first) *
             pow(0.2-error/4, moves[1].first) *
             pow(0.2-error/4, moves[1].first) *
             pow(error/4, moves[3].first);
  mls[10] = poss_counts *
             pow(0.4-error/4, moves[0].first) *
             pow(0.4-error/4, moves[1].first) *
             pow(0.2-error/4, moves[1].first) *
             pow(error/4, moves[3].first);

  /* Tetra-allelic pentaploid */
  mls[11] = poss_counts *
              pow(0.4-error/4, moves[0].first) *
              pow(0.2-error/4, moves[1].first) *
              pow(0.2-error/4, moves[2].first) *
              pow(0.2-error/4, moves[3].first);

  /* Bi-allelic hexaploid */
  mls[12] = poss_counts *
              pow(0.8333333-error/4, moves[0].first) *
              pow(0.1666667-error/4, moves[1].first) *
              pow(error/4, moves[2].first+moves[3].first);

  /* Tri-allelic hexaploid */
  mls[13] = poss_counts *
              pow(0.6666667-error/4, moves[0].first) *
              pow(0.1666667-error/4, moves[1].first) *
              pow(0.1666667-error/4, moves[1].first) *
              pow(error/4, moves[3].first);
  mls[14] = poss_counts *
              pow(0.5-error/4, moves[0].first) *
              pow(0.3333333-error/4, moves[1].first) *
              pow(0.1666667-error/4, moves[1].first) *
              pow(error/4, moves[3].first);

  /* Tetra-allelic hexaploid */
  mls[15] = poss_counts *
              pow(0.5-error/4, moves[0].first) *
              pow(0.1666667-error/4, moves[1].first) *
              pow(0.1666667-error/4, moves[2].first) *
              pow(0.1666667-error/4, moves[3].first);
  mls[16] = poss_counts *
              pow(0.3333333-error/4, moves[0].first) *
              pow(0.3333333-error/4, moves[1].first) *
              pow(0.1666667-error/4, moves[2].first) *
              pow(0.1666667-error/4, moves[3].first);

  /* Bi-allelic septaploid */
  mls[17] = poss_counts *
              pow(0.8571429-error/4, moves[0].first) *
              pow(0.1428571-error/4, moves[1].first) *
              pow(error/4, moves[2].first+moves[3].first);
  mls[18] = poss_counts *
              pow(0.7142857-error/4, moves[0].first) *
              pow(0.2857143-error/4, moves[1].first) *
              pow(error/4, moves[2].first+moves[3].first);
  mls[19] = poss_counts *
              pow(0.5714286-error/4, moves[0].first) *
              pow(0.4285714-error/4, moves[1].first) *
              pow(error/4, moves[2].first+moves[3].first);

  /* Tri-allelic septaploid */
  mls[20] = poss_counts *
              pow(0.7142857-error/4, moves[0].first) *
              pow(0.1428571-error/4, moves[1].first) *
              pow(0.1428571-error/4, moves[1].first) *
              pow(error/4, moves[3].first);
  mls[21] = poss_counts *
              pow(0.5714286-error/4, moves[0].first) *
              pow(0.2857143-error/4, moves[1].first) *
              pow(0.1428571-error/4, moves[1].first) *
              pow(error/4, moves[3].first);
  mls[22] = poss_counts *
              pow(0.3333333-error/4, moves[0].first) *
              pow(0.3333333-error/4, moves[1].first) *
              pow(0.1428571-error/4, moves[1].first) *
              pow(error/4, moves[3].first);

  /* Tetra-allelic septaploid */
  mls[23] = poss_counts *
              pow(0.5714286-error/4, moves[0].first) *
              pow(0.1428571-error/4, moves[1].first) *
              pow(0.1428571-error/4, moves[2].first) *
              pow(0.1428571-error/4, moves[3].first);
  mls[24] = poss_counts *
              pow(0.4285714-error/4, moves[0].first) *
              pow(0.2857143-error/4, moves[1].first) *
              pow(0.1428571-error/4, moves[2].first) *
              pow(0.1428571-error/4, moves[3].first);
  mls[25] = poss_counts *
              pow(0.2857143-error/4, moves[0].first) *
              pow(0.2857143-error/4, moves[1].first) *
              pow(0.2857143-error/4, moves[2].first) *
              pow(0.1428571-error/4, moves[3].first);


  /* Phred scale the likelihoods. */
//  for (int j = 0; j < 51; j++){
  for (int j = 0; j < 26; j++){
    if (mls[j] == 0){
      mls[j] = 9999;
    } else {
      mls[j] = trunc(-10 * log10(mls[j]));
      if(mls[j] == -0){mls[j] = 0;}
    }
  }

  /* Convert int to nucleotides */
  vector <string> gts[26];
//  int_to_nucs(gts, moves);

//  debug = 1;
  if(debug == 1){
    cout << "\n";
    cout << "*** Debug counts_2_phredlh ***\n";
    cout << nuc_cnt[0] << "," << nuc_cnt[1] << "," << nuc_cnt[2]<< "," << nuc_cnt[3] << ":";
    cout << "poss_counts=" << poss_counts;
    cout << "\n";
  }

  /* Transfer likelihoods to parent array. */
//  for (int j = 0; j < 51; j++){
  for (int j = 0; j < 26; j++){
//    pls[i][j] = int(mls[j]);
    mlhs[j] = int(mls[j]);
  }

}


/* Likelihoods */
//void mult_pl(int pls[][26], int nsamp, float err, int nuc_cnts[][8], int min_cnt){

/* Parse each site to samples */
void mult_pl(int pls[][26], int nsamp, float err, int nuc_cnts[][8]){
//  cout <<  "##### New site (row) #####\n";
  for(int i=0; i<nsamp; i++){
//    counts_2_plh(pls[i], nuc_cnts[i], err, min_cnt);
    counts_2_plh(pls[i], nuc_cnts[i], err);
  }
}



/* ----- ----- ***** ----- ----- */
/*        Print functions        */
/* ----- ----- ***** ----- ----- */

void print_usage(){
  cerr << "\n";
  cerr << "  -c print allele counts in genotpye section.\n";
  cerr << "  -e allowable genotyping error [default = 1e-9]; must not be zero.\n";
  cerr << "  -h print this help message.\n";
  cerr << "  -m print vcf header (meta) information.\n";
  cerr << "  -p print phred scaled likelihoods in genotype section.\n";
  cerr << "  -s file with sample names in same order as in\n     the s/bam file, one name per line.\n";
//  cerr << "  -t minimum threshold for calling an allele [default = 0].\n";
  cerr << "\n";
}


//void print_header(float error, int min_cnt, string sfile){
void print_header(float error, string sfile){
  cout << "##fileformat=VCFv4.2\n";
  cout << "##source=gtV0.0.0\n";
//  cout << "##FILTER=<ID=NA,Description=\"Per nucleotide minimum threshold of " << min_cnt << "\">\n";
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
//  int min_cnt = 0; // Minimum threshold.
  string sfile = "NA";

  /* Parse command line options. */
//  while ((opt = getopt(argc, argv, "ce:hmps:t:")) != -1) {
  while ((opt = getopt(argc, argv, "ce:hmps:")) != -1) {
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
      case 's': // file containing sample names
        sfile = optarg;
        break;
//      case 't': // Minimum threshold
//        min_cnt = atoi(optarg);
//        break;
      default: /* '?' */
//        fprintf(stderr, "Usage: %s [-ce:hmpst:]\n", argv[0]);
        fprintf(stderr, "Usage: %s [-ce:hmps]\n", argv[0]);
        print_usage();
        exit(EXIT_FAILURE);
    }
  }

  /* Print header */
//  if(header == 1){print_header(error, min_cnt, sfile);}
  if(header == 1){print_header(error, sfile);}

  /* Parse line by line or site by site. */
  while (getline(cin,lineInput)) {
//    cout << "##### ----- New variant (row) ----- #####\n";
//    cout << lineInput << "\n";
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
//        cout << "Processing sample " << sampn << "\n";
        for(int j=0; j<8; j++){nuc_cnts[sampn][j] = 0;}
        if(fields[i] == "*"){
//          cout << fields[i] << ": no data string\n";
          // No data.
          //for(int j=0; j<8; j++){nuc_cnts[sampn][j] = 0;}
        } else {
          // Count each nucleotide.
//          cout << fields[i] << ": count string\n";
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
    mult_pl(pls, nsamp, error, nuc_cnts);

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
