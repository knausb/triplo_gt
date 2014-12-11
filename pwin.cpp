// Compile with:
// g++ -std=c++0x pwin.cpp -lgmpxx -lgmp -o pwin
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <regex>

//#include <math.h> /* log10 */
//#include <gmp.h> // GNU multiple precision arithmetic library, not standard: libgmp3-dev.
#include <boost/algorithm/string.hpp> // Not standard on Macs or Ubuntu: libboost1.46-dev.
//#include <boost/regex/v4/regex.hpp>

using namespace std;
using namespace boost;

/* ----- ----- ***** ----- ----- */
/*           Functions           */
/* ----- ----- ***** ----- ----- */

void count_gts_regex(int nGT[], int P1[], int P2[], int P3[],
                     int P4[], int P5[], int P6[], int P7[],
                     int sampn, string GT){

  cout << "Got here.\n";
  cout << "GT: " << GT << "\n";

  /* Homozygotes. */
  std::regex query("(A/A|C/C|G/G|T/T)");
  if (std::regex_match (GT, query) ){
    cout << "Found a homozygote.\n";

  }

//  if (std::regex_match (GT, std::regex("(^A/A$|^C/C$|^G/G$|^T/T$)") )){
//  if (std::regex_match (GT, std::regex("A/A|C/C|G/G|T/T", std::regex_constants::ECMAScript) )){
//  boost::basic_regex e("A/A|C/C|G/G|T/T");
//  boost::regex e("[0-9]");
//  boost::regex e("(\\d{4}[- ]){3}\\d{4}");
//  boost::regex e("A");

//  if(regex_match(GT, e)){
/*    cout << "Found a homozygote.\n";
    nGT[sampn]++;
    P1[sampn]++;
  }
*/

  /* Heterozygotes. */
/*
  if (std::regex_match (GT, std::regex("(^A/C$|^A/G$|^A/T$|^C/G$|^C/T$|^G/T$)") )){
    nGT[sampn]++;
    P2[sampn]++;
  }*/

  /* Triploids. */
/*
  if (std::regex_match (GT, std::regex ("[ACGT]/[ACGT]/[ACGT]", std::regex_constants::basic) )){
//  if (std::regex_match (GT, std::regex("^[ACGT]/[ACGT]/[ACGT]$") )){
    cout << "Found a triploid.\n";
    nGT[sampn]++;
    P3[sampn]++;
  }*/

  /* Tetraploids. */
/*
  if (std::regex_match (GT, std::regex("^[ACGT]/[ACGT]/[ACGT]/[ACGT]$", std::regex_constants::basic) )){
    cout << "Found a tetraploid.\n";
    nGT[sampn]++;
    P4[sampn]++;
  }*/

  /* Pentaploids. */
/*  if (std::regex_match (GT, std::regex("^[ACGT]/[ACGT]/[ACGT]/[ACGT]/[ACGT]$", std::regex_constants::basic) )){
    cout << "Found a pentaploid.\n";
    nGT[sampn]++;
    P5[sampn]++;
  }*/

  /* Hexaploids. */
/*  if (std::regex_match (GT, std::regex("^[ACGT]/[ACGT]/[ACGT]/[ACGT]/[ACGT]/[ACGT]$", std::regex_constants::basic) )){
    nGT[sampn]++;
    P6[sampn]++;
  }*/

  /* Septaploids. */
/*  if (std::regex_match (GT, std::regex("^[ACGT]/[ACGT]/[ACGT]/[ACGT]/[ACGT]/[ACGT]/[ACGT]$", std::regex_constants::basic) )){
    nGT[sampn]++;
    P7[sampn]++;
  }*/
}


void count_gts(int nGT[], int P1[], int P2[], int P3[], int P4[], int sampn, string GT){

  /* Homozygotes. */
  if(GT == "A/A"){
    nGT[sampn]++;
    P1[sampn]++;
  }
  if(GT == "C/C"){
    nGT[sampn]++;
    P1[sampn]++;
  }
  if(GT == "G/G"){
    nGT[sampn]++;
    P1[sampn]++;
  }
  if(GT == "T/T"){
    nGT[sampn]++;
    P1[sampn]++;
  }

  /* Heterozygotes. */
  if(GT == "A/C"){
    nGT[sampn]++;
    P2[sampn]++;
  }
  if(GT == "A/G"){
    nGT[sampn]++;
    P2[sampn]++;
  }
  if(GT == "A/T"){
    nGT[sampn]++;
    P2[sampn]++;
  }
  if(GT == "C/G"){
    nGT[sampn]++;
    P2[sampn]++;
  }
  if(GT == "C/T"){
    nGT[sampn]++;
    P2[sampn]++;
  }
  if(GT == "G/T"){
    nGT[sampn]++;
    P2[sampn]++;
  }

  /* Bi-allelic triploids. */
  if(GT == "A/A/C"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "A/A/G"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "A/A/T"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "C/C/A"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "C/C/G"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "C/C/T"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "G/G/A"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "G/G/C"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "G/G/T"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "T/T/A"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "T/T/C"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "T/T/G"){
    nGT[sampn]++;
    P3[sampn]++;
  }

  /* Tri-allelic triploids. */
  if(GT == "A/C/G"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "A/C/T"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "A/G/T"){
    nGT[sampn]++;
    P3[sampn]++;
  }
  if(GT == "C/G/T"){
    nGT[sampn]++;
    P3[sampn]++;
  }

  /* Bi-allelic tetraploids. */
  if(GT == "A/A/A/C"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "A/A/A/G"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "A/A/A/T"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "C/C/C/A"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "C/C/C/G"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "C/C/C/T"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "G/G/G/A"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "G/G/G/C"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "G/G/G/T"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "T/T/T/A"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "T/T/T/C"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "T/T/T/G"){
    nGT[sampn]++;
    P4[sampn]++;
  }

  /* Tri-allelic tetraploids. */
  if(GT == "A/A/C/G"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "A/A/C/T"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "A/A/G/T"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "C/C/A/G"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "C/C/A/T"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "C/C/G/T"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "G/G/A/C"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "G/G/A/T"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "G/G/C/T"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "T/T/A/C"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "T/T/A/G"){
    nGT[sampn]++;
    P4[sampn]++;
  }
  if(GT == "T/T/C/G"){
    nGT[sampn]++;
    P4[sampn]++;
  }

  /* Tetra-allelic tetraploids. */
  if(GT == "A/C/G/T"){
    nGT[sampn]++;
    P4[sampn]++;
  }

}




void proc_win(string chromo, int start, int stop, vector <string> lines){

  /* Initialize on first row (variant). */

  vector <string> fields;
  split( fields, lines[0], is_any_of( "\t" ) );
  vector <string> format;
  vector <string> sample;

//  cout << fields[8] << "\n";

  int nsamps = fields.size() - 9;

  int nGTs[nsamps];
  int RDs[nsamps];
  int P1[nsamps];
  int P2[nsamps];
  int P3[nsamps];
  int P4[nsamps];
  int P5[nsamps];
  int P6[nsamps];
  int P7[nsamps];

  /* Initialize arrays to zero. */
  for(int i=0; i<nsamps; i++){
    nGTs[i] = 0;
    RDs[i] = 0;
    P1[i] = 0;
    P2[i] = 0;
    P3[i] = 0;
    P4[i] = 0;
    P5[i] = 0;
    P6[i] = 0;
    P7[i] = 0;
  }

  /* Determine format positions */
  int RD = 0, CT = 0, GT = 0;
  split( format, fields[8], is_any_of( ":" ) );
  for(int i=0; i<format.size(); i++){
    if(format[i] == "RD"){RD = i;}
    if(format[i] == "CT"){CT = i;}
    if(format[i] == "GT"){GT = i;}
  }

  /* Iterate over samples. */
  for(int i=9; i<fields.size(); i++){
    split(sample, fields[i], is_any_of(":"));
    RDs[i-9] = RDs[i-9] + atoi(sample[RD].c_str());
//    count_gts(nGTs, P1, P2, P3, P4, i-9, sample[GT]);
    count_gts_regex(nGTs, P1, P2, P3, P4, P5, P6, P7, i-9, sample[GT]);
  }


  /* Iterate over rows (variants). */
  for(int i=1; i<lines.size(); i++){
    split( fields, lines[i], is_any_of( "\t" ) );
    /* Determine format positions */
    split( format, fields[8], is_any_of( ":" ) );
    for(int j=0; j<format.size(); j++){
      if(format[j] == "RD"){RD = j;}
      if(format[j] == "CT"){CT = j;}
      if(format[j] == "GT"){GT = j;}
    }

    /* Iterate over samples. */
    for(int j=9; j<fields.size(); j++){
      split(sample, fields[j], is_any_of(":"));
      RDs[j-9] = RDs[j-9] + atoi(sample[RD].c_str());
//      count_gts(nGTs, P1, P2, P3, P4, j-9, sample[GT]);
      count_gts_regex(nGTs, P1, P2, P3, P4, P5, P6, P7, i-9, sample[GT]);
    }

  }


  /* Print. */
  cout << chromo << "\t" << start << "\t" << stop;
  cout << "\t" << lines.size() << "\t";
  cout << "." << "\t" << "." << "\t" << "." << "\t" << "." << "\t";
  cout << "GT:RD:PD:MP";
//  cout << "GT:RD:P1,P2,P3,P4";


  for(int i=0; i<nsamps; i++){
    cout << "\t" << nGTs[i] << ":" << RDs[i];
    cout << ":" << P1[i] << "," << P2[i] << "," << P3[i] << "," << P4[i];
    cout << "," << P5[i] << "," << P6[i] << "," << P7[i];

    /* Find majority rule ploidy. */
    int mjploid = 0;
    if(P1[i] + P2[i] + P3[i] + P4[i] + P5[i] + P6[i] + P7[i] > 0){mjploid=1;}
    int ploids[] = {P1[i], P2[i], P3[i], P4[i], P5[i], P6[i], P7[i]};
    if(P2[i] > ploids[mjploid]){mjploid=2;}
    if(P3[i] > ploids[mjploid]){mjploid=2;}
    if(P4[i] > ploids[mjploid]){mjploid=2;}
    if(P5[i] > ploids[mjploid]){mjploid=2;}
    if(P6[i] > ploids[mjploid]){mjploid=2;}
    if(P7[i] > ploids[mjploid]){mjploid=2;}

    cout << ":" << mjploid;

//    cout << "\t";
  }

  cout << "\n";
}

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
  cout << "GT:RD:PD";
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
    cout << "Error, this does not appear to be a properly formatted file:\n";
    cout << "Expecting a header beginning with '#'.\n";
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
  int last = start;

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
  std::vector <string> lines;
  lines.push_back(line);
  while(getline (cin,line)){
    split( fields, line, is_any_of( "\t" ) );

    if(atoi(fields[1].c_str()) <= stop){
      /* Position is within current window */
      lines.push_back(line);
    } else {
      /* Position is beyond current window */
      /* Process window. */
      proc_win(fields[0], start, stop, lines);

      /* Begin the next window. */
      last = stop;
      start = stop + 1;
      stop = start + win;

      /* Do we have empty windows? */
      if(atoi(fields[1].c_str()) > stop){
        int i = (atoi(fields[1].c_str())-stop)/(win+1);
        for(int j=0; j<=i; j++){
          print_null(fields[0], start, stop, snames);
          last = stop;
          start = stop + 1;
          stop = start + win;
        }
      }
      lines.clear();
      lines.push_back(line);
    }
  }
  return 0;
}

