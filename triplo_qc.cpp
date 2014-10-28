// Compile with:
// g++ -std=c++0x qc.cpp -o qc
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <boost/algorithm/string.hpp> // Not standard on Macs or Ubuntu: libboost1.46-dev.

using namespace std;
using namespace boost;



void cntA(long long AAg[], long long AAe[], int samp, string counts){
  vector <string> cnts;
  split(cnts, counts, is_any_of(","));
  AAg[samp] = AAg[samp] + stoi(cnts[0]) + stoi(cnts[1]);
  AAe[samp] = AAe[samp] + stoi(cnts[2]) + stoi(cnts[3]) +
                          stoi(cnts[4]) + stoi(cnts[5]) +
                          stoi(cnts[6]) + stoi(cnts[7]);
}

void cntC(long long CCg[], long long CCe[], int samp, string counts){
  vector <string> cnts;
  split(cnts, counts, is_any_of(","));
  CCg[samp] = CCg[samp] + stoi(cnts[2]) + stoi(cnts[3]);
  CCe[samp] = CCe[samp] + stoi(cnts[0]) + stoi(cnts[1]) +
                          stoi(cnts[4]) + stoi(cnts[5]) +
                          stoi(cnts[6]) + stoi(cnts[7]);
}

void cntG(long long GGg[], long long GGe[], int samp, string counts){
  vector <string> cnts;
  split(cnts, counts, is_any_of(","));
  GGg[samp] = GGg[samp] + stoi(cnts[4]) + stoi(cnts[5]);
  GGe[samp] = GGe[samp] + stoi(cnts[0]) + stoi(cnts[1]) +
                          stoi(cnts[2]) + stoi(cnts[3]) +
                          stoi(cnts[6]) + stoi(cnts[7]);
}

void cntT(long long TTg[], long long TTe[], int samp, string counts){
  vector <string> cnts;
  split(cnts, counts, is_any_of(","));
  TTg[samp] = TTg[samp] + stoi(cnts[6]) + stoi(cnts[7]);
  TTe[samp] = TTe[samp] + stoi(cnts[0]) + stoi(cnts[1]) +
                          stoi(cnts[2]) + stoi(cnts[3]) +
                          stoi(cnts[4]) + stoi(cnts[5]);
}

void cnt_homo(long long nHo[], int samp, string GT){
  if(GT == "A/A"){nHo[samp]++;}
  if(GT == "C/C"){nHo[samp]++;}
  if(GT == "G/G"){nHo[samp]++;}
  if(GT == "T/T"){nHo[samp]++;}
}

void cnt_het(long long nHe[], int samp, string GT){
  if(GT == "A/C"){nHe[samp]++;}
  if(GT == "A/G"){nHe[samp]++;}
  if(GT == "A/T"){nHe[samp]++;}
  if(GT == "C/G"){nHe[samp]++;}
  if(GT == "C/T"){nHe[samp]++;}
  if(GT == "G/T"){nHe[samp]++;}
}


void cnt_tri(long long int ntrip[], int samp, string GT){
    if(GT == "A/A/C"){ntrip[samp]++;}
    if(GT == "A/A/G"){ntrip[samp]++;}
    if(GT == "A/A/T"){ntrip[samp]++;}
    if(GT == "C/C/A"){ntrip[samp]++;}
    if(GT == "C/C/G"){ntrip[samp]++;}
    if(GT == "C/C/T"){ntrip[samp]++;}
    if(GT == "G/G/A"){ntrip[samp]++;}
    if(GT == "G/G/C"){ntrip[samp]++;}
    if(GT == "G/G/T"){ntrip[samp]++;}
    if(GT == "T/T/A"){ntrip[samp]++;}
    if(GT == "T/T/C"){ntrip[samp]++;}
    if(GT == "T/T/G"){ntrip[samp]++;}

    if(GT == "A/C/G"){ntrip[samp]++;}
    if(GT == "A/C/T"){ntrip[samp]++;}
    if(GT == "A/G/T"){ntrip[samp]++;}
    if(GT == "C/G/T"){ntrip[samp]++;}
}

void print_out(int nsamp, long long rds[], long long nGT[], long long nNA[], long long AAg[], long long AAe[], long long CCg[], long long CCe[], long long GGg[], long long GGe[], long long TTg[], long long TTe[], long long nHo[], long long nHe[], long long ntrip[]){
  /* Print header */

/*
  cout << "RD" << "\t";
  cout << "nGT" << "\t" << "nNA" << "\t" << "nHo" << "\t" << "nHe";
  cout << "\t" << "AAc" << "\t" << "AAe" << "\t" << "CCc" << "\t" << "CCe";
  cout << "\t" << "GGc" << "\t" << "GGe" << "\t" << "TTc" << "\t" << "TTe";
  cout << "\t" << "nTrip";
  cout << "\n";
*/

  cout << "RD";
  for(int i=0; i<nsamp; i++){cout << "\t" << rds[i];}
  cout << "\n";

  cout << "nGT";
  for(int i=0; i<nsamp; i++){cout << "\t" << nGT[i];}
  cout << "\n";

  cout << "nNA";
  for(int i=0; i<nsamp; i++){cout << "\t" << nNA[i];}
  cout << "\n";

  cout << "AAg";
  for(int i=0; i<nsamp; i++){cout << "\t" << AAg[i];}
  cout << "\n";

  cout << "AAe";
  for(int i=0; i<nsamp; i++){cout << "\t" << AAe[i];}
  cout << "\n";

  cout << "CCg";
  for(int i=0; i<nsamp; i++){cout << "\t" << CCg[i];}
  cout << "\n";

  cout << "CCe";
  for(int i=0; i<nsamp; i++){cout << "\t" << CCe[i];}
  cout << "\n";

  cout << "GGg";
  for(int i=0; i<nsamp; i++){cout << "\t" << GGg[i];}
  cout << "\n";

  cout << "GGe";
  for(int i=0; i<nsamp; i++){cout << "\t" << GGe[i];}
  cout << "\n";

  cout << "TTg";
  for(int i=0; i<nsamp; i++){cout << "\t" << TTg[i];}
  cout << "\n";

  cout << "TTe";
  for(int i=0; i<nsamp; i++){cout << "\t" << TTe[i];}
  cout << "\n";

  cout << "nHo";
  for(int i=0; i<nsamp; i++){cout << "\t" << nHo[i];}
  cout << "\n";

  cout << "nHe";
  for(int i=0; i<nsamp; i++){cout << "\t" << nHe[i];}
  cout << "\n";

  cout << "ntrip";
  for(int i=0; i<nsamp; i++){cout << "\t" << ntrip[i];}
  cout << "\n";

}


/* Main */

int main(){
  string lineInput;
  string test;
  vector <string> fields;
  vector <string> format;
  int nsites = 0;

/*
  while(getline(cin, lineInput)){
    std::string test = lineInput.substr(0, 1);
    if (!(test == "#"))
    {
//        ++numlines;
    }
    else break;
  }
*/




  /* Initialize with the first line. */
  getline(cin, lineInput);
  split( fields, lineInput, is_any_of( "\t " ) );
  int nsamp = fields.size()-9;

//  long long int stats [nsamp][12]; // nGT, nNA, nHo, nHe, AAg, AAe, CCg, CCe, GGg, GGe, TTg, TTe.
  long long rds[nsamp];
  long long ntrip[nsamp];
  long long nGT[nsamp]; 
  long long nNA[nsamp];
  long long nHo[nsamp];
  long long nHe[nsamp];
  long long AAg[nsamp];
  long long AAe[nsamp];
  long long CCg[nsamp];
  long long CCe[nsamp];
  long long GGg[nsamp];
  long long GGe[nsamp];
  long long TTg[nsamp];
  long long TTe[nsamp];

  /* Initialize to zero */
  for(int i=0; i<nsamp; i++){
    rds[i] = 0;
    ntrip[i] = 0;
    nGT[i] = 0;
    nNA[i] = 0;
    nHo[i] = 0;
    nHe[i] = 0;
    AAg[i] = 0;
    AAe[i] = 0;
    CCg[i] = 0;
    CCe[i] = 0;
    GGg[i] = 0;
    GGe[i] = 0;
    TTg[i] = 0;
    TTe[i] = 0;
  }

  /* Determine format positions */
  int RD = 0, CT = 0, GT = 0;
  split( format, fields[8], is_any_of( ":" ) );
  for(int i=0; i<format.size(); i++){
    if(format[i] == "RD"){RD = i;}
    if(format[i] == "CT"){CT = i;}
    if(format[i] == "GT"){GT = i;}
  }

  /* Iterate over samples */
  /* Samples start on line 9 of VCF files */

  for(int i=9; i < fields.size(); i++){ // Samples.
//    for(int j=0; j<12; j++){stats[i-9][j] = 0;} // Initialize to zero.
    split( format, fields[i], is_any_of( ":" ) );
    cnt_homo(nHo, i-9, format[GT]);
    cnt_het(nHe, i-9, format[GT]);
    if(format[GT] == "./."){nNA[i-9]++;}
//    if(format[GT] == "./."){nNA[i-9]++;}
    cnt_tri(ntrip, i-9, format[GT]);
    rds[i-9] = rds[i-9] + stoi(format[RD]);
    if(format[GT] == "A/A"){cntA(AAg, AAe, i-9, format[CT]);}
    if(format[GT] == "C/C"){cntC(CCg, CCe, i-9, format[CT]);}
    if(format[GT] == "G/G"){cntG(GGg, GGe, i-9, format[CT]);}
    if(format[GT] == "T/T"){cntT(TTg, TTe, i-9, format[CT]);}


//);}
  }


  /* Parse line by line. */
  while (getline(cin, lineInput)) {
    split( fields, lineInput, is_any_of( "\t " ) );

    /* Iterate over samples */
    for(int i=9; i < fields.size(); i++){
      split( format, fields[i], is_any_of( ":" ) );
      cnt_homo(nHo, i-9, format[GT]);
      cnt_het(nHe, i-9, format[GT]);
      if(format[GT] == "./."){nNA[i-9]++;}
      if(format[GT] == "A/A"){cntA(AAg, AAe, i-9, format[CT]);}
      if(format[GT] == "C/C"){cntC(CCg, CCe, i-9, format[CT]);}
      if(format[GT] == "G/G"){cntG(GGg, GGe, i-9, format[CT]);}
      if(format[GT] == "T/T"){cntT(TTg, TTe, i-9, format[CT]);}
      cnt_tri(ntrip, i-9, format[GT]);
      rds[i-9] = rds[i-9] + stoi(format[RD]);
      }
  }


  for(int i=0; i<nsamp; i++){
    nGT[i] = nHo[i] + nHe[i] + ntrip[i];
  }

  /* Print output */
  print_out(nsamp, rds, nGT, nNA, AAg, AAe, CCg, CCe, GGg, GGe, TTg, TTe, nHo, nHe, ntrip);

  /*
  for(int i=0; i<nsamp; i++){ // Sample
    cout << rds[i] << "\t";
    cout << stats[i][0];
    for(int j=1; j<12; j++){ // Stat
      cout << "\t" << stats[i][j];
    }
    cout << "\t" << ntrip[i];
//    cout << "\t" << stats[i][5] / (stats[i][4] + stats[i][5]);
    cout << "\n";
  }
  */

  return 0;
}


