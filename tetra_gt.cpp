// Compile with:
// g++ -std=c++0x gt.cpp -lgmpxx -lgmp -o gt
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <math.h> /* log10 */
#include <gmp.h> // GNU multiple precision arithmetic library, not standard: libgmp3-dev.
#include <boost/algorithm/string.hpp> // Not standard on Macs or Ubuntu: libboost1.46-dev.

using namespace std;
using namespace boost;


/* Functions */

void refA(string nucs, int sampn, int nuc_cnt[][8]){
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

void refC(string nucs, int sampn, int nuc_cnt[][8]){
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

void refG(string nucs, int sampn, int nuc_cnt[][8]){
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

void refT(string nucs, int sampn, int nuc_cnt[][8]){
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



/* Likelihoods */
void mult_pl(int pls[][26], int nsamp, float err, int nuc_cnts[][8], int min_cnt){
//  cout << "\nmult_pl nsamp: " << nsamp << "\n";
  for(int i=0; i<nsamp; i++){
    float mls [26]; // Maximum likelihoods.
    mpz_t fac [5]; // Factorials n, A, C, G, T.
    mpf_t pos [3]; // Floats for num, den, pos;
    int nuc_cnt[4];
    nuc_cnt[0] = nuc_cnts[i][0] + nuc_cnts[i][1];
    nuc_cnt[1] = nuc_cnts[i][2] + nuc_cnts[i][3];
    nuc_cnt[2] = nuc_cnts[i][4] + nuc_cnts[i][5];
    nuc_cnt[3] = nuc_cnts[i][6] + nuc_cnts[i][7];
    int rd = nuc_cnt[0] + nuc_cnt[1] + nuc_cnt[2] + nuc_cnt[3];

    /* Subtract minimum threshold before 
       genotype calling.  */
    if(nuc_cnt[0] >= min_cnt){
      nuc_cnt[0] = nuc_cnt[0] - min_cnt;
    } else {
      nuc_cnt[0] = 0;
    }
    if(nuc_cnt[1] >= min_cnt){
      nuc_cnt[1] = nuc_cnt[1] - min_cnt;
    } else {
      nuc_cnt[1] = 0;
    }
    if(nuc_cnt[2] >= min_cnt){
      nuc_cnt[2] = nuc_cnt[2] - min_cnt;
    } else {
      nuc_cnt[2] = 0;
    }
    if(nuc_cnt[3] >= min_cnt){
      nuc_cnt[3] = nuc_cnt[3] - min_cnt;
    } else {
      nuc_cnt[3] = 0;
    }

//    cout << "\nrd = " << rd << "\n";

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

    /* Count section. */
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
    double posd = mpf_get_d(pos[0]);

    /* Homozygotes. */
    // AA, CC, GG, TT.
    mls[0] = posd * pow(1-(3*err)/4, nuc_cnt[0]) * pow(err/4, nuc_cnt[1]+nuc_cnt[2]+nuc_cnt[3]);
    mls[1] = posd * pow(1-(3*err)/4, nuc_cnt[1]) * pow(err/4, nuc_cnt[0]+nuc_cnt[2]+nuc_cnt[3]);
    mls[2] = posd * pow(1-(3*err)/4, nuc_cnt[2]) * pow(err/4, nuc_cnt[0]+nuc_cnt[1]+nuc_cnt[3]);
    mls[3] = posd * pow(1-(3*err)/4, nuc_cnt[3]) * pow(err/4, nuc_cnt[0]+nuc_cnt[1]+nuc_cnt[2]);

    // Biallelic heterozygotes.
    // AC, AG, AT, CG, CT, GT.
    mls[4] = posd * pow(0.5-err/4, nuc_cnt[0]+nuc_cnt[1]) * pow(err/4, nuc_cnt[2]+nuc_cnt[3]);
    mls[5] = posd * pow(0.5-err/4, nuc_cnt[0]+nuc_cnt[2]) * pow(err/4, nuc_cnt[1]+nuc_cnt[3]);
    mls[6] = posd * pow(0.5-err/4, nuc_cnt[0]+nuc_cnt[3]) * pow(err/4, nuc_cnt[1]+nuc_cnt[2]);
    mls[7] = posd * pow(0.5-err/4, nuc_cnt[1]+nuc_cnt[2]) * pow(err/4, nuc_cnt[0]+nuc_cnt[3]);
    mls[8] = posd * pow(0.5-err/4, nuc_cnt[1]+nuc_cnt[3]) * pow(err/4, nuc_cnt[0]+nuc_cnt[2]);
    mls[9] = posd * pow(0.5-err/4, nuc_cnt[2]+nuc_cnt[3]) * pow(err/4, nuc_cnt[0]+nuc_cnt[1]);

    // Biallelic triploid loci.
    // AAC, AAG, AAT, CCA, CCG, CCT, GGA, GGC, GGT, TTA, TTC, TTG.
    float prop3 = 0.333333;
    float prop6 = 0.666667;
    mls[10] = posd * pow(prop6-err/4, nuc_cnt[0]) * pow(prop3-err/4, nuc_cnt[1]) * pow(err/4, nuc_cnt[2]+nuc_cnt[3]);
    mls[11] = posd * pow(prop6-err/4, nuc_cnt[0]) * pow(prop3-err/4, nuc_cnt[2]) * pow(err/4, nuc_cnt[1]+nuc_cnt[3]);
    mls[12] = posd * pow(prop6-err/4, nuc_cnt[0]) * pow(prop3-err/4, nuc_cnt[3]) * pow(err/4, nuc_cnt[1]+nuc_cnt[2]);

    mls[13] = posd * pow(prop6-err/4, nuc_cnt[1]) * pow(prop3-err/4, nuc_cnt[0]) * pow(err/4, nuc_cnt[2]+nuc_cnt[3]);
    mls[14] = posd * pow(prop6-err/4, nuc_cnt[1]) * pow(prop3-err/4, nuc_cnt[2]) * pow(err/4, nuc_cnt[0]+nuc_cnt[3]);
    mls[15] = posd * pow(prop6-err/4, nuc_cnt[1]) * pow(prop3-err/4, nuc_cnt[3]) * pow(err/4, nuc_cnt[0]+nuc_cnt[2]);

    mls[16] = posd * pow(prop6-err/4, nuc_cnt[2]) * pow(prop3-err/4, nuc_cnt[0]) * pow(err/4, nuc_cnt[1]+nuc_cnt[3]);
    mls[17] = posd * pow(prop6-err/4, nuc_cnt[2]) * pow(prop3-err/4, nuc_cnt[1]) * pow(err/4, nuc_cnt[0]+nuc_cnt[3]);
    mls[18] = posd * pow(prop6-err/4, nuc_cnt[2]) * pow(prop3-err/4, nuc_cnt[3]) * pow(err/4, nuc_cnt[0]+nuc_cnt[1]);

    mls[19] = posd * pow(prop6-err/4, nuc_cnt[3]) * pow(prop3-err/4, nuc_cnt[0]) * pow(err/4, nuc_cnt[1]+nuc_cnt[2]);
    mls[20] = posd * pow(prop6-err/4, nuc_cnt[3]) * pow(prop3-err/4, nuc_cnt[1]) * pow(err/4, nuc_cnt[0]+nuc_cnt[2]);
    mls[21] = posd * pow(prop6-err/4, nuc_cnt[3]) * pow(prop3-err/4, nuc_cnt[2]) * pow(err/4, nuc_cnt[0]+nuc_cnt[1]);

    // Triallelic triploid loci.
    // ACG, ACT, AGT, CGT.
    mls[22] = posd * pow(prop3-err/4, nuc_cnt[0]+nuc_cnt[1]+nuc_cnt[2]) * pow(err/4, nuc_cnt[3]);
    mls[23] = posd * pow(prop3-err/4, nuc_cnt[0]+nuc_cnt[1]+nuc_cnt[3]) * pow(err/4, nuc_cnt[2]);
    mls[24] = posd * pow(prop3-err/4, nuc_cnt[0]+nuc_cnt[2]+nuc_cnt[3]) * pow(err/4, nuc_cnt[1]);
    mls[25] = posd * pow(prop3-err/4, nuc_cnt[1]+nuc_cnt[2]+nuc_cnt[3]) * pow(err/4, nuc_cnt[0]);

    /* Biallelic tetraploid loci */
    /* AAAC, AAAG, AAAT, CCCA, CCCG, CCCT, GGGA, GGGC, GGGT, TTTA, TTTC, TTTG. */


    /* Triallelic tetraploid loci */
    /* AACG, AACT, AAGT, CCAG, CCAT, CCGT, GGAC, GGAT, GGCT, TTAC, TTAG, TTCG */


    /* Tetra-allelic tetraploid loci */
    /* ACGT */




/*
    cout << "\n";
    cout << "\n";
    cout << "posd: " << posd << "\n";
    cout << "CC: " << pow(1-(3*err)/4, nuc_cnt[1]);
    cout << "\n";
    cout << "CCe: " << pow(err/4, nuc_cnt[0]+nuc_cnt[2]+nuc_cnt[3]);
    cout << "\n";

    cout << "Nuc cnt:";
    cout << nuc_cnt[0];
    for(int j=1; j<4; j++){cout << "," << nuc_cnt[j];}
    cout << "\n";
    cout << "A,C,G,T,AC,AG,AT,CG,CT,GT,AAC,AAG,AAT,CCA,CCG,CCT,GGA,GGC,GGT,TTA,TTC,TTG,ACG,ACT,AGT,CGT\n";
    cout << "MLS:\n";
    cout << mls[0];
    for (int j=1; j<26; j++){
      cout << "," << mls[j];
    }
*/

    // Phred scale the likelihoods.
    for (int j = 0; j < 26; j++){
      if (mls[j] == 0){
        mls[j] = 9999;
      } else {
        mls[j] = trunc(-10 * log10(mls[j]));
        if(mls[j] == -0){mls[j] = 0;}
      }
    }

//    cout << "\nGot here.\n";
    // Set parent array values.
//    cout << "\n";
    for (int j = 0; j < 26; j++){
      pls[i][j] = int(mls[j]);
//      cout << pls[i][j] << "\t" << mls[j] << "\n";
    }
    // Clean up the mpz_t handles or else we will leak memory
    for (int j=0; j<5; j++){mpz_clear(fac[j]);}
    for (int j=0; j<3; j++){mpf_clear(pos[j]);}
  }
}


/* Determine the genotype. */
void det_gt(string gts[], int nsamp, int rds[], int nuc_cnts[][8], int pls[][26]){
  for(int i=0; i<nsamp; i++){
    if(rds[i] > 0){
      int k = 0;
      for (int j = 1; j < 26; j++){
        if(pls[i][j] < pls[i][k]){k = j;}
      }

      if(k==0){gts[i] = "A/A";}
      if(k==1){gts[i] = "C/C";}
      if(k==2){gts[i] = "G/G";}
      if(k==3){gts[i] = "T/T";}
      if(k==4){gts[i] = "A/C";}
      if(k==5){gts[i] = "A/G";}
      if(k==6){gts[i] = "A/T";}
      if(k==7){gts[i] = "C/G";}
      if(k==8){gts[i] = "C/T";}
      if(k==9){gts[i] = "G/T";}

      if(k==10){gts[i] = "A/A/C";}
      if(k==11){gts[i] = "A/A/G";}
      if(k==12){gts[i] = "A/A/T";}
      if(k==13){gts[i] = "C/C/A";}
      if(k==14){gts[i] = "C/C/G";}
      if(k==15){gts[i] = "C/C/T";}
      if(k==16){gts[i] = "G/G/A";}
      if(k==17){gts[i] = "G/G/C";}
      if(k==18){gts[i] = "G/G/T";}
      if(k==19){gts[i] = "T/T/A";}
      if(k==20){gts[i] = "T/T/C";}
      if(k==21){gts[i] = "T/T/G";}

      if(k==22){gts[i] = "A/C/G";}
      if(k==23){gts[i] = "A/C/T";}
      if(k==24){gts[i] = "A/G/T";}
      if(k==25){gts[i] = "C/G/T";}
    }
  }
}

int cnt_gts(int nsamp, string gts[]){
  int unigt [26] = {};
  for(int i=0; i<nsamp; i++){
    if(gts[i]=="A/A"){unigt[0] = 1;}
    if(gts[i]=="C/C"){unigt[1] = 1;}
    if(gts[i]=="G/G"){unigt[2] = 1;}
    if(gts[i]=="T/T"){unigt[3] = 1;}

    if(gts[i]=="A/C"){unigt[4] = 1;}
    if(gts[i]=="A/G"){unigt[5] = 1;}
    if(gts[i]=="A/T"){unigt[6] = 1;}
    if(gts[i]=="C/G"){unigt[7] = 1;}
    if(gts[i]=="C/T"){unigt[8] = 1;}
    if(gts[i]=="G/T"){unigt[9] = 1;}

    if(gts[i]=="A/A/C"){unigt[10] = 1;}
    if(gts[i]=="A/A/G"){unigt[11] = 1;}
    if(gts[i]=="A/A/T"){unigt[12] = 1;}
    if(gts[i]=="C/C/A"){unigt[13] = 1;}
    if(gts[i]=="C/C/G"){unigt[14] = 1;}
    if(gts[i]=="C/C/T"){unigt[15] = 1;}
    if(gts[i]=="G/G/A"){unigt[16] = 1;}
    if(gts[i]=="G/G/C"){unigt[17] = 1;}
    if(gts[i]=="G/G/T"){unigt[18] = 1;}
    if(gts[i]=="T/T/A"){unigt[19] = 1;}
    if(gts[i]=="T/T/C"){unigt[20] = 1;}
    if(gts[i]=="T/T/G"){unigt[21] = 1;}

    if(gts[i]=="A/C/G"){unigt[22] = 1;}
    if(gts[i]=="A/C/T"){unigt[23] = 1;}
    if(gts[i]=="A/G/T"){unigt[24] = 1;}
    if(gts[i]=="C/G/T"){unigt[25] = 1;}
  }

  int cnt=0;
  for(int i=0; i<26; i++){cnt = cnt+unigt[i];}
  return(cnt);
}

void min_tresh(string gts[], int nsamp, int min_cnt, int nuc_cnts[][8]){
  for(int i=0; i<nsamp; i++){
    int nuc_cnt[4];
    nuc_cnt[0] = nuc_cnts[i][0] + nuc_cnts[i][1];
    nuc_cnt[1] = nuc_cnts[i][2] + nuc_cnts[i][3];
    nuc_cnt[2] = nuc_cnts[i][4] + nuc_cnts[i][5];
    nuc_cnt[3] = nuc_cnts[i][6] + nuc_cnts[i][7];

    /* Homozygotes */
    if(gts[i] == "A/A"){if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}}
    if(gts[i] == "C/C"){if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}}
    if(gts[i] == "G/G"){if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}}
    if(gts[i] == "T/T"){if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}}

    /* Diploid heterozygotes */
    if(gts[i] == "A/C"){
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "A/G"){
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "A/T"){
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "C/G"){
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "C/T"){
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "G/T"){
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
    }

    /* Triploid biallelic heterozygotes */
    if(gts[i] == "A/A/C"){
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "A/A/G"){
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "A/A/T"){
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "C/C/A"){
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "C/C/G"){
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "C/C/T"){
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "G/G/A"){
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "G/G/C"){
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "G/G/T"){
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "T/T/A"){
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "T/T/C"){
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "T/T/G"){
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
    }

    /* Triploid triallelic heterozygotes */
    if(gts[i] == "A/C/G"){
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "A/C/T"){
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "A/G/T"){
      if(nuc_cnt[0] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
    }
    if(gts[i] == "C/G/T"){
      if(nuc_cnt[1] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[2] < min_cnt){gts[i] = "./.";}
      if(nuc_cnt[3] < min_cnt){gts[i] = "./.";}
    }
//    cout << "\n";
//    cout << gts[i] << ":" << min_cnt << ":" << nuc_cnt[0] << "," << nuc_cnt[1] << "," << nuc_cnt[2] << "," << nuc_cnt[3] << "\n";
  }
}


/* Print functions */

void print_header(float error, int min_cnt, string sfile){
//  cout << sfile << "\n";
  cout << "##fileformat=VCFv4.2\n";
  cout << "##source=gtV0.0.0\n";
  cout << "##FILTER=<ID=NA,Description=\"Per nucleotide minimum threshold of " << min_cnt << "\">\n";
  cout << "##FILTER=<ID=NA,Description=\"Error rate of " << error << " used for likelihood calculation\">\n";
  cout << "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Read depth\">\n";
  cout << "##FORMAT=<ID=CT,Number=4,Type=Integer,Description=\"Count of each nucleotide (A,C,G,T)\">\n";
  cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
  cout << "##FORMAT=<ID=PL,Number=26,Type=Integer,Description=\"Phred scaled likelihood for 26 genotpyes\">\n";
  cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

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

  cout << "\n";
}

void print_fix(vector <string> fields){
  cout << fields[0] << "\t" << fields[1] << "\t" << ".";
  cout << "\t" << fields[2] << "\t" << ".";
  cout << "\t" << "." << "\t" << ".";
  cout << "\t" << ".";
}

void print_locus(vector <string> fields, int counts, int phred, int nsamp, int rds[], int nuc_cnts[][8], int pls[][26], string gts[]){
  /* Print fixed portion of locus. */
  print_fix(fields);

  cout << "\t";
  cout << "RD";
  if(counts == 1){cout << ":CT";}
  if(phred == 1){cout << ":PL";}
  cout << ":GT";

  for(int i=0; i<nsamp; i++){
    cout << "\t";
    // Read depth.
    cout << rds[i];
    // Counts.
    if(counts == 1){
      cout << ":" << nuc_cnts[i][0];
      for(int j=1; j<8; j++){cout << "," << nuc_cnts[i][j];}
    }
    // Phred-scaled likelihoods.
    if(phred == 1){
      cout << ":" << pls[i][0];
      for(int j=1; j<25; j++){cout << "," << pls[i][j];}
    }
    // Genotype.
    cout << ":" << gts[i];
  }
  cout << "\n";
}



/* Main */

int main(int argc, char **argv) {
  string lineInput;
  vector <string> fields;
  int opt, header = 0, counts = 0, phred = 0;
//  float error = 0.0;  // Can not be zero!!!
//  float error = 0.01;
//  float error = 0.000001;
  float error = 0.000000001;
  int min_cnt = 3; // Minimum threshold.
  string sfile = "NA";

  // Parse command line options.
  while ((opt = getopt(argc, argv, "ce:hm:ps:")) != -1) {
    switch (opt) {
      case 'c':
        counts = 1;
        break;
      case 'e':
        error = atof(optarg);
        break;
      case 'h':
        header = 1;
        break;
      case 'm':
        min_cnt = atoi(optarg);
        break;
      case 'p':
        phred = 1;
        break;
      case 's':
        sfile = optarg;
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s [-h]\n", argv[0]);
        exit(EXIT_FAILURE);
    }
  }

  /* Header. */
  if(header == 1){print_header(error, min_cnt, sfile);}

  /* Parse line by line. */
  while (getline(cin,lineInput)) {
    split( fields, lineInput, is_any_of( "\t " ) );
    int nsamp = (fields.size()-3)/3;  // Determine the number of samples.
    int nuc_cnts [nsamp][8]; // A,a,C,c,G,g,T,t.
    int rds [nsamp]; // Read depth.
    for(int i=0; i<nsamp; i++){rds[i] = 0;}
    string gts [nsamp]; // Genotypes.
//    for(int i=0; i<nsamp; i++){gts[i] = "NA";}
    for(int i=0; i<nsamp; i++){gts[i] = "./.";}

    /* Diploid PLs. */
//    int pls[] = {0,0,0,0,0,0,0,0,0,0}; // AA,CC,GG,TT,AC,AG,AT,CG,CT,GT.

    /* Triploid PLs. */
    /*
    AA,CC,GG,TT,AC,AG,AT,CG,CT,GT, : 10
    AAC,AAG,AAT,CCA,CCG,CCT, : 6
    GGA,GGC,GGT,TTA,TTC,TTG, : 6
    ACG,ACT,AGT,CGT. : 4
    */
    int pls [nsamp][26];

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
          if(fields[2] == "T"){refG(fields[i], sampn, nuc_cnts);}
        }
      }
    }

    /* Get read depths */
    get_rd(rds, nsamp, nuc_cnts);
//    get_rd(rds, sampn, nuc_cnts);

    /* Calculate Phred-scaled likelihoods */
    mult_pl(pls, nsamp, error, nuc_cnts, min_cnt);
//    mult_pl(pls, sampn, error, nuc_cnts);

    /* Determine a genotype */
    det_gt(gts, nsamp, rds, nuc_cnts, pls);

    /* Minimum threshold. */
    min_tresh(gts, nsamp, min_cnt, nuc_cnts);

    /* Print locus. */
//    vector <string> vals;
//    vals = unique(gts);

    int unique_gt = cnt_gts(nsamp, gts);
    if(unique_gt > 1){
      print_locus(fields, counts, phred, nsamp, rds, nuc_cnts, pls, gts);
    }

  }
  return 0;
}


