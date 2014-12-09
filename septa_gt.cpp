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

void int_to_nucs(string gts[27], vector<pair<int,int>> moves){

  /* Convert pairs to array of ints. */
  int nucs[4];
  for(int i=0; i<4; i++){
    nucs[i] = moves[i].second;
  }

//  cout << nucs[0] << "," << nucs[1] << "," << nucs[2] << "," << nucs[3];
//  cout << "\n";
//  cout << nucs[0] << "," << nucs[1] << "," << nucs[2] << "," << nucs[3];
//  cout << "\n";

  /* Array indices */
  /* 0 = Homozyote 1/0/0/0 */
  /* 1 = Bi-allelic heterozyote 1/1 */
  /* 2 = Bi-allelic triploid 2/1 */
  /* 3 = Tri-allelic triploid 1/1/1 */

  /* 4 = Bi-allelic tetraploid loci 3/1/0/0 */
  /* 5 = Tri-allelic tetraploid loci 2/1/1/0 */
  /* 6 = Tetra-allelic tetraploid loci 1/1/1/1 */

  /* 7 = Bi-allelic pentaploid 4/1/0/0 */
  /* 8 = Bi-allelic pentaploid 3/2/0/0 */
  /* 9 = Tri-allelic pentaploid 3/1/1/0 */
  /* 10 = Tri-allelic pentaploid 2/2/1/0 */
  /* 11 = Tetra-allelic pentaploid 2/1/1/1 */

  /* 12 = Bi-allelic hexaploid loci */
  /* 13 = Tri-allelic hexaploid loci 4/1/1 */
  /* 14 = Tri-allelic hexaploid loci 3/2/1 */
  /* 15 = Tetra-allelic hexaploid loci 3/1/1/1 */
  /* 16 = Tetra-allelic hexaploid loci 2/2/1/1 */

  /* 17 = Bi-allelic septaploid loci 6/1 */
  /* 18 = Bi-allelic septaploid loci 5/2 */
  /* 19 = Bi-allelic septaploid loci 4/3 */
  /* 20 = Tri-allelic septaploid loci 5/1/1 */
  /* 21 = Tri-allelic septaploid loci 4/2/1 */
  /* 22 = Tri-allelic septaploid loci 3/3/1 */
  /* 23 = Tri-allelic septaploid loci 3/2/2 */

  /* 24 = Tetra-allelic septaploid loci 4/1/1/1 */
  /* 25 = Tetra-allelic septaploid loci 3/2/1/1 */
  /* 26 = Tetra-allelic septaploid loci 2/2/2/1 */

  if(nucs[0] == 0){
    gts[0] = "A/A";
    if(nucs[1] == 1){
      /* A/C Bi-allelic */
      gts[1]  = "A/C";
      gts[2]  = "A/A/C";
      gts[4]  = "A/A/A/C";
      gts[7]  = "A/A/A/A/C";
      gts[8]  = "A/A/A/C/C";
      gts[12] = "A/A/A/A/A/C";
      gts[17] = "A/A/A/A/A/A/C";
      gts[18] = "A/A/A/A/A/C/C";
      gts[19] = "A/A/A/A/C/C/C";
      if(nucs[2] == 2){
        /* A/C/G Tri-allelic */
        gts[3]  = "A/C/G";
        gts[5]  = "A/A/C/G";
        gts[9]  = "A/A/A/C/G";
        gts[10] = "A/A/C/C/G";
        gts[13] = "A/A/A/A/C/G";
        gts[14] = "A/A/A/C/C/G";
        gts[20] = "A/A/A/A/A/C/G";
        gts[21] = "A/A/A/A/C/C/G";
        gts[22] = "A/A/A/C/C/C/G";
        gts[23] = "A/A/A/C/C/G/G";
//        if(nucs[3]==3){
          /* A/C/G/T Tetra-allelic */
          gts[6]  = "A/C/G/T";
          gts[11] = "A/A/C/G/T";
          gts[15] = "A/A/A/C/G/T";
          gts[16] = "A/A/C/C/G/T";
          gts[24] = "A/A/A/A/C/G/T";
          gts[25] = "A/A/A/C/C/G/T";
          gts[26] = "A/A/C/C/G/G/T";
//        }
      } else if(nucs[2] == 3){
        /* A/C/T Tri-allelic */
        gts[3]  = "A/C/T";
        gts[5]  = "A/A/C/T";
        gts[9]  = "A/A/A/C/T";
        gts[10] = "A/A/C/C/T";
        gts[13] = "A/A/A/A/C/T";
        gts[14] = "A/A/A/C/C/T";
        gts[20] = "A/A/A/A/A/C/T";
        gts[21] = "A/A/A/A/C/C/T";
        gts[22] = "A/A/A/C/C/C/T";
        gts[23] = "A/A/A/C/C/T/T";
          /* A/C/T/G Tetra-allelic */
          gts[6]  = "A/C/T/G";
          gts[11] = "A/A/C/T/G";
          gts[15] = "A/A/A/C/T/G";
          gts[16] = "A/A/C/C/T/G";
          gts[24] = "A/A/A/A/C/T/G";
          gts[25] = "A/A/A/C/C/T/G";
          gts[26] = "A/A/C/C/T/T/G";
      }
    } else if(nucs[1] == 2){
      /* Bi-allelic */
      gts[1]  = "A/G";
      gts[2]  = "A/A/G";
      gts[4]  = "A/A/A/G";
      gts[7]  = "A/A/A/A/G";
      gts[8]  = "A/A/A/G/G";
      gts[12] = "A/A/A/A/A/G";
      gts[17] = "A/A/A/A/A/A/G";
      gts[18] = "A/A/A/A/A/G/G";
      gts[19] = "A/A/A/A/G/G/G";
      if(nucs[2] == 1){
        /* Tri-allelic */
        gts[3]  = "A/G/C";
        gts[5]  = "A/A/G/C";
        gts[9]  = "A/A/A/G/C";
        gts[10] = "A/A/G/G/C";
        gts[13] = "A/A/A/A/G/C";
        gts[14] = "A/A/A/G/G/C";
        gts[20] = "A/A/A/A/A/G/C";
        gts[21] = "A/A/A/A/G/G/C";
        gts[22] = "A/A/A/G/G/G/C";
        gts[23] = "A/A/A/G/G/C/C";
          /* Tetra-allelic */
          gts[6]  = "A/G/C/T";
          gts[11] = "A/A/G/C/T";
          gts[15] = "A/A/A/G/C/T";
          gts[16] = "A/A/G/G/C/T";
          gts[24] = "A/A/A/A/G/C/T";
          gts[25] = "A/A/A/G/G/C/T";
          gts[26] = "A/A/G/G/C/C/T";
      } else if(nucs[2] == 3){
        /* Tri-allelic */
        gts[3]  = "A/G/T";
        gts[5]  = "A/A/G/T";
        gts[9]  = "A/A/A/G/T";
        gts[10] = "A/A/G/G/T";
        gts[13] = "A/A/A/A/G/T";
        gts[14] = "A/A/A/G/G/T";
        gts[20] = "A/A/A/A/A/G/T";
        gts[21] = "A/A/A/A/G/G/T";
        gts[22] = "A/A/A/G/G/G/T";
        gts[23] = "A/A/A/G/G/T/T";
          /* Tetra-allelic */
          gts[6]  = "A/G/T/C";
          gts[11] = "A/A/G/T/C";
          gts[15] = "A/A/A/G/T/C";
          gts[16] = "A/A/G/G/T/C";
          gts[24] = "A/A/A/A/G/T/C";
          gts[25] = "A/A/A/G/G/T/C";
          gts[26] = "A/A/G/G/T/T/C";
      }
    } else if(nucs[1] == 3){
      /* Bi-allelic */
      gts[1]  = "A/T";
      gts[2]  = "A/A/T";
      gts[4]  = "A/A/A/T";
      gts[7]  = "A/A/A/A/T";
      gts[8]  = "A/A/A/T/T";
      gts[12] = "A/A/A/A/A/T";
      gts[17] = "A/A/A/A/A/A/T";
      gts[18] = "A/A/A/A/A/T/T";
      gts[19] = "A/A/A/A/T/T/T";
      if(nucs[2] == 1){
        /* Tri-allelic */
        gts[3]  = "A/T/C";
        gts[5]  = "A/A/T/C";
        gts[9]  = "A/A/A/T/C";
        gts[10] = "A/A/T/T/C";
        gts[13] = "A/A/A/A/T/C";
        gts[14] = "A/A/A/T/T/C";
        gts[20] = "A/A/A/A/A/T/C";
        gts[21] = "A/A/A/A/T/T/C";
        gts[22] = "A/A/A/T/T/T/C";
        gts[23] = "A/A/A/T/T/C/C";
          /* Tetra-allelic */
          gts[6]  = "A/T/C/G";
          gts[11] = "A/A/T/C/G";
          gts[15] = "A/A/A/T/C/G";
          gts[16] = "A/A/T/T/C/G";
          gts[24] = "A/A/A/A/T/C/G";
          gts[25] = "A/A/A/T/T/C/G";
          gts[26] = "A/A/T/T/C/C/G";
      } else if(nucs[2] == 2){
        /* Tri-allelic */
        gts[3]  = "A/T/G";
        gts[5]  = "A/A/T/G";
        gts[9]  = "A/A/A/T/G";
        gts[10] = "A/A/T/T/G";
        gts[13] = "A/A/A/A/T/G";
        gts[14] = "A/A/A/T/T/G";
        gts[20] = "A/A/A/A/A/T/G";
        gts[21] = "A/A/A/A/T/T/G";
        gts[22] = "A/A/A/T/T/T/G";
        gts[23] = "A/A/A/T/T/G/G";
          /* Tetra-allelic */
          gts[6]  = "A/T/G/C";
          gts[11] = "A/A/T/G/C";
          gts[15] = "A/A/A/T/G/C";
          gts[16] = "A/A/T/T/G/C";
          gts[24] = "A/A/A/A/T/G/C";
          gts[25] = "A/A/A/T/T/G/C";
          gts[26] = "A/A/T/T/G/G/C";
      }
    }
  } else if (nucs[0] == 1) {
    gts[0] = "C/C";
    if(nucs[1] == 0){
      /* Bi-allelic */
      gts[1]  = "C/A";
      gts[2]  = "C/C/A";
      gts[4]  = "C/C/C/A";
      gts[7]  = "C/C/C/C/A";
      gts[8]  = "C/C/C/A/A";
      gts[12] = "C/C/C/C/C/A";
      gts[17] = "C/C/C/C/C/C/A";
      gts[18] = "C/C/C/C/C/A/A";
      gts[19] = "C/C/C/C/A/A/A";
      if(nucs[2] == 2){
        /* Tri-allelic */
        gts[3]  = "C/A/G";
        gts[5]  = "C/C/A/G";
        gts[9]  = "C/C/C/A/G";
        gts[10] = "C/C/A/A/G";
        gts[13] = "C/C/C/C/A/G";
        gts[14] = "C/C/C/A/A/G";
        gts[20] = "C/C/C/C/C/A/G";
        gts[21] = "C/C/C/C/A/A/G";
        gts[22] = "C/C/C/A/A/A/G";
        gts[23] = "C/C/C/A/A/G/G";
          /* Tetra-allelic */
          gts[6]  = "C/A/G/T";
          gts[11] = "C/C/A/G/T";
          gts[15] = "C/C/C/A/G/T";
          gts[16] = "C/C/A/A/G/T";
          gts[24] = "C/C/C/C/A/G/T";
          gts[25] = "C/C/C/A/A/G/T";
          gts[26] = "C/C/A/A/G/G/T";
      } else if(nucs[2] == 3){
        /* Tri-allelic */
        gts[3]  = "C/A/T";
        gts[5]  = "C/C/A/T";
        gts[9]  = "C/C/C/A/T";
        gts[10] = "C/C/A/A/T";
        gts[13] = "C/C/C/C/A/T";
        gts[14] = "C/C/C/A/A/T";
        gts[20] = "C/C/C/C/C/A/T";
        gts[21] = "C/C/C/C/A/A/T";
        gts[22] = "C/C/C/A/A/A/T";
        gts[23] = "C/C/C/A/A/T/T";
          /* Tetra-allelic */
          gts[6]  = "C/A/T/G";
          gts[11] = "C/C/A/T/G";
          gts[15] = "C/C/C/A/T/G";
          gts[16] = "C/C/A/A/T/G";
          gts[24] = "C/C/C/C/A/T/G";
          gts[25] = "C/C/C/A/A/T/G";
          gts[26] = "C/C/A/A/T/T/G";
      }
    } else if(nucs[1] == 2){
      /* Bi-allelic */
      gts[1]  = "C/G";
      gts[2]  = "C/C/G";
      gts[4]  = "C/C/C/G";
      gts[7]  = "C/C/C/C/G";
      gts[8]  = "C/C/C/G/G";
      gts[12] = "C/C/C/C/C/G";
      gts[17] = "C/C/C/C/C/C/G";
      gts[18] = "C/C/C/C/C/G/G";
      gts[19] = "C/C/C/C/G/G/G";
      if(nucs[2] == 0){
        /* Tri-allelic */
        gts[3]  = "C/G/A";
        gts[5]  = "C/C/G/A";
        gts[9]  = "C/C/C/G/A";
        gts[10] = "C/C/G/G/A";
        gts[13] = "C/C/C/C/G/A";
        gts[14] = "C/C/C/G/G/A";
        gts[20] = "C/C/C/C/C/G/A";
        gts[20] = "C/C/C/C/C/G/A";
        gts[21] = "C/C/C/C/G/G/A";
        gts[22] = "C/C/C/G/G/G/A";
        gts[23] = "C/C/C/G/G/A/A";
          /* Tetra-allelic */
          gts[6]  = "C/G/A/T";
          gts[11] = "C/C/G/A/T";
          gts[15] = "C/C/C/G/A/T";
          gts[16] = "C/C/G/G/A/T";
          gts[24] = "C/C/C/C/G/A/T";
          gts[25] = "C/C/C/G/G/A/T";
          gts[26] = "C/C/G/G/A/A/T";
      } else if(nucs[2] == 3){
        /* Tri-allelic */
        gts[3]  = "C/G/T";
        gts[5]  = "C/C/G/T";
        gts[9]  = "C/C/C/G/T";
        gts[10] = "C/C/G/G/T";
        gts[13] = "C/C/C/C/G/T";
        gts[14] = "C/C/C/G/G/T";
        gts[20] = "C/C/C/C/C/G/T";
        gts[21] = "C/C/C/C/G/G/T";
        gts[22] = "C/C/C/G/G/G/T";
        gts[23] = "C/C/C/G/G/T/T";
          /* Tetra-allelic */
          gts[6]  = "C/G/T/A";
          gts[11] = "C/C/G/T/A";
          gts[15] = "C/C/C/G/T/A";
          gts[16] = "C/C/G/G/T/A";
          gts[24] = "C/C/C/C/G/T/A";
          gts[25] = "C/C/C/G/G/T/A";
          gts[26] = "C/C/G/G/T/T/A";
      }
    } else if(nucs[1] == 3){
      /* Bi-allelic */
      gts[1]  = "C/T";
      gts[2]  = "C/C/T";
      gts[4]  = "C/C/C/T";
      gts[7]  = "C/C/C/C/T";
      gts[8]  = "C/C/C/T/T";
      gts[12] = "C/C/C/C/C/T";
      gts[17] = "C/C/C/C/C/C/T";
      gts[18] = "C/C/C/C/C/T/T";
      gts[19] = "C/C/C/C/T/T/T";
      if(nucs[2] == 0){
        /* Tri-allelic */
        gts[3]  = "C/T/A";
        gts[5]  = "C/C/T/A";
        gts[9]  = "C/C/C/T/A";
        gts[10] = "C/C/T/T/A";
        gts[13] = "C/C/C/C/T/A";
        gts[14] = "C/C/C/T/T/A";
        gts[20] = "C/C/C/C/C/T/A";
        gts[21] = "C/C/C/C/T/T/A";
        gts[22] = "C/C/C/T/T/T/A";
        gts[23] = "C/C/C/T/T/A/A";
          /* Tetra-allelic */
          gts[6]  = "C/T/A/G";
          gts[11] = "C/C/T/A/G";
          gts[15] = "C/C/C/T/A/G";
          gts[16] = "C/C/T/T/A/G";
          gts[24] = "C/C/C/C/T/A/G";
          gts[25] = "C/C/C/T/T/A/G";
          gts[26] = "C/C/T/T/A/A/G";
      } else if(nucs[2] == 2){
        /* Tri-allelic */
        gts[3]  = "C/T/G";
        gts[5]  = "C/C/T/G";
        gts[9]  = "C/C/C/T/G";
        gts[10] = "C/C/T/T/G";
        gts[13] = "C/C/C/C/T/G";
        gts[14] = "C/C/C/T/T/G";
        gts[20] = "C/C/C/C/C/T/G";
        gts[21] = "C/C/C/C/T/T/G";
        gts[22] = "C/C/C/T/T/T/G";
        gts[23] = "C/C/C/T/T/G/G";
          /* Tetra-allelic */
          gts[6]  = "C/T/G/A";
          gts[11] = "C/C/T/G/A";
          gts[15] = "C/C/C/T/G/A";
          gts[16] = "C/C/T/T/G/A";
          gts[24] = "C/C/C/C/T/G/A";
          gts[25] = "C/C/C/T/T/G/A";
          gts[26] = "C/C/T/T/G/G/A";
      }
    }
  } else if (nucs[0] == 2) {
    gts[0] = "G/G";
    if(nucs[1] == 0){
      /* Bi-allelic */
      gts[1]  = "G/A";
      gts[2]  = "G/G/A";
      gts[4]  = "G/G/G/A";
      gts[7]  = "G/G/G/G/A";
      gts[8]  = "G/G/G/A/A";
      gts[12] = "G/G/G/G/G/A";
      gts[17] = "G/G/G/G/G/G/A";
      gts[18] = "G/G/G/G/G/A/A";
      gts[19] = "G/G/G/G/A/A/A";
      if(nucs[2] == 1){
        /* Tri-allelic */
        gts[3]  = "G/A/C";
        gts[5]  = "G/G/A/C";
        gts[9]  = "G/G/G/A/C";
        gts[10] = "G/G/A/A/C";
        gts[13] = "G/G/G/G/A/C";
        gts[14] = "G/G/G/A/A/C";
        gts[20] = "G/G/G/G/G/A/C";
        gts[21] = "G/G/G/G/A/A/C";
        gts[22] = "G/G/G/A/A/A/C";
        gts[23] = "G/G/G/A/A/C/C";
          /* Tetra-allelic */
          gts[6]  = "G/A/C/T";
          gts[11] = "G/G/A/C/T";
          gts[15] = "G/G/G/A/C/T";
          gts[16] = "G/G/A/A/C/T";
          gts[24] = "G/G/G/G/A/C/T";
          gts[25] = "G/G/G/A/A/C/T";
          gts[26] = "G/G/A/A/C/C/T";
      } else if(nucs[2] == 3){
        /* Tri-allelic */
        gts[3]  = "G/A/T";
        gts[5]  = "G/G/A/T";
        gts[9]  = "G/G/G/A/T";
        gts[10] = "G/G/A/A/T";
        gts[13] = "G/G/G/G/A/T";
        gts[14] = "G/G/G/A/A/T";
        gts[20] = "G/G/G/G/G/A/T";
        gts[21] = "G/G/G/G/A/A/T";
        gts[22] = "G/G/G/A/A/A/T";
        gts[23] = "G/G/G/A/A/T/T";
          /* Tetra-allelic */
          gts[6]  = "G/A/T/C";
          gts[11] = "G/G/A/T/C";
          gts[15] = "G/G/G/A/T/C";
          gts[16] = "G/G/A/A/T/C";
          gts[24] = "G/G/G/G/A/T/C";
          gts[25] = "G/G/G/A/A/T/C";
          gts[26] = "G/G/A/A/T/T/C";
      }
    } else if(nucs[1] == 1){
      /* Bi-allelic */
      gts[1]  = "G/C";
      gts[2]  = "G/G/C";
      gts[4]  = "G/G/G/C";
      gts[7]  = "G/G/G/G/C";
      gts[8]  = "G/G/G/C/C";
      gts[12] = "G/G/G/G/G/C";
      gts[17] = "G/G/G/G/G/G/C";
      gts[18] = "G/G/G/G/G/C/C";
      gts[19] = "G/G/G/G/C/C/C";
      if(nucs[2] == 0){
        /* Tri-allelic */
        gts[3]  = "G/C/A";
        gts[5]  = "G/G/C/A";
        gts[9]  = "G/G/G/C/A";
        gts[10] = "G/G/C/C/A";
        gts[13] = "G/G/G/G/C/A";
        gts[14] = "G/G/G/C/C/A";
        gts[20] = "G/G/G/G/G/C/A";
        gts[21] = "G/G/G/G/C/C/A";
        gts[22] = "G/G/G/C/C/C/A";
        gts[23] = "G/G/G/C/C/A/A";
          /* Tetra-allelic */
          gts[6]  = "G/C/A/T";
          gts[11] = "G/G/C/A/T";
          gts[15] = "G/G/G/C/A/T";
          gts[16] = "G/G/C/C/A/T";
          gts[24] = "G/G/G/G/C/A/T";
          gts[25] = "G/G/G/C/C/A/T";
          gts[26] = "G/G/C/C/A/A/T";
      } else if(nucs[2] == 3){
        /* Tri-allelic */
        gts[3]  = "G/C/T";
        gts[5]  = "G/G/C/T";
        gts[9]  = "G/G/G/C/T";
        gts[10] = "G/G/C/C/T";
        gts[13] = "G/G/G/G/C/T";
        gts[14] = "G/G/G/C/C/T";
        gts[20] = "G/G/G/G/G/C/T";
        gts[21] = "G/G/G/G/C/C/T";
        gts[22] = "G/G/G/C/C/C/T";
        gts[23] = "G/G/G/C/C/T/T";
          /* Tetra-allelic */
          gts[6]  = "G/C/T/A";
          gts[11] = "G/G/C/T/A";
          gts[15] = "G/G/G/C/T/A";
          gts[16] = "G/G/C/C/T/A";
          gts[24] = "G/G/G/G/C/T/A";
          gts[25] = "G/G/G/C/C/T/A";
          gts[26] = "G/G/C/C/T/T/A";
      }
    } else if(nucs[1] == 3){
      /* Bi-allelic */
      gts[1]  = "G/T";
      gts[2]  = "G/G/T";
      gts[4]  = "G/G/G/T";
      gts[7]  = "G/G/G/G/T";
      gts[8]  = "G/G/G/T/T";
      gts[12] = "G/G/G/G/G/T";
      gts[17] = "G/G/G/G/G/G/T";
      gts[18] = "G/G/G/G/G/T/T";
      gts[19] = "G/G/G/G/T/T/T";
      if(nucs[2] == 0){
        /* Tri-allelic */
        gts[3]  = "G/T/A";
        gts[5]  = "G/G/T/A";
        gts[9]  = "G/G/G/T/A";
        gts[10] = "G/G/T/T/A";
        gts[13] = "G/G/G/G/T/A";
        gts[14] = "G/G/G/T/T/A";
        gts[20] = "G/G/G/G/G/T/A";
        gts[21] = "G/G/G/G/T/T/A";
        gts[22] = "G/G/G/T/T/T/A";
        gts[23] = "G/G/G/T/T/A/A";
          /* Tetra-allelic */
          gts[6]  = "G/T/A/C";
          gts[11] = "G/G/T/A/C";
          gts[15] = "G/G/G/T/A/C";
          gts[16] = "G/G/T/T/A/C";
          gts[24] = "G/G/G/G/T/A/C";
          gts[25] = "G/G/G/T/T/A/C";
          gts[26] = "G/G/T/T/A/A/C";
      } else if(nucs[2] == 1){
        /* Tri-allelic */
        gts[3]  = "G/T/C";
        gts[5]  = "G/G/T/C";
        gts[9]  = "G/G/G/T/C";
        gts[10] = "G/G/T/T/C";
        gts[13] = "G/G/G/G/T/C";
        gts[14] = "G/G/G/T/T/C";
        gts[20] = "G/G/G/G/G/T/C";
        gts[21] = "G/G/G/G/T/T/C";
        gts[22] = "G/G/G/T/T/T/C";
        gts[23] = "G/G/G/T/T/C/C";
          /* Tetra-allelic */
          gts[6]  = "G/T/C/A";
          gts[11] = "G/G/T/C/A";
          gts[15] = "G/G/G/T/C/A";
          gts[16] = "G/G/T/T/C/A";
          gts[24] = "G/G/G/G/T/C/A";
          gts[25] = "G/G/G/T/T/C/A";
          gts[26] = "G/G/T/T/C/C/A";
      }
    }
  } else if (nucs[0] == 3) {
    gts[0] = "T/T";
    if(nucs[1] == 0){
      /* Bi-allelic */
      gts[1]  = "T/A";
      gts[2]  = "T/T/A";
      gts[4]  = "T/T/T/A";
      gts[7]  = "T/T/T/T/A";
      gts[8]  = "T/T/T/A/A";
      gts[12] = "T/T/T/T/T/A";
      gts[17] = "T/T/T/T/T/T/A";
      gts[18] = "T/T/T/T/T/A/A";
      gts[19] = "T/T/T/T/A/A/A";
      if(nucs[2] == 1){
        /* Tri-allelic */
        gts[3]  = "T/A/C";
        gts[5]  = "T/T/A/C";
        gts[9]  = "T/T/T/A/C";
        gts[10] = "T/T/A/A/C";
        gts[13] = "T/T/T/T/A/C";
        gts[14] = "T/T/T/A/A/C";
        gts[20] = "T/T/T/T/T/A/C";
        gts[21] = "T/T/T/T/A/A/C";
        gts[22] = "T/T/T/A/A/A/C";
        gts[23] = "T/T/T/A/A/C/C";
          /* Tetra-allelic */
          gts[6]  = "T/A/C/G";
          gts[11] = "T/T/A/C/G";
          gts[15] = "T/T/T/A/C/G";
          gts[16] = "T/T/A/A/C/G";
          gts[24] = "T/T/T/T/A/C/G";
          gts[25] = "T/T/T/A/A/C/G";
          gts[26] = "T/T/A/A/C/C/G";
      } else if(nucs[2] == 2){
        /* Tri-allelic */
        gts[3]  = "T/A/G";
        gts[5]  = "T/T/A/G";
        gts[9]  = "T/T/T/A/G";
        gts[10] = "T/T/A/A/G";
        gts[13] = "T/T/T/T/A/G";
        gts[14] = "T/T/T/A/A/G";
        gts[20] = "T/T/T/T/T/A/G";
        gts[21] = "T/T/T/T/A/A/G";
        gts[22] = "T/T/T/A/A/A/G";
        gts[23] = "T/T/T/A/A/G/G";
          /* Tetra-allelic */
          gts[6]  = "T/A/G/C";
          gts[11] = "T/T/A/G/C";
          gts[15] = "T/T/T/A/G/C";
          gts[16] = "T/T/A/A/G/C";
          gts[24] = "T/T/T/T/A/G/C";
          gts[25] = "T/T/T/A/A/G/C";
          gts[26] = "T/T/A/A/G/G/C";
      }
    } else if(nucs[1] == 1){
      /* Bi-allelic */
      gts[1]  = "T/C";
      gts[2]  = "T/T/C";
      gts[4]  = "T/T/T/C";
      gts[7]  = "T/T/T/T/C";
      gts[8]  = "T/T/T/C/C";
      gts[12] = "T/T/T/T/T/C";
      gts[17] = "T/T/T/T/T/T/C";
      gts[18] = "T/T/T/T/T/C/C";
      gts[19] = "T/T/T/T/C/C/C";
      if(nucs[2] == 0){
        /* Tri-allelic */
        gts[3]  = "T/C/A";
        gts[5]  = "T/T/C/A";
        gts[9]  = "T/T/T/C/A";
        gts[10] = "T/T/C/C/A";
        gts[13] = "T/T/T/T/C/A";
        gts[14] = "T/T/T/C/C/A";
        gts[20] = "T/T/T/T/T/C/A";
        gts[21] = "T/T/T/T/C/C/A";
        gts[22] = "T/T/T/C/C/C/A";
        gts[23] = "T/T/T/C/C/A/A";
          /* Tetra-allelic */
          gts[6]  = "T/C/A/G";
          gts[11] = "T/T/C/A/G";
          gts[15] = "T/T/T/C/A/G";
          gts[16] = "T/T/C/C/A/G";
          gts[24] = "T/T/T/T/C/A/G";
          gts[25] = "T/T/T/C/C/A/G";
          gts[26] = "T/T/C/C/A/A/G";
      } else if(nucs[2] == 2){
        /* Tri-allelic */
        gts[3]  = "T/C/G";
        gts[3]  = "T/C/G";
        gts[5]  = "T/T/C/G";
        gts[9]  = "T/T/T/C/G";
        gts[10] = "T/T/C/C/G";
        gts[13] = "T/T/T/T/C/G";
        gts[14] = "T/T/T/C/C/G";
        gts[20] = "T/T/T/T/T/C/G";
        gts[21] = "T/T/T/T/C/C/G";
        gts[22] = "T/T/T/C/C/C/G";
        gts[23] = "T/T/T/C/C/G/G";
          /* Tetra-allelic */
          gts[6]  = "T/C/G/A";
          gts[11] = "T/T/C/G/A";
          gts[15] = "T/T/T/C/G/A";
          gts[16] = "T/T/C/C/G/A";
          gts[24] = "T/T/T/T/C/G/A";
          gts[25] = "T/T/T/C/C/G/A";
          gts[26] = "T/T/C/C/G/G/A";
      }
    } else if(nucs[1] == 2){
      /* Bi-allelic */
      gts[1]  = "T/G";
      gts[2]  = "T/T/G";
      gts[4]  = "T/T/T/G";
      gts[7]  = "T/T/T/T/G";
      gts[8]  = "T/T/T/G/G";
      gts[12] = "T/T/T/T/T/G";
      gts[17] = "T/T/T/T/T/T/G";
      gts[18] = "T/T/T/T/T/G/G";
      gts[19] = "T/T/T/T/G/G/G";
      if(nucs[2] == 0){
        /* Tri-allelic */
        gts[3]  = "T/G/A";
        gts[5]  = "T/T/G/A";
        gts[9]  = "T/T/T/G/A";
        gts[10] = "T/T/G/G/A";
        gts[13] = "T/T/T/T/G/A";
        gts[14] = "T/T/T/G/G/A";
        gts[20] = "T/T/T/T/T/G/A";
        gts[21] = "T/T/T/T/G/G/A";
        gts[22] = "T/T/T/G/G/G/A";
        gts[23] = "T/T/T/G/G/A/A";
          /* Tetra-allelic */
          gts[6]  = "T/G/A/C";
          gts[11] = "T/T/G/A/C";
          gts[15] = "T/T/T/G/A/C";
          gts[16] = "T/T/G/G/A/C";
          gts[24] = "T/T/T/T/G/A/C";
          gts[25] = "T/T/T/G/G/A/C";
          gts[26] = "T/T/G/G/A/A/C";
      } else if(nucs[2] == 1){
        /* Tri-allelic */
        gts[3]  = "T/G/C";
        gts[3]  = "T/G/C";
        gts[5]  = "T/T/G/C";
        gts[9]  = "T/T/T/G/C";
        gts[10] = "T/T/G/G/C";
        gts[13] = "T/T/T/T/G/C";
        gts[14] = "T/T/T/G/G/C";
        gts[20] = "T/T/T/T/T/G/C";
        gts[21] = "T/T/T/T/G/G/C";
        gts[22] = "T/T/T/G/G/G/C";
        gts[23] = "T/T/T/G/G/C/C";
          /* Tetra-allelic */
          gts[6]  = "T/G/C/A";
          gts[11] = "T/T/G/C/A";
          gts[15] = "T/T/T/G/C/A";
          gts[16] = "T/T/G/G/C/A";
          gts[24] = "T/T/T/T/G/C/A";
          gts[25] = "T/T/T/G/G/C/A";
          gts[26] = "T/T/G/G/C/C/A";
      }
    }
  }

/*
  cout << "*** Debug int_to_nucs ***";
  cout << "\n";
  cout << nucs[0] << "," << nucs[1] << "," << nucs[2]<< "," << nucs[3];
  cout << "\n";
  cout << gts[0] << "," << gts[1];
  cout << "\n";
*/

//  cout << "\n";
}


/* Create type and function to help sort nucleotides. */
typedef std::pair<int,int> mypair;
bool comparator ( const mypair& l, const mypair& r)
    { return l.first > r.first; }


//void counts_2_plh(int mlhs[51], int nuc_cnts[8], float error, int min_cnt, int debug=0){
//void counts_2_plh(int mlhs[27], int nuc_cnts[8], float error, int debug=0){
void counts_2_plh(int& mlhs, string& gt, int nuc_cnts[8], float error, 
                  int min_count, int max_count, int debug=0){

  /* Add forward and reverse counts. */
  int nuc_cnt[4];
  nuc_cnt[0] = nuc_cnts[0] + nuc_cnts[1];
  nuc_cnt[1] = nuc_cnts[2] + nuc_cnts[3];
  nuc_cnt[2] = nuc_cnts[4] + nuc_cnts[5];
  nuc_cnt[3] = nuc_cnts[6] + nuc_cnts[7];

  /* Apply minimum threshold */
  for(int i=0; i < 4; i++){
    if(nuc_cnt[i] < min_count){nuc_cnt[i] = 0;}
  }

  /* Sort nucleotide counts. */
  vector<pair<int,int>> moves = {
    {nuc_cnt[0], 0},
    {nuc_cnt[1], 1},
    {nuc_cnt[2], 2},
    {nuc_cnt[3], 3}
  };

  std::sort(moves.begin(), moves.end(), comparator);

  /* Possible ways to obtain observed counts */
  double poss_counts = possible_counts(nuc_cnt);


  /* --*-- --*-- --*-- */
  /*       Models.     */
  /* --*-- --*-- --*-- */

  /* Likelihoods. */
  float mls [27];
  for(int j=0; j<27; j++){mls[j]=0;}

  /* Homozygote. */
  mls[0] = poss_counts * 
             pow(1-(3*error)/4, moves[0].first) * 
             pow(error/4, moves[1].first+moves[2].first+moves[3].first);

  /* Bi-allelic heterozygote. */
  mls[1] = poss_counts *
             pow(0.5-error/4, moves[0].first+moves[1].first) *
             pow(error/4, moves[2].first+moves[3].first);

  /* Bi-allelic triploid. */
  float prop3 = 0.3333333;
  float prop6 = 0.6666667;
  mls[2] = poss_counts * 
             pow(0.6666667-error/4, moves[0].first) * 
             pow(0.3333333-error/4, moves[1].first) *
             pow(error/4, moves[2].first+moves[3].first);

  /* Tri-allelic triploid loci. */
  mls[3] = poss_counts *
             pow(0.3333333-error/4, moves[0].first+moves[1].first+moves[2].first) *
             pow(error/4, moves[3].first);

  /* Bi-allelic tetraploid loci */
  mls[4] = poss_counts *
             pow(0.75-error/4, moves[0].first) *
             pow(0.25-error/4, moves[1].first) *
             pow(error/4, moves[2].first+moves[3].first);

  /* Tri-allelic tetraploid loci */
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
              pow(0.4285714-error/4, moves[0].first) *
              pow(0.4285714-error/4, moves[1].first) *
              pow(0.1428571-error/4, moves[1].first) *
              pow(error/4, moves[3].first);
  mls[23] = poss_counts *
              pow(0.4285714-error/4, moves[0].first) *
              pow(0.2857143-error/4, moves[1].first) *
              pow(0.2857143-error/4, moves[1].first) *
              pow(error/4, moves[3].first);

  /* Tetra-allelic septaploid */
  mls[24] = poss_counts *
              pow(0.5714286-error/4, moves[0].first) *
              pow(0.1428571-error/4, moves[1].first) *
              pow(0.1428571-error/4, moves[2].first) *
              pow(0.1428571-error/4, moves[3].first);
  mls[25] = poss_counts *
              pow(0.4285714-error/4, moves[0].first) *
              pow(0.2857143-error/4, moves[1].first) *
              pow(0.1428571-error/4, moves[2].first) *
              pow(0.1428571-error/4, moves[3].first);
  mls[26] = poss_counts *
              pow(0.2857143-error/4, moves[0].first) *
              pow(0.2857143-error/4, moves[1].first) *
              pow(0.2857143-error/4, moves[2].first) *
              pow(0.1428571-error/4, moves[3].first);


  /* Convert ints to nucleotides */
  string gts[27];
  int_to_nucs(gts, moves);


  /* Identify genotype with the maximum likelihood. */
  int maxl = 0;
  float suml = mls[0];
  for (int j = 1; j < 27; j++){
    if(mls[j] > mls[maxl]){
      maxl = j;
    }
    suml = suml + mls[j];
  }

  /* Phred scale the likelihoods. */

  for (int j = 0; j < 27; j++){
    if (mls[j] == 0){
      mls[j] = 9999;
    } else {
      mls[j] = trunc(-10 * log10(mls[j]));
      if(mls[j] == -0){mls[j] = 0;}
    }
  }


  /* Debuging */
//  debug = 1;
  if(debug == 1){
    cout << "\n";
    cout << "*** Debug counts_2_phredlh ***\n";
    cout << "Counts: " << nuc_cnt[0] << "," << nuc_cnt[1] << "," << nuc_cnt[2]<< "," << nuc_cnt[3];
    cout << "\n";
    cout << "possible counts=" << poss_counts;
    cout << "\n";

    cout << 0 << "\t" << mls[0] << "\t" << gts[0];
    cout << "\n";
    for (int j = 1; j < 27; j++){
      cout << j << "\t" << mls[j] << "\t" << gts[j];
      cout << "\n";
    }
    cout << "\n";

    cout << "Most likely genotype: " << gts[maxl];
    cout << "\n";
    if(suml >= 27){
      cout << "All genotypes equally likely: sum = " << suml;
      cout << "\n";
      cout << "Set genotype to null: ./.";
      cout << "\n";
    }
  }

  /* Transfer values to parent. */
  if(suml >= 27){
    /* All likelihoods are equal. */
    mlhs = int(mls[maxl]);
    gt = "./.";
  } else if(maxl >= 17){
    /* Called genotype is septaploid */
    mlhs = int(mls[maxl]);
    gt = "./.";
  } else if(nuc_cnt[0] + nuc_cnt[1] + nuc_cnt[2] + nuc_cnt[3] >= max_count){
    /* Called genotype has coverage above max threshold */
    gt = "./.";
  } else {
    mlhs = int(mls[maxl]);
    gt = gts[maxl];
  }
}


/* Parse each site to samples */
//void mult_pl(int pls[][27], int nsamp, float err, int nuc_cnts[][8]){
void mult_pl(int pls[], string gts[], int nsamp, float err, int nuc_cnts[][8], vector <int> min_count, vector <int> max_count){
  for(int i=0; i<nsamp; i++){
    counts_2_plh(pls[i], gts[i], nuc_cnts[i], err, min_count[i], max_count[i]);
  }



}



/* Determine if a site is polymorphic */
int site_polymorphic(int nsamp, string gts[]){
//  cout << gts[0];
  int comp = 0;
  if(gts[0] == gts[1]){comp = 1;}
  for(int i=1; i<nsamp; i++){
//    cout << "," << gts[i];
    if(gts[0] == gts[i]){comp = comp + 1;}
//    comp = comp + gts[0] == gts[i];
  }
//  cout << "\n";
//  cout << "Comparison = " << comp;
//  cout << "; Number of samples = " << nsamp;
//  cout << "\n";
//  cout << "\n";
  if(comp == nsamp){
    return(0);
  } else {
    return(1);
  }
}


//int diff_ref(int nsamp, string gts[]){
int diff_ref(string ref, int nsamp, string gts[]){
  vector <string> haps;
  int gtsum = 0;

  for(int i=0; i<nsamp; i++){
    if(gts[i] != "./."){
      split( haps, gts[i], boost::algorithm::is_any_of( "/" ) );
      for(int j=1; j<haps.size(); j++){
        if(haps[j] != ref){gtsum++;}
      }
    }
  }

  if(gtsum > 0){
    return(1);
  } else {
    return (0);
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
  cerr << "  -s file with lower and upper count thresholds for calling genotypes\n";
  cerr << "     on row one and two, same number of columns as samples.\n";
  cerr << "\n";
  cerr << "Expects piped output from SAMTools::mpileup where\n";
  cerr << "each sample consists of three columns starting at\n";
  cerr << "the fourth column.\n";
  cerr << "\n";
}


void print_header(float error, string sfile, string tfile){
  vector <string> fields;
  cout << "##fileformat=VCFv4.2\n";
  cout << "##source=gtV0.0.0\n";

  /* Threshold file */
  if(tfile != "NA"){
    string line;
    ifstream myfile (tfile);
    if (myfile.is_open())
    {
      cout << "##FILTER=<ID=NA,Description=\"Lower count threshold for calling genotypes = ";
      getline (myfile,line);
      split( fields, line, boost::algorithm::is_any_of( "\t " ) );
      cout << fields[0];
      for(int i=1; i<fields.size(); i++){
        cout << "," << fields[i];
      }
      cout << "\">";
      cout << "\n";
      cout << "##FILTER=<ID=NA,Description=\"Upper count threshold for calling genotypes = ";
      getline (myfile,line);
      split( fields, line, boost::algorithm::is_any_of( "\t " ) );
      cout << fields[0];
      for(int i=1; i<fields.size(); i++){
        cout << "," << fields[i];
      }
      cout << "\">";
      cout << "\n";
    myfile.close();
    }
    else cout << "Unable to open file"; 
  }
//
//  cout << "\n";

  cout << "##FILTER=<ID=NA,Description=\"Error rate of " << error << " used for likelihood calculation\">\n";
  cout << "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Read depth\">\n";
  cout << "##FORMAT=<ID=CT,Number=8,Type=Integer,Description=\"Count of each nucleotide (A,a,C,c,G,g,T,t)\">\n";
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


void print_fix(vector <string> fields){
  cout << fields[0] << "\t" << fields[1] << "\t" << ".";
  cout << "\t" << fields[2] << "\t" << ".";
  cout << "\t" << "." << "\t" << ".";
  cout << "\t" << ".";
}


//void print_locus(vector <string> fields, int counts, int phred, int nsamp, int rds[], int nuc_cnts[][8], int pls[][51], string gts[]){
void print_locus(vector <string> fields, int counts, int phred, int nsamp, int rds[], int nuc_cnts[][8], int pls[], string gts[]){
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
      cout << ":" << pls[i];
//      for(int j=1; j<25; j++){cout << "," << pls[i][j];}
    }
    // Genotype.
    cout << ":" << gts[i];
  }
  cout << "\n";
}




/* ----- ----- ***** ----- ----- */
/*             Main              */
/* ----- ----- ***** ----- ----- */


int main(int argc, char **argv) {
  string lineInput;
  vector <string> fields;
  vector <int> min_count;
  vector <int> max_count;

  int opt; // options
  int header = 0; // print header/meta data
  int counts = 0; // print counts 
  int phred = 0; // Scale likelihoods as phred values
  float error = 0.000000001; // Can not be zero!!!
  string sfile = "NA";
  string tfile = "NA";

  /* Parse command line options. */
//  while ((opt = getopt(argc, argv, "ce:hmps:t:")) != -1) {
  while ((opt = getopt(argc, argv, "ce:hmps:t:")) != -1) {
    switch (opt) {
      case 'c': // print counts
        counts = 1;
        break;
      case 'e': // permissible error
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
      case 't': // file containing upper and lower count thresholds
        tfile = optarg;
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s [-ce:hmps]\n", argv[0]);
        print_usage();
        exit(EXIT_FAILURE);
    }
  }

  /* Print header */
  if(header == 1){print_header(error, sfile, tfile);}

  /* Read in thresholds */
  if(tfile != "NA"){
    string line;
    ifstream myfile (tfile);
    if (myfile.is_open()){

      /* Minimum threshold */
      getline (myfile,line);
      split( fields, line, boost::algorithm::is_any_of( "\t " ) );
      for(int i=0; i<fields.size(); i++){
        min_count.push_back(stoi(fields[i]));
      }

      /* Maximum threshold */
      getline (myfile,line);
      split( fields, line, boost::algorithm::is_any_of( "\t " ) );
      for(int i=0; i<fields.size(); i++){
        max_count.push_back(stoi(fields[i]));
      }

    myfile.close();
    }
    else cout << "Unable to open file"; 
  }

/*
  for(int i=0; i < min_count.size(); i++){
    cout << min_count[i] << "\t" << max_count[i] << "\n";
  }
*/


//  }

  /* Parse each line (site) by sample (column). */
  while (getline(cin,lineInput)) {
    split( fields, lineInput, boost::algorithm::is_any_of( "\t " ) );

    /* Declare variables */
    int nsamp = (fields.size()-3)/3;  // Determine the number of samples.
//    int nsamp = (fields.size()-2)/2;  // Determine the number of samples.
    int nuc_cnts [nsamp][8]; // A,a,C,c,G,g,T,t.
    int rds [nsamp]; // Read depth.

    string gts [nsamp]; // Genotypes.
    int pls [nsamp];  // Phred scaled likelihoods. 
//    int pls [nsamp][27];  // Phred scaled likelihoods. 

    /* Initialize variables. */
    for(int i=0; i<nsamp; i++){
      rds[i] = 0;
      gts[i] = "./.";
      pls[i] = 9999;
      if(tfile == "NA"){
        min_count.push_back(0);
        max_count.push_back(0);
      }
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
//    string gts[27];
//    mult_pl(pls, nsamp, error, nuc_cnts);
    mult_pl(pls, gts, nsamp, error, nuc_cnts, min_count, max_count);


/*
    cout << gts[0];
    for(int j=1; j<nsamp; j++){
      cout << "," << gts[j];
    }
    cout << "\n";

    cout << pls[0];
    for(int j=1; j<nsamp; j++){
      cout << "," << pls[j];
    }
    cout << "\n";
*/


    /* Print locus. */
//    int unique_gt = site_polymorphic(nsamp, gts);
    int unique_gt = diff_ref(fields[2], nsamp, gts);
//    int unique_gt = diff_ref(nsamp, gts);

    if(unique_gt == 1){
      print_locus(fields, counts, phred, nsamp, rds, nuc_cnts, pls, gts);
      /* Debug */
//      debug1(fields, counts, phred, nsamp, rds, nuc_cnts, pls, gts, error, min_cnt);
    }
  }
  return 0;
}



/* ----- ----- ***** ----- ----- */
// EOF.
