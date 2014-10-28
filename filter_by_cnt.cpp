// Compile with:
// g++ -std=c++0x filter_by_cnt.cpp -o filter_by_cnt
#include <iostream>
#include <string>

using namespace std;


int main(int argc, char *argv[]){
  int opt, thin=1000, rec=0;
  string lineInput;

  /* Parse command line options. */
  while ((opt = getopt(argc, argv, "n:")) != -1) {
    switch (opt) {
      case 'n':
        thin = atoi(optarg);
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s [-n]\n", argv[0]);
        exit(EXIT_FAILURE);
    }
  }

  /* Parse line by line. */
  while (getline(cin,lineInput)) {
/*
    if(lineInput.find("#") == 0){
      cout << lineInput;
    }
*/
    rec++;
    if(rec % thin == 0){
      cout << lineInput;
      cout << "\n";
    }
  }

  return 0;
}

