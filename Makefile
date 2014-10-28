
all:
	g++ -std=c++0x gt.cpp -lgmpxx -lgmp -o gt
	g++ -std=c++0x qc.cpp -o qc
	g++ -std=c++0x filter_by_cnt.cpp -o filter_by_cnt
