
all:
	g++ -std=c++0x triplo_gt.cpp -lgmpxx -lgmp -o triplo_gt
	g++ -std=c++0x triplo_qc.cpp -o triplo_qc
	g++ -std=c++0x filter_by_cnt.cpp -o filter_by_cnt
