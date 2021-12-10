 # build an executable
all: Case_1_P_negative.cpp Case_1_P.cpp Case_2.1_P_Prime.cpp Case_2.2_P_Double_Prime.cpp
	g++ -O3 -Wall -o Case_1_P_negative.out Case_1_P_negative.cpp -lm
	g++ -O3 -Wall -o Case_1_P.out Case_1_P.cpp -lm
	g++ -O3 -Wall -o Case_2.1_P_Prime.out Case_2.1_P_Prime.cpp -lm
	g++ -O3 -Wall -o Case_2.2_P_Double_Prime.out Case_2.2_P_Double_Prime.cpp -lm


