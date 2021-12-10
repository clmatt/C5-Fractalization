#include <iostream>
#include <cmath>

using namespace std;

const int minNum = 9;
const int maxNum = 1000;
const double eps = 0.0000001; //Tolerance

//This program check the integer program P to see if there are any negative values
//It does this by simply checking all possible values of inputs
//In the case where there is a negative value (also listed in the appendix of the paper) it is outputted, otherwise nothing is output if the inputs have negative size
//The maxNum can increase past 1000, but the program slows significantly as n grows, and 1000 suffices for the next results 


int main() {
	long long int i1, i2, i3, i4, i, j, k, n;
	double f, d, val, temp;

	int x[6]; //x[i] = x_i n in P
	double C5Vec[maxNum+1], C5sVec[maxNum+1], bin[maxNum+1];

	double C5(long long int);
	long long int bin5(long long int);
	double C5s(long long int);

	//Initialize Variables
	for (n = 0; n <= 5; n++) {
		x[n] = 0;
	}

	//Precalculate values
	for (i = 0; i <= maxNum; ++i) {
		C5Vec[i] = C5(i);
		C5sVec[i] = (double)C5s(i);
		bin[i] = bin5(i);
	}

	for (n = minNum; n <= maxNum; n++) {
		cout << "Checking n = " << n << endl;

		//Run through possibilities of x_i's
		for (i1 = 0; i1 <= n; i1++) {
			for (i2 = i1; i2 <= n - i1; i2++) {
				for (i3 = i2; i3 <= n - i1 - i2; i3++) {
					for (i4 = i3; 2 * i4 <= n - i1 - i2 - i3; i4++) { //This 2 guarantees x[5] big enough so ordering is correct
						x[1] = i1;
						x[2] = i2;
						x[3] = i3;
						x[4] = i4;
						x[5] = n - i1 - i2 - i3 - i4;

						//Sum from bound with f
						temp = 0.0;
						for (i = 1; i <= 5; i++) {
							for (j = i + 1; j <= 5; j++) {
								temp = temp + (x[i] / ((double)n)) * (x[j] / ((double)n));
							}
						}

						f = (2 * n) / ((double)n - 1) * (temp - 2.0 * (-0.175431374077117 + 8.75407592662244 * C5sVec[n]) / (21.0 * C5sVec[n]));

						//Make sure we don't have a negative f
						if (f > -eps) {
							d = (f * n * (n - 1) - 2) / ((double)(2 * n));
							k = floor((double)n / (double)5);
							j = n % 5;

							//Objective function
							val = f * n * (n - 1) / ((double)2) * (
								x[3] * x[4] * x[5]
								- 3 * n * d * x[3] * x[4] / ((double)8)
								- f * x[3] * n * n / ((double)8)
								- (f * n * n - (f + d) * n - 1) / ((double)8) * (2 * x[1] + 2 * x[2] + x[3] + x[4] + x[5])
								- 9 * x[1] * x[1] * (d * n + 2) / ((double)32)
								)
								+ pow(k, 5 - j) * pow(k + 1, j)
								+ (5 - j) * C5Vec[k] * bin[k]
								+ j * C5Vec[k + 1] * bin[k + 1]
								- x[1] * x[2] * x[3] * x[4] * x[5]
								- C5Vec[x[1]] * bin[x[1]]
								- C5Vec[x[2]] * bin[x[2]]
								- C5Vec[x[3]] * bin[x[3]]
								- C5Vec[x[4]] * bin[x[4]]
								- C5Vec[x[5]] * bin[x[5]];

							if (val < eps) {
								cout << "Found a negative value, with part sizes: " << endl;
								cout << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << endl;
							}
						}
					}
				}
			}
		}
	}


	return 0;
}

//Compute n choose 5
long long int bin5(long long int n) {

	if (n < 5) return 0;

	else return n * (n - 1) * (n - 2) * (n - 3) * (n - 4) * (n - 5) / 120;
}

//Computes C5
double C5(long long int n) {
	long long int nf, nc;
	int a;

	if (n < 0) return -1;

	else if (n < 5) return 0;

	else if (n == 5) return 1;

	else {
		nf = floor(n / ((double)5));
		nc = ceil(n / ((double)5));
		a = n % 5;

		if (nf < 5 && nc < 5) {
			return  120 * (pow(nf, 5 - a) / ((double)n) * pow(nc, a) / ((double)(n - 1))) / ((double)(n - 2) * (n - 3) * (n - 4));
		}
		else if (nf < 5 && nc == 5) {
			return 120 * (pow(nf, 5 - a) / ((double)n) * pow(nc, a) / ((double)(n - 1))) / ((double)(n - 2) * (n - 3) * (n - 4))
				+ ((a)*C5(nc)) * (nc / ((double)n)) * (nc - 1) / ((double)(n - 1)) * (nc - 2) / ((double)(n - 2)) * (nc - 3) / ((double)(n - 3)) * (nc - 4) / ((double)(n - 4));
		}
		else {
			return 120 * (pow(nf, 5 - a) / ((double)n) * pow(nc, a) / ((double)(n - 1))) / ((double)(n - 2) * (n - 3) * (n - 4))
				+ ((5 - a) * C5(nf)) * (nf / ((double)n)) * (nf - 1) / ((double)(n - 1)) * (nf - 2) / ((double)(n - 2)) * (nf - 3) / ((double)(n - 3)) * (nf - 4) / ((double)(n - 4))
				+ ((a)*C5(nc)) * (nc / ((double)n)) * (nc - 1) / ((double)(n - 1)) * (nc - 2) / ((double)(n - 2)) * (nc - 3) / ((double)(n - 3)) * (nc - 4) / ((double)(n - 4));
		}
	}
}

//Computes C5(n^*)
double C5s(long long int n) {
	if (n < 5) return 0;

	else return (n / ((double)26) + n * (n - 1) * (n - 2) * (n - 3) * (n - 4) * C5(n)) / ((double)n * n * n * n * n);
}


