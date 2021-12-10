#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

//This program checks possible counterexamples to the first nonlinear program P
//The program takes possible values of part sizes and outputs whether these part sizes could possibly give more C5s then the balanced iterated blowup of a C5
//The program bounds the possible number of funky edges using Lemma 2.2

//The only things that should need to be changed is this n and what the x_i's are (a few lines below)
//n is the sum of the x_i's
const int n = 9; 

int main() {
	int i, j, k, f, numPairs, counter, C5count, maxC5;
	double temp, fdble;
	bool repeat1, repeat2, C5check, myVal;

	long long int numIt;

	int x[6], funky[9], fiveSet[5], a[3],b[5]; //Index is max funky + 1 for index, funky, possible pairs is large enough
	bool insideEdges[n + 1][n + 1], edges[n + 1][n + 1];

	vector < vector<int> > possiblePairs(1000, vector<int>(3));

	double C5(int);
	double C5s(int);
	void nextSubset(int[], int, int, bool&);

	//Values taken from Medium Code
	//Input these manually 
	//Order of these matters!
	x[0] = 0; //Always zero helps with off by one errors
	x[1] = 1;
	x[2] = 2;
	x[3] = 2;
	x[4] = 2;
	x[5] = 2;
	
	cout << "Checking a graph with part sizes: " << x[1] << ", " << x[2] << ", " << x[3] << ", " << x[4] << ", " << x[5] << "." << endl << endl;

	if (x[1] + x[2] + x[3] + x[4] + x[5] != n) {
		cout << "FIX N!!!!" << endl;
		return -2;
	}

	numIt = 0;

	for (i = 0; i <= n; ++i) {
		for (j = 0; j <= n; ++j) {
			insideEdges[i][j] = false;
			edges[i][j] = false;
		}
	}

	//Determine which edges are inside parts
	counter = 1;
	for (i = 1; i <= 5; ++i) {
		for (j = counter; j < counter + x[i] - 1; ++j) {
			for (k = j + 1; k <= counter + x[i] - 1; ++k) {
				insideEdges[j][k] = true;
				insideEdges[k][j] = true;
			}
		}
		counter = counter + x[i];
	}

	//Calculate max number of funky edges
	temp = 0.0;
	for (i = 1; i < 5; i++) {
		for (j = i + 1; j <= 5; j++) {
			temp = temp + (x[i] / ((double)n)) * (x[j] / ((double)n));
		}
	}
	fdble = ((2 * n) / ((double)n - 1) * (temp - 2.0 * (-0.175431374077117 + 8.75407592662244 * C5s(n)) / (21.0 * C5s(n)))); //Lemma 2.2
	f = floor(n * (n - 1) * fdble / ((double)2));

	cout << "The maximum number of funky edges is: " << f << endl << endl;


	//Create edge set without any funky edges
	numPairs = 0;
	for (i = 1; i <= x[1]; i++) {
		for (j = 1; j <= x[2]; j++) {
			possiblePairs[numPairs][1] = i;
			possiblePairs[numPairs][2] = j+x[1];
			edges[i][j+x[1]] = true;
			numPairs = numPairs + 1;
		}
	}
	for (i = 1; i <= x[2]; i++) {
		for (j = 1; j <= x[3]; j++) {
			possiblePairs[numPairs][1] = i + x[1];
			possiblePairs[numPairs][2] = j + x[1] + x[2];
			edges[i+x[1]][j+x[1]+x[2]] = true;
			numPairs = numPairs + 1;
		}
	}
	for (i = 1; i <= x[3]; i++) {
		for (j = 1; j <= x[4]; j++) {
			possiblePairs[numPairs][1] = i + x[1] + x[2];
			possiblePairs[numPairs][2] = j + x[1] + x[2] + x[3];
			edges[i+x[1]+x[2]][j+x[1]+x[2]+x[3]] = true;
			numPairs = numPairs + 1;
		}
	}
	for (i = 1; i <= x[4]; i++) {
		for (j = 1; j <= x[5]; j++) {
			possiblePairs[numPairs][1] = i + x[1] + x[2] + x[3];
			possiblePairs[numPairs][2] = j + x[1] + x[2] + x[3] + x[4];
			edges[i + x[1] + x[2] + x[3]][j + x[1] + x[2] + x[3] + x[4]] = true;
			numPairs = numPairs + 1;
		}
	}
	for (i = 1; i <= x[1]; i++) {
		for (j = 1; j <= x[5]; j++) {
			possiblePairs[numPairs][1] = i;
			possiblePairs[numPairs][2] = j + x[1] + x[2] + x[3] + x[4];
			edges[i][j + x[1] + x[2] + x[3] + x[4]] = true;
			numPairs = numPairs + 1;
		}
	}


	maxC5 = 0;

	//Iterate through all possible numbers of funky edges
	for (i = 1; i <= f; i++) {
		//Intialize indices of funky edges
		//Funky is indices of possiblePairs
		for (k = 0; k < i; k++) {
			funky[k] = k;
		}

		repeat1 = false; //Check whether we've gone through all possible subsets of funky edges

		//Just assuming that inside edges are there if we need them and not if we don't
		while (repeat1 == false) {
			//Create edge set (flip funky edges)
			for (k = 0; k < i; k++) {
				edges[possiblePairs[funky[k]][1]][possiblePairs[funky[k]][2]] = !edges[possiblePairs[funky[k]][1]][possiblePairs[funky[k]][2]];
			}

			//Make edges symmetric
			for (k = 1; k < n; k++) {
				for (j = k + 1; j <= n; j++) {
					edges[j][k] = edges[k][j];
				}
			}

			//Intialize 5 set to count # of C5s
			//fiveSet is which vertices are considered for a C5
			for (k = 0; k < 5; k++) {
				fiveSet[k] = k + 1;
			}

			C5count = 0;
			repeat2 = false;

			//Count nuber of C5s where 4->4 3->2
			while (repeat2 == false) {
				//Need to check all possible permutations of the 5 points to check if any of them make a C5
				a[0] = 0;
				a[1] = 1;
				a[2] = 2;

				sort(a, a + 3);

				myVal = false;

				do {
					/*
					if (numIt % 10000000 == 0) {
						cout << "The number of iterations is: " << numIt << endl;
					}*/

					b[0] = a[0];
					b[1] = a[1];
					b[2] = a[2];
					b[3] = 3;
					b[4] = 4;

					//Hardcoding check if there is a c5
					C5check = true;

					if (insideEdges[fiveSet[b[1]]][fiveSet[b[2]]] == false && edges[fiveSet[b[1]]][fiveSet[b[2]]] == false) {
						C5check = false;
					}
					else if (insideEdges[fiveSet[b[2]]][fiveSet[b[3]]] == false && edges[fiveSet[b[2]]][fiveSet[b[3]]] == false) {
						C5check = false;
					}
					else if (insideEdges[fiveSet[b[3]]][fiveSet[b[4]]] == false && edges[fiveSet[b[3]]][fiveSet[b[4]]] == false) {
						C5check = false;
					}
					else if (insideEdges[fiveSet[b[0]]][fiveSet[b[4]]] == false && edges[fiveSet[b[0]]][fiveSet[b[4]]] == false) {
						C5check = false;
					}
					else if (insideEdges[fiveSet[b[0]]][fiveSet[b[1]]] == false && edges[fiveSet[b[0]]][fiveSet[b[1]]] == false) {
						C5check = false;
					}
					else if (insideEdges[fiveSet[b[1]]][fiveSet[b[3]]] == false && edges[fiveSet[b[1]]][fiveSet[b[3]]] == true) {
						C5check = false;
					}
					else if (insideEdges[fiveSet[b[1]]][fiveSet[b[4]]] == false && edges[fiveSet[b[1]]][fiveSet[b[4]]] == true) {
						C5check = false;
					}
					else if (insideEdges[fiveSet[b[2]]][fiveSet[b[4]]] == false && edges[fiveSet[b[2]]][fiveSet[b[4]]] == true) {
						C5check = false;
					}
					else if (insideEdges[fiveSet[b[0]]][fiveSet[b[2]]] == false && edges[fiveSet[b[0]]][fiveSet[b[2]]] == true) {
						C5check = false;
					}
					else if (insideEdges[fiveSet[b[0]]][fiveSet[b[3]]] == false && edges[fiveSet[b[0]]][fiveSet[b[3]]] == true) {
						C5check = false;
					}

					numIt = numIt + 1;

					myVal = myVal || C5check;
				} while (next_permutation(a, a + 3) && !myVal);

				if (myVal == true) ++C5count;

				else {
					//Count nuber of C5s where 4->4 3->3
					//Need to check all possible permutations of the 5 points to check if any of them make a C5
					a[0] = 0;
					a[1] = 1;
					a[2] = 2;

					sort(a, a + 3);

					myVal = false;

					do {
						//If interested in run time, can uncheck this
						/*
						if (numIt % 10000000 == 0) {
							cout << "The number of iterations is: " << numIt << endl;
						}*/

						b[0] = a[0];
						b[1] = a[1];
						b[2] = 3;
						b[3] = a[2];
						b[4] = 4;

						//Hardcoding check if there is a c5
						C5check = true;

						if (insideEdges[fiveSet[b[1]]][fiveSet[b[2]]] == false && edges[fiveSet[b[1]]][fiveSet[b[2]]] == false) {
							C5check = false;
						}
						else if (insideEdges[fiveSet[b[2]]][fiveSet[b[3]]] == false && edges[fiveSet[b[2]]][fiveSet[b[3]]] == false) {
							C5check = false;
						}
						else if (insideEdges[fiveSet[b[3]]][fiveSet[b[4]]] == false && edges[fiveSet[b[3]]][fiveSet[b[4]]] == false) {
							C5check = false;
						}
						else if (insideEdges[fiveSet[b[0]]][fiveSet[b[4]]] == false && edges[fiveSet[b[0]]][fiveSet[b[4]]] == false) {
							C5check = false;
						}
						else if (insideEdges[fiveSet[b[0]]][fiveSet[b[1]]] == false && edges[fiveSet[b[0]]][fiveSet[b[1]]] == false) {
							C5check = false;
						}
						else if (insideEdges[fiveSet[b[1]]][fiveSet[b[3]]] == false && edges[fiveSet[b[1]]][fiveSet[b[3]]] == true) {
							C5check = false;
						}
						else if (insideEdges[fiveSet[b[1]]][fiveSet[b[4]]] == false && edges[fiveSet[b[1]]][fiveSet[b[4]]] == true) {
							C5check = false;
						}
						else if (insideEdges[fiveSet[b[2]]][fiveSet[b[4]]] == false && edges[fiveSet[b[2]]][fiveSet[b[4]]] == true) {
							C5check = false;
						}
						else if (insideEdges[fiveSet[b[0]]][fiveSet[b[2]]] == false && edges[fiveSet[b[0]]][fiveSet[b[2]]] == true) {
							C5check = false;
						}
						else if (insideEdges[fiveSet[b[0]]][fiveSet[b[3]]] == false && edges[fiveSet[b[0]]][fiveSet[b[3]]] == true) {
							C5check = false;
						}

						numIt = numIt + 1;

						myVal = myVal || C5check;
					} while (next_permutation(a, a + 3) && !myVal);

					if (myVal == true) ++C5count;
				}

				nextSubset(fiveSet, n+1, 5, repeat2);
			}

			maxC5 = max(C5count, maxC5);

			//Flip Funky edges back
			for (k = 0; k < i; k++) {
				edges[possiblePairs[funky[k]][1]][possiblePairs[funky[k]][2]] = !edges[possiblePairs[funky[k]][1]][possiblePairs[funky[k]][2]];
			}

			//Make edges symmetric
			for (k = 1; k < n; k++) {
				for (j = k + 1; j <= n; j++) {
					edges[j][k] = edges[k][j];
				}
			}
			nextSubset(funky, numPairs, i, repeat1);
		}
		cout << "With at least " << i << " funky edges there are at most " << maxC5 << " C5s." << endl;
	}

	cout << endl << "The max number of C5s in the iterated blow up is: " << round(C5(n) * n * (n - 1) * (n - 2) * (n - 3) * (n - 4) / ((double)120)) << endl << endl;

	return 0;
}

double C5(int n) {
	int nf, nc;
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

double C5s(int n) {
	if (n < 5) return 0;

	else return (n / ((double)26) + n * (n - 1) * (n - 2) * (n - 3) * (n - 4) * C5(n)) / ((double)n * n * n * n * n);
}

void nextSubset(int subset[], int n, int k, bool& repeat) {
	int index,i;
	bool val;

	val = true;
	index = k-1;

	while (val == true) {

		if (index < 0) {
			repeat = true;
			val = false;
		}

		else if (subset[index] < n - k + index) {
			val = false;

			for (i = index; i < k; ++i) {
				if (i == index) {
					++subset[i];
				}
				else {
					subset[i] = subset[i - 1] + 1;
				}
			}
		}

		else if (subset[index] > n - k + index) {
			cout << "Something went wrong" << endl;
		}

		else {
			index = index - 1;
		}
	}
}

