#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

//Shows that the program P' doesn't have any negative values
//Does an iterated grid search - take some portion of the search space and bound the value the objective function can have if it's positive move on to the next portion otherwise split into smaller pieces and repeat
//The program won't terminate if the solution of the program could be negative, otherwise it will termninate 

int main() {
	long long int i, j, k, n, counter, temp, bin,numIt;
	double C5, obj, f, d, eps;
	bool val,cons; //cons tells whether constraints are satisfied

	double xlTemp[6], xuTemp[6];

	std::vector < vector<double> > xu(1000, vector<double>(6)); //xu[.,0] = 0, xu[.,i] = x_i + epsilon
	std::vector < vector<double> > xl(1000, vector<double>(6));

	n = 1000;
	C5 = 0.0384609;
	val = true;
	counter = 0;
	numIt = 0;
	eps = 0.0000001; //Tolerance

	f = 2 * n / ((double)n - 1) * (10 / ((double)25) - 2.0 * (-0.175431374077116 + 8.75407592662244 * C5) / (21 * C5));
	d = 0.2; //Could slightly increase to get a better bound on d, but not necessary 

	for (i = 1; i < 6; ++i) {
		xl[counter][i] = 0.166;
		xu[counter][i] = 0.234;
	}

	while (val == true) {
		//Gives updates on running status
		//Uncomment if interested
		/*
		if (numIt % 100000 == 0) {
			cout << "The iteration number is: " << numIt << endl;
			cout << "Counter = " << counter << endl;
			cout << "The vector xl at the current step is: " << xl[counter][1] << " " << xl[counter][2] << " " << xl[counter][3] << " " << xl[counter][4] << " " << xl[counter][5] << endl;
			cout << "The vector xu at the current step is: " << xu[counter][1] << " " << xu[counter][2] << " " << xu[counter][3] << " " << xu[counter][4] << " " << xu[counter][5] << endl << endl;
		}*/
		

		numIt = numIt + 1;

		//Objective Values
		obj = xl[counter][3] * xl[counter][4] * xl[counter][5]
			- 3.0*d*xu[counter][3]*xu[counter][4]/((double)8)
			- f * xu[counter][3] / ((double)8)
			- f * xu[counter][3] / ((double)8)
			- f * (2.0 * xu[counter][1] + 2.0 * xu[counter][2] + xu[counter][3] + xu[counter][4] + xu[counter][5]) / ((double)8)
			- 9.0 * d * xu[counter][1] * xu[counter][1] / ((double)32)
			- 9.0 * xu[counter][1] * xu[counter][1] / ((double)16 * n);


		if (obj < eps) {
			for (i = 0; i < 6; ++i) {
				xlTemp[i] = xl[counter][i];
				xuTemp[i] = xu[counter][i];
			}

			for (i = 0; i < 16; ++i) {
				temp = i;
				cons = true;

				/*
				Bin is the binary representation of i
				Depending whether the ith digit is 0 or 1 we will either
				keep xl[i] the same and make xu[i] = (xu[i]+xl[i})/2 or vice versa
				*/
				for (j = 0; j < 4; ++j) {
					bin = temp % 2;
					temp = temp / 2;

					if (bin == 0) {
						xl[counter][j + 1] = xlTemp[j + 1];
						xu[counter][j + 1] = (xuTemp[j + 1] + xlTemp[j + 1]) / ((double)2) + eps; //
					}

					else {
						xl[counter][j + 1] = (xuTemp[j + 1] + xlTemp[j + 1]) / ((double)2) - eps; //
						xu[counter][j + 1] = xuTemp[j + 1];
					}
				}

				//Only 4 dimensions, can define x_5 based on other x_i's
				xl[counter][5] = std::max(0.0, 1 - xu[counter][1] - xu[counter][2] - xu[counter][3] - xu[counter][4]) - eps; //
				xu[counter][5] = 1 - xl[counter][1] - xl[counter][2] - xl[counter][3] - xl[counter][4] + eps; //

				//Check all constraints are satisfied
				//Proper ordering
				for (j = 1; j < 5; ++j) {
					for (k = j + 1; k < 6; ++k) {
						if (xu[counter][j] < xl[counter][k] - eps) cons = false; 
					}
				}

				//Check whether x5 has any feasible points
				if (xu[counter][5] < 0) cons = false;

				else if (xu[counter][1] * xu[counter][2] + xu[counter][1] * xu[counter][3] + xu[counter][1] * xu[counter][4] + xu[counter][1] * xu[counter][5]
					+ xu[counter][2] * xu[counter][3] + xu[counter][2] * xu[counter][4] + xu[counter][2] * xu[counter][5] + xu[counter][3] * xu[counter][4] + xu[counter][3] * xu[counter][5] + xu[counter][4] * xu[counter][5]
					< 2.0 * (-0.175431374077117 + 8.75407592662244 * C5) / (21.0 * C5) - eps) cons = false;

				//Come from application of inequality above
				else if (xl[counter][1] > 0.234 + eps) cons = false;

				else if (xu[counter][5] < 0.166 - eps) cons = false;

				if (cons == true) counter = counter + 1;
			}
		}

		counter = counter - 1;

		if (counter == -1) val = false;

	}

	cout << "Program finished, no negative values found." << endl;

	return 0;
}
