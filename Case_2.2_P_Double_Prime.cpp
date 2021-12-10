#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

//Shows that the program P'' doesn't have any negative values
//Does an iterated grid search - take some portion of the search space and bound the value the objective function can have if it's positive move on to the next portion otherwise split into smaller pieces and repeat
//The program won't terminate if the solution of the program could be negative, otherwise it will termninate 

int main() {

	int i, j, k, n, counter, numIt, temp, bin;
	double C5, eps, f, M1, M2, M3, temp1, temp2, minVal, d;
	bool val, cons;

	double xlTemp[6], xuTemp[6], blTemp[6], buTemp[6];

	//ruT is between [0,1] but ru is between [0,x] 
	//ruT_i x_i = ru
	//The T means total
	//End up need needing more than 5000 vectors
	std::vector < vector<double> > xu(5000, vector<double>(6)); //xu[.,0] = 0, xu[.,i] = x_i + epsilon
	std::vector < vector<double> > xl(5000, vector<double>(6));
	std::vector < vector<double> > ru(5000, vector<double>(6));
	std::vector < vector<double> > rl(5000, vector<double>(6));
	std::vector < vector<double> > bu(5000, vector<double>(6));
	std::vector < vector<double> > bl(5000, vector<double>(6));
	std::vector < vector<double> > buT(5000, vector<double>(6));
	std::vector < vector<double> > blT(5000, vector<double>(6));

	n = 1000;
	C5 = 0.0384609;
	d = 0.2;
	val = true;
	counter = 0; //Determines where in stack (list of vectors above) we are
	numIt = 0;
	eps = 0.00001; //Tolerance
	f = (2 * n) / (double)(n - 1) * (10 / ((double)25) - 2.0 * (-0.175431374077117 + 8.75407592662244 * C5) / (21 * C5));

	for (i = 1; i < 6; ++i) {
		xl[counter][i] = 0.166 - eps;
		rl[counter][i] = 0.0;
		bl[counter][i] = 0.0;
		blT[counter][i] = 0.0;

		xu[counter][i] = 0.234 + eps;
		ru[counter][i] = 0.234 + eps;
		bu[counter][i] = 0.234 + eps;
		buT[counter][i] = 1.0;
	}

	minVal = 0.0;

	while (val == true) {
		//Gives updates on running status
		//Comment out if you don't want them
		/*if (numIt % 10000 == 0) {
			cout << "The iteration number is: " << numIt << endl;
			cout << "Counter = " << counter << endl;
			cout << "The vector xl at the current step is: " << xl[counter][1] << " " << xl[counter][2] << " " << xl[counter][3] << " " << xl[counter][4] << " " << xl[counter][5] << endl;
			cout << "The vector xu at the current step is: " << xu[counter][1] << " " << xu[counter][2] << " " << xu[counter][3] << " " << xu[counter][4] << " " << xu[counter][5] << endl;
			cout << "The vector bl at the current step is: " << bl[counter][1] << " " << bl[counter][2] << " " << bl[counter][3] << " " << bl[counter][4] << " " << bl[counter][5] << endl;
			cout << "The vector bu at the current step is: " << bu[counter][1] << " " << bu[counter][2] << " " << bu[counter][3] << " " << bu[counter][4] << " " << bu[counter][5] << endl;
			cout << "The vector rl at the current step is: " << rl[counter][1] << " " << rl[counter][2] << " " << rl[counter][3] << " " << rl[counter][4] << " " << rl[counter][5] << endl;
			cout << "The vector ru at the current step is: " << ru[counter][1] << " " << ru[counter][2] << " " << ru[counter][3] << " " << ru[counter][4] << " " << ru[counter][5] << endl << endl;
		}*/

		numIt = numIt + 1;


		//Calculates the objective function
		M1 = 0.0;
		for (i = 1; i < 6; ++i) {
			temp1 = 1;
			for (j = 1; j < 6; ++j) {
				if (i != j) temp1 = temp1 * xl[counter][j];
			}

			temp2 = 0;
			for (j = 1; j < 6; ++j) {
				for (k = 1; k < 6; ++k) {
					if ((i != j) && (i != k) && (j != k)) temp2 = std::max(temp2, xu[counter][j] * xu[counter][k]);
				}
			}

			M1 = std::max(M1, temp1 - f * temp2 / ((double)2));
		}

		M2 = 0.0;
		for (i = 1; i < 6; ++i) {
			M2 = M2 + ru[counter][(1 + i) % 5] * bu[counter][(2 + i) % 5] * bu[counter][(3 + i) % 5] * ru[counter][(4 + i) % 5];
		}

		for (i = 2; i < 6; ++i) {
			M2 = M2 + ru[counter][i] * ru[counter][i] * bu[counter][i] * bu[counter][i] / ((double)16);
		}

		M3 = 0;
		for (i = 0; i < 5; ++i) {
			if (bu[counter][(1 + i) % 5] * bu[counter][(2 + i) % 5] + bu[counter][(1 + i) % 5] * bu[counter][(5 + i) % 5] + bu[counter][(3 + i) % 5] * bu[counter][(2 + i) % 5] > M3) {
				M3 = bu[counter][(1 + i) % 5] * bu[counter][(2 + i) % 5] + bu[counter][(1 + i) % 5] * bu[counter][(5 + i) % 5] + bu[counter][(3 + i) % 5] * bu[counter][(2 + i) % 5];
			}

			if (ru[counter][(5 + i) % 5] * bu[counter][(2 + i) % 5] + ru[counter][(2 + i) % 5] * bu[counter][(2 + i) % 5] > M3) {
				M3 = ru[counter][(5 + i) % 5] * bu[counter][(2 + i) % 5] + ru[counter][(2 + i) % 5] * bu[counter][(2 + i) % 5];
			}

			if (bu[counter][(1 + i) % 5] * ru[counter][(3 + i) % 5] + bu[counter][(1 + i) % 5] * ru[counter][(1 + i) % 5] > M3) {
				M3 = bu[counter][(1 + i) % 5] * ru[counter][(3 + i) % 5] + bu[counter][(1 + i) % 5] * ru[counter][(1 + i) % 5];
			}

			if (ru[counter][(1 + i) % 5] * ru[counter][(3 + i) % 5] + ru[counter][(5 + i) % 5] * ru[counter][(3 + i) % 5] + ru[counter][(1 + i) % 5] * ru[counter][(4 + i) % 5] > M3) {
				M3 = ru[counter][(1 + i) % 5] * ru[counter][(3 + i) % 5] + ru[counter][(5 + i) % 5] * ru[counter][(3 + i) % 5] + ru[counter][(1 + i) % 5] * ru[counter][(4 + i) % 5];
			}

			if (bu[counter][(4 + i) % 5] * ru[counter][(3 + i) % 5] + bu[counter][(3 + i) % 5] * ru[counter][(3 + i) % 5] > M3) {
				M3 = bu[counter][(4 + i) % 5] * ru[counter][(3 + i) % 5] + bu[counter][(3 + i) % 5] * ru[counter][(3 + i) % 5];
			}

			if (bu[counter][(5 + i) % 5] * ru[counter][(1 + i) % 5] + bu[counter][(1 + i) % 5] * ru[counter][(1 + i) % 5] > M3) {
				M3 = bu[counter][(5 + i) % 5] * ru[counter][(1 + i) % 5] + bu[counter][(1 + i) % 5] * ru[counter][(1 + i) % 5];
			}
		}

		minVal = M1 - M2 - (1 + 4 * M3) * f / 16.0;

		if (minVal < -10000) {
			cout << "Found a negative value." << endl;
			cout << "The minimum value is: " << minVal << endl;
			cout << "The iteration number is: " << numIt << endl;
			cout << "Counter = " << counter << endl;
			cout << "The vector xl at the current step is: " << xl[counter][1] << " " << xl[counter][2] << " " << xl[counter][3] << " " << xl[counter][4] << " " << xl[counter][5] << endl;
			cout << "The vector xu at the current step is: " << xu[counter][1] << " " << xu[counter][2] << " " << xu[counter][3] << " " << xu[counter][4] << " " << xu[counter][5] << endl;
			cout << "The vector bl at the current step is: " << bl[counter][1] << " " << bl[counter][2] << " " << bl[counter][3] << " " << bl[counter][4] << " " << bl[counter][5] << endl;
			cout << "The vector bu at the current step is: " << bu[counter][1] << " " << bu[counter][2] << " " << bu[counter][3] << " " << bu[counter][4] << " " << bu[counter][5] << endl;
			cout << "The vector rl at the current step is: " << rl[counter][1] << " " << rl[counter][2] << " " << rl[counter][3] << " " << rl[counter][4] << " " << rl[counter][5] << endl;
			cout << "The vector ru at the current step is: " << ru[counter][1] << " " << ru[counter][2] << " " << ru[counter][3] << " " << ru[counter][4] << " " << ru[counter][5] << endl << endl;
		}

		//Grid Stuff 
		if (minVal < eps) {
			for (i = 1; i < 6; ++i) {
				xlTemp[i] = xl[counter][i];
				xuTemp[i] = xu[counter][i];
				blTemp[i] = blT[counter][i];
				buTemp[i] = buT[counter][i];
			}

			/*
			Bin is the binary representation of i
			Depending whether the ith digit is 0 or 1 we will either
			keep xl[i] the same and make xu[i] = (xu[i]+xl[i})/2 or vice versa, then base b_i off that
			*/
			for (i = 0; i < pow(2, 9); ++i) {
				temp = i;
				cons = true;

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

				for (j = 4; j < 9; ++j) {
					bin = temp % 2;
					temp = temp / 2;

					if (bin == 0) {
						blT[counter][j - 3] = blTemp[j - 3];
						buT[counter][j - 3] = (buTemp[j - 3] + blTemp[j - 3]) / ((double)2);
					}

					else {
						blT[counter][j - 3] = (buTemp[j - 3] + blTemp[j - 3]) / ((double)2);
						buT[counter][j - 3] = buTemp[j - 3];
					}
				}

				//Can define x_5 based on other x_i's
				xl[counter][5] = std::max(0.0, 1 - xu[counter][1] - xu[counter][2] - xu[counter][3] - xu[counter][4]) -eps; //
				xu[counter][5] = 1 - xl[counter][1] - xl[counter][2] - xl[counter][3] - xl[counter][4] + eps; //

				//Can define r_i's based on b_i's and do without T at same time
				for (j = 1; j < 6; ++j) {
					bl[counter][j] = blT[counter][j] * xl[counter][j] - eps; //
					bu[counter][j] = buT[counter][j] * xu[counter][j] + eps; //
					rl[counter][j] = std::max(0.0, xl[counter][j] - bu[counter][j]) - eps; //
					ru[counter][j] = std::min(1.0, xu[counter][j] - bl[counter][j]) + eps; //
				}

				//Check all constraints are satisfied
				//Check whether x5 has any feasible points
				if (xu[counter][5] < eps) cons = false;

				else if (xu[counter][1] * xu[counter][2] + xu[counter][1] * xu[counter][3] + xu[counter][1] * xu[counter][4] + xu[counter][1] * xu[counter][5]
					+ xu[counter][2] * xu[counter][3] + xu[counter][2] * xu[counter][4] + xu[counter][2] * xu[counter][5] + xu[counter][3] * xu[counter][4] + xu[counter][3] * xu[counter][5] + xu[counter][4] * xu[counter][5]
					< 2.0 * (-0.175431374077117 + 8.75407592662244 * C5) / (21.0 * C5) - eps) cons = false;

				//Come from application of inequality above
				else if (xl[counter][1] > 0.234 + eps) cons = false;

				else if (xu[counter][5] < 0.166 - eps) cons = false;

				//Comes from funky degree bounds
				else if (bu[counter][2] + ru[counter][3] + ru[counter][4] + bu[counter][5] < d / 2.0 + eps) cons = false;

				//Proper ordering
				if (cons == true) {
					for (j = 1; j < 5; ++j) {
						for (k = j + 1; k < 6; ++k) {
							if (xu[counter][j] < xl[counter][k] - eps)  cons = false; 
						}
					}
				}

				//Comes from move to x_1 being best possible
				if (cons == true) {
					temp1 = 0.0;
					for (j = 2; j < 6; ++j) {
						temp1 = std::max(temp1, rl[counter][(1 + j) % 5] + bl[counter][(j + 2) % 5] + bl[counter][(j + 3) % 5] + rl[counter][(j + 4) % 5]);
					}

					if (ru[counter][2] + bu[counter][3] + bu[counter][4] + ru[counter][5] < temp1 + eps) cons = false;
				}

				if (cons == true) counter = counter + 1;
			}
		}

		counter = counter - 1;

		if (counter == -1) val = false;
	}

	cout << "Program finished, no negative values found." << endl;
	return 0;
}
