#include<bits/stdc++.h>
#include<fstream>
#include<sstream>
using namespace std;

// reading a CSV file

class csv_reader {
	protected:
	int rows = 0;
	map<int, vector<double>> data;
	public :

	// method to read data from a csv file
	bool readCSV (const string& filename) {
		ifstream file(filename);
		if (!file.is_open()) {
			cerr << "Error opening file! " << endl;
			return false;
		}

		string line;
		while (getline(file, line))
		{
			stringstream ss(line);
			string cell;
			rows++;

			int col = 0;
			while (getline(ss, cell, ',')) {
				double val = stod(cell);
				data[col].push_back(val);
				col++;
			}
		}

		return true;
	}
	
	void printData()
	{
		for (int i = 0; i < rows; i++)
		{
			cout << data[0][i] << " " << data[1][i] << endl;
		}
		
	}

	int getRows () {
		return rows;
	} 

	int getVal (int key, int i) {
		return data[key][i];
	}
};

class preProcessing : protected csv_reader {
	protected :
	vector<vector<bool>> vec;
	public :

	int checkWindow (int s, map<int, vector<double>> map, int day_size, int n, int key) {
		if (s + day_size >= n) return s;
		while (s < s + day_size)
		{
			if (data[key][s] <= 0) return s;
			s++;
		}
		return s;
	}

	void preprocessedWindows (int day_size) {
		int n = getRows();
		
        vec.resize(4, vector<bool>(n - day_size + 1, false));
		for (int key = 1; key < 4; key++)
		{
			int i = 0;
			bool fswitch = true;
			while (i + day_size < n)
			{
				int p = i;
				if (fswitch == true) {
					i = checkWindow(p, data, day_size, n, key);
					fswitch = false;
				}

				if (i == p + day_size)
				{
					while (i < n)
					{
						if (data[1][i] > 0) {
							vec[1][i-day_size] = true; 
							i++;
						}
						else {
							i++;
							fswitch = true;
						}
					}
				}
				else {
					i++;
					fswitch = true;
				}
				
			}
		}
		
		
	}
};

class preChecker : protected preProcessing {
	public :
    
	// Key --> meaning
	// 0 --> time (t)
	// 1 --> displacment (d)
	// 2 --> velocity (v)
	// 3 --> change in velocity (dv)

	bool stageI (int begin,int window) {
		for (int i = begin; i < begin + window; i++)
		{
			if (vec[1][i] <= 0) return false;
		}
		
		return true;
	}

	bool stageII (int begin, int window) {
		for (int i = begin; i < begin + window; i++)
		{
			if (vec[2][i] <= 0) return false;
		}
		
		return true;
	}

	bool stageIII (int begin, int window) {
		double count = 0;
		for (int i = begin; i < begin + window; i++)
		{
			if (vec[3][i] > 0) count = count + 1.00;
		}
		
		if ((count*100/window) >= 75.00) {
			return true;
		}

		return false;
	}
};

class curveFitting : protected preChecker {
	protected :
	// coefficients
    vector<double> prefix_sum_y;
	vector<double> prefix_sum_xy;
	vector<double> prefix_sum_x2y;
	
	public :

	void preComputation(vector<double> velocities)
	{
		int n = velocities.size();
		prefix_sum_y.resize(n, 0.00);
		prefix_sum_xy.resize(n, 0.00);
		prefix_sum_x2y.resize(n, 0.00);

		prefix_sum_y[0] = velocities[0];
		
		for (int i = 1; i < n; i++)
		{
			prefix_sum_y[i] = prefix_sum_y[i-1] + (velocities[i]);
			prefix_sum_xy[i] = prefix_sum_xy[i-1] + (velocities[i] * i);
			prefix_sum_x2y[i] = prefix_sum_x2y[i-1] + (velocities[i] * i * i);
		}
	}

	double gaussElimination (int s, int e, vector<double> velocities) {
		double sum_x = ((e*(e+1)) - (s*(s+1)))/2;
		double sum_x2 = ((e*(e+1)*(2*e + 1)) - (s*(s+1)*(2*s + 1)))/6;
		double sum_x3 = ((e*e * (e+1)*(e+1)) - (s*s * (s+1)*(s+1)))/4;
		double sum_x4 = ((e*(e+1)*(2*e + 1)*(3*e*e + 3*e - 1)) - (s*(s+1)*(2*s + 1)*(3*s*s + 3*s - 1)))/30;

		double sum_y = prefix_sum_y[e] - prefix_sum_y[s] + velocities[s];
		double sum_xy = prefix_sum_xy[e] - prefix_sum_xy[s] + velocities[s];
		double sum_x2y = prefix_sum_x2y[e] - prefix_sum_x2y[s] + velocities[s];
        
		// Finding solution to system of linear equations.
		vector<vector<double>> augmentedMatrix = {
			{(double) e-s+1, sum_x, sum_x2, sum_y},
			{sum_x, sum_x2, sum_x3, sum_xy},
			{sum_x2, sum_x3, sum_x4, sum_x2y}
		};

        for (int i = 0; i < 3; i++)
		{
			// 1. Pivoting
            int max_element = fabs(augmentedMatrix[i][i]);
			int max_row = i;
			for (int j = i+1; j < 3; j++)
			{
				if (max_element < fabs(augmentedMatrix[j][i]))
				{
					max_row = j;
					max_element = fabs(augmentedMatrix[j][i]);
				}
			}
			
			// Swapping current row with max_row

			for (int k = i; k <= 3; k++)
			{
				swap(augmentedMatrix[max_row][k], augmentedMatrix[i][k]);
			}

			// Making all rows below pivoted row 0

			for (int k = i+1; k < 3; k++)
			{
				int mul = -1 * (augmentedMatrix[k][i]/augmentedMatrix[i][i]);
				for (int j = i; j <= 3; j++)
				{
					if (i == j) {
						augmentedMatrix[i][j] = 0;
					}
					else {
						augmentedMatrix[k][j] += mul * augmentedMatrix[i][j];
					}
				}
			}
		}
		
		// Back-substitution
		vector<double> result(3);
		for (int i = 2; i >= 0; i--) {
			result[i] = augmentedMatrix[i][3] / augmentedMatrix[i][i];
			for (int k = i - 1; k >= 0; k--) {
				augmentedMatrix[k][3] -= augmentedMatrix[k][i] * result[i];
			}
		}

		return result[0];
	}
};

int main()
{
	string s;
	cout << "Enter Filename : ";
	cin >> s;

	csv_reader read;
	if (read.readCSV(s)){
		read.printData();

        int n = read.getRows();
		
		// 1 day = 1440 min
		// 1 day = 1440/15 = 96 rows.

		int _size = 96;

		// Set Window
		int window = 5;

        preChecker stages;
		vector<double> velocities;
		for (int i = 0; i + window < n; i++)
		{
			if (stages.stageI (i, window)) {
				if (stages.stageII (i, window-1)) {
					if (stages.stageIII (i, window-1)) {
						for (int k = 0; k < window; i++)
						{
							double avg_velocity = (read.getVal(1, k + _size - 1) - read.getVal(1, k)) / ((double) _size * 15);
							velocities.push_back(avg_velocity);
						}
						
					}
				}
			}
		}

		// curve fitting

		curveFitting parabola;
		int window_2 = 10;
		vector<double> A;
		for (int i = 0; i + window_2 < velocities.size(); i++)
		{
			double a = parabola.gaussElimination(i, i+window_2, velocities);	

			A.push_back(a);
		}
		
		int window_3 = 4;
		for (int i = 0; i + window_3 < A.size(); i++)
		{
			double count = 0.00;
			for (int j = i; j < j + window_3; j++)
			{
				if (A[j] > 0) count++;
			}
			
			if ((count*100)/((double)window_3) >= 75.00) {
				double count2 = 0.00;
				for (int j = i; j < j + window_3 - 1; j++)
				{
					if (A[j] < A[j+1]) count2++;
				}

				if ((count*100)/((double)window_3) >= 75.00) {
					// Employing IVM
				}
			}
		}
		
	}
	else {
		cout << "File already opened !" << endl;
	}


	return 0;
}