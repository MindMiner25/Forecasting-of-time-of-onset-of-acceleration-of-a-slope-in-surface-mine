#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void gaussianElimination(vector<vector<double>> &augmentedMatrix, int n) {
    // Forward elimination process
    for (int i = 0; i < n; ++i) {
        // Search for maximum in this column
        double maxElement = fabs(augmentedMatrix[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (fabs(augmentedMatrix[k][i]) > maxElement) {
                maxElement = fabs(augmentedMatrix[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (pivoting)
        for (int k = i; k <= n; ++k) {
            swap(augmentedMatrix[maxRow][k], augmentedMatrix[i][k]);
        }

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; ++k) {
            double c = -augmentedMatrix[k][i] / augmentedMatrix[i][i];
            for (int j = i; j <= n; ++j) {
                if (i == j) {
                    augmentedMatrix[k][j] = 0;
                } else {
                    augmentedMatrix[k][j] += c * augmentedMatrix[i][j];
                }
            }
        }
    }

    // Back-substitution
    vector<double> result(n);
    for (int i = n - 1; i >= 0; i--) {
        result[i] = augmentedMatrix[i][n] / augmentedMatrix[i][i];
        for (int k = i - 1; k >= 0; k--) {
            augmentedMatrix[k][n] -= augmentedMatrix[k][i] * result[i];
        }
    }

    // Output the result
    cout << "Solution: " << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << result[i] << endl;
    }
}

int main() {
    int n = 3;  // Number of variables (or equations)
    vector<vector<double>> augmentedMatrix = {
        {5, 10, 30, 17.1},
        {10, 30, 100, 47.5},
        {30, 100, 354, 157.5}
    };

    gaussianElimination(augmentedMatrix, n);
    return 0;
}
