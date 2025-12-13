#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

ofstream fout("output.txt");

void printMatrix(double a[][20], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            fout << setw(10) << fixed << setprecision(4) << a[i][j] << " ";
        }
        fout << endl;
    }
}

int main() {
  ifstream fin("input.txt");

    if (!fin) {
        fout << "Error: Cannot open input.txt"<<endl;
        return 1;
    }
    if (!fout) {
        fout << "Error: Cannot open output.txt"<<endl;
        return 1;
    }

    int n;
    fin >> n;

    double a[20][20], x[20];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            fin >> a[i][j];
        }
    }

    fout << "\nInitial Matrix:"<<endl;
    printMatrix(a, n);

    for (int i = 0; i < n - 1; i++) {

        for (int k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[i][i])) {
                for (int j = 0; j <= n; j++) {
                    swap(a[i][j], a[k][j]);
                }
                fout << "\nAfter swapping row " << i + 1 << " and row " << k + 1 << ":"<<endl;
                printMatrix(a, n);
            }
        }

        if (fabs(a[i][i]) < 1e-9) continue;

        for (int k = i + 1; k < n; k++) {
            double factor = a[k][i] / a[i][i];
            fout << "\nR" << k + 1 << " = R" << k + 1
                 << " - (" << factor << ") * R" << i + 1 << endl;

            for (int j = i; j <= n; j++) {
                a[k][j] -= factor * a[i][j];
            }
            printMatrix(a, n);
        }
    }

    fout << "\nMatrix after Forward Elimination (Upper Triangular):"<< endl;
    printMatrix(a, n);

    bool noSolution = false, infinite = false;

    for (int i = 0; i < n; i++) {
        bool allZero = true;
        for (int j = 0; j < n; j++) {
            if (fabs(a[i][j]) > 1e-9) {
                allZero = false;
                break;
            }
        }

        if (allZero && fabs(a[i][n]) > 1e-9) {
            noSolution = true;
            break;
        }
        else if (allZero && fabs(a[i][n]) < 1e-9) {
            infinite = true;
        }
    }

    if (noSolution) {
        fout << "\nThe system has No Solution (Inconsistent system)." << endl;
        return 0;
    }
    else if (infinite) {
        fout << "\nThe system has Infinite Solutions (Dependent equations)." << endl;
        return 0;
    }

    fout << "\nBack Substitution Steps:"<< endl;

    for (int i = n - 1; i >= 0; i--) {
        double sum = a[i][n];
        fout << "x" << i + 1 << " = ( " << a[i][n];

        for (int j = i + 1; j < n; j++) {
            sum -= a[i][j] * x[j];
            fout << " - (" << a[i][j] << " * x" << j + 1 << ")";
        }

        x[i] = sum / a[i][i];
        fout << " ) / " << a[i][i] << endl;
        fout << "x" << i + 1 << " = " << fixed << setprecision(6) << x[i] << endl << endl;
    }

    fout << "\nUnique Solution:"<<endl;
    for (int i = 0; i < n; i++) {
        fout << "x" << i + 1 << " = " << fixed << setprecision(6) << x[i] << endl;
    }

    return 0;
}
