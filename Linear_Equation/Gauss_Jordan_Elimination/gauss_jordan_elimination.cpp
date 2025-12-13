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

    double a[20][20];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            fin >> a[i][j];
        }
    }

    fout << "\nInitial Matrix:"<<endl;
    printMatrix(a, n);

    for (int i = 0; i < n; i++) {
        for (int k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[i][i])) {
                for (int j = 0; j <= n; j++) {
                    swap(a[i][j], a[k][j]);
                }
                fout << "\nAfter swapping row " << i + 1 << " and row " << k + 1 << ":"<<endl;
                printMatrix(a, n);
            }
        }

        double pivot = a[i][i];
        if (fabs(pivot) < 1e-9) continue;

        for (int j = 0; j <= n; j++) {
            a[i][j] /= pivot;
        }
        fout << "\nAfter dividing row " << i + 1 << " by pivot (" << pivot << "):"<<endl;
        printMatrix(a, n);

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = a[k][i];
                for (int j = 0; j <= n; j++) {
                    a[k][j] -= factor * a[i][j];
                }
            }
        }
        fout << "\nAfter eliminating column " << i + 1 << ":"<<endl;
        printMatrix(a, n);
    }

    fout << "\nFinal Matrix (Reduced Row Echelon Form):"<<endl;
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
    }
    else if (infinite) {
        fout << "\nThe system has Infinite Solutions (Dependent equations)." << endl;
    }
    else {
        fout << "\nThe system has a Unique Solution:"<<endl;
        for (int i = 0; i < n; i++) {
            fout << "x" << i + 1 << " = " << fixed << setprecision(6) << a[i][n] << endl;
        }
    }

    return 0;
}