#include <bits/stdc++.h>
using namespace std;

double predict(double a, double b, double x) {
    return a + b * x;
}

bool computeRegression(
    const vector<double>& x,
    const vector<double>& y,
    double& a,
    double& b
) {
    int n = x.size();
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    double denominator = (n * sumX2) - (sumX * sumX);
    if (denominator == 0) return false;

    b = (n * sumXY - sumX * sumY) / denominator;
    a = (sumY - b * sumX) / n;

    return true;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin.is_open()) {
        cout << "Error: Could not open input file.\n";
        return 0;
    }

    if (!fout.is_open()) {
        cout << "Error: Could not open output file.\n";
        return 0;
    }
int n, testCase = 1;

while (true) {
    if (!(fin >> n)) break;

    vector<double> x(n), y(n);

    for (int i = 0; i < n; i++) {
        fin >> x[i] >> y[i];
    }

    double X;
    fin >> X;

    cout << "Test Case " << testCase << "\n";
    fout << "Test Case " << testCase << "\n";

    double a, b;
    if (!computeRegression(x, y, a, b)) {
        cout << "Error: Regression cannot be computed.\n\n";
        fout << "Error: Regression cannot be computed.\n\n";
        testCase++;
        continue;
    }

    cout << "y = " << a << " + " << b << "x\n";
    fout << "y = " << a << " + " << b << "x\n";

    double Y = predict(a, b, X);

    cout << "Predicted Y = " << Y << " for X = " << X << "\n\n";
    fout << "Predicted Y = " << Y << " for X = " << X << "\n\n";

    testCase++;
}


    fin.close();
    fout.close();

    return 0;
}
