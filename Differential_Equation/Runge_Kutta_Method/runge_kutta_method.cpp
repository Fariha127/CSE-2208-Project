#include <bits/stdc++.h>
using namespace std;

// Differential equation: dy/dx = x^2 - y
double f(double x, double y) {
    return x*x - y;
}

double rungeKutta(double x0, double y0, double h, double x) {
    int n = (x - x0) / h;
    double k1, k2, k3, k4;

    for (int i = 0; i < n; i++) {
        k1 = h * f(x0, y0);
        k2 = h * f(x0 + h/2, y0 + k1/2);
        k3 = h * f(x0 + h/2, y0 + k2/2);
        k4 = h * f(x0 + h, y0 + k3);

        y0 = y0 + (k1 + 2*k2 + 2*k3 + k4) / 6;
        x0 = x0 + h;
    }

    return y0;
}

int main() {
    ifstream inFile("input.txt");
    ofstream outFile("output.txt");

    if (!inFile) {
        cout << "Error: Could not open Input.txt" << endl;
        return 1;
    }
    int T;
    inFile >> T;
     for (int tc = 1; tc <= T; tc++) {

    double x0, y0, h, x;

    inFile >> x0 >> y0 >> h >> x;

    double result = rungeKutta(x0, y0, h, x);

    cout << fixed << setprecision(4);
    cout << "Result of y at x = " << x << "is : " << result << endl;

    outFile << fixed << setprecision(6);
    outFile << "Result of y at x = " << x << "is : " << result << endl;
     }

    inFile.close();
    outFile.close();

    return 0;
}

