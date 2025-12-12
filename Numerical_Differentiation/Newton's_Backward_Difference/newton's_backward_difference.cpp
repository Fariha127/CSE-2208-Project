#include <bits/stdc++.h>
using namespace std;

int factorial(int n)
{
    int f = 1;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

int main()
{
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

    int n;
    fin >> n;

    float x[n], y[n][n];

    for (int i = 0; i < n; i++) {
        fin >> x[i] >> y[i][0];
    }

    float value;
    fin >> value;

    float h = x[1] - x[0];

    for (int i = 1; i < n; i++) {
        for (int j = n - 1; j >= i; j--)
            y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
    }

    int i = n - 1;
    float v = (value - x[i]) / h;

    float derivative = y[i][1] / h;
    
    if (n >= 3) {
        derivative += ((2 * v + 1) / 2.0) * (y[i][2] / h);
    }
    
    if (n >= 4) {
        derivative += ((3 * v * v + 6 * v + 2) / 6.0) * (y[i][3] / h);
    }
    
    if (n >= 5) {
        derivative += ((4 * v * v * v + 18 * v * v + 22 * v + 6) / 24.0) * (y[i][4] / h);
    }

    fout << "\nAt x = " << value << ", Derivative = " << derivative << endl;

    fin.close();
    fout.close();

    return 0;
}
