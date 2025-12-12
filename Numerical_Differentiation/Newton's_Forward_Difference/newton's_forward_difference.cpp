#include <bits/stdc++.h>
using namespace std;

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
        for (int j = 0; j < n - i; j++)
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
    }

    int i = 0;
    float u = (value - x[i]) / h;
    
    float derivative = y[i][1] / h;

    if (n >= 3) {
        derivative += ((2 * u - 1) / 2.0) * (y[i][2] / h);
    }

    if (n >= 4) {
        derivative += ((3 * u * u - 6 * u + 2) / 6.0) * (y[i][3] / h);
    }

    if (n >= 5) {
        derivative += ((4 * u * u * u - 18 * u * u + 22 * u - 6) / 24.0) * (y[i][4] / h);
    }

    fout << "\nAt x = " << value << ", Derivative = " << derivative << endl;

    fin.close();
    fout.close();

    return 0;
}
