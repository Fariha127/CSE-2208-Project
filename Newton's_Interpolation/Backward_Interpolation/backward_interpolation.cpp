#include <bits/stdc++.h>
using namespace std;


double u_term(double u, int n)
{
    double temp = u;
    for (int i = 1; i < n; i++)
        temp = temp * (u + i);
    return temp;
}


int factorial(int n)
{
    int f = 1;
    for (int i = 2; i <= n; i++)
        f = f * i;
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

    double value;
    fin >> value;
    
    for (int i = 1; i < n; i++) {
        for (int j = n - 1; j >= i; j--)
            y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
    }

    fout << "Backward Difference Table:\n\n";
    for (int i = 0; i < n; i++) {
        fout << setw(8) << x[i] << "\t";
        for (int j = 0; j <= i; j++)
            fout << setw(10) << y[i][j] << "\t";
        fout << endl;
    }

    double sum = y[n - 1][0];
    double u = (value - x[n - 1]) / (x[1] - x[0]);
    
    for (int i = 1; i < n; i++) {
        sum = sum + (u_term(u, i) * y[n - 1][i]) / factorial(i);
    }

    fout << "\nInterpolated value at X = " << value << " is Y = " << sum << endl;

    fin.close();
    fout.close();

    return 0;
}
