#include <bits/stdc++.h>
using namespace std;


float v_term(float u, int n)
{
    float temp = u;
    for (int i = 1; i < n; i++)
        temp = temp * (u - i);
    return temp;
}


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

   
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n - i; j++)
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
    }

    
    fout << "Forward Difference Table:\n\n";
    for (int i = 0; i < n; i++) {
        fout << setw(8) << x[i] << "\t";
        for (int j = 0; j < n - i; j++)
            fout << setw(10) << y[i][j] << "\t";
        fout << endl;
    }

    
    float sum = y[0][0];
    float u = (value - x[0]) / (x[1] - x[0]);
    
    for (int i = 1; i < n; i++) {
        sum = sum + (v_term(u, i) * y[0][i]) / factorial(i);
    }

    fout << "\nInterpolated value at X = " << value << " is Y = " << sum << endl;

    fin.close();
    fout.close();

    return 0;
}
