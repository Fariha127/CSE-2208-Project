#include <bits/stdc++.h>
using namespace std;

double f(double x) {
    return x*x + 2*x + 1;
}

double simpson38(double a, double b, int n) {
    if (n % 3 != 0) return -1;

    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        double x = a + i*h;
        if (i % 3 == 0)
            sum += 2 * f(x);
        else
            sum += 3 * f(x);
    }

    return (3.0 * h / 8.0) * sum;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin.is_open() || !fout.is_open()) {
        cout << "Error opening file.\n";
        return 0;
    }

    int testCases;
    fin >> testCases;

    for (int t = 1; t <= testCases; t++) {
        double a, b;
        int n;

        fin >> a >> b >> n;

        double result = simpson38(a, b, n);

        cout << fixed << setprecision(3);
        fout << fixed << setprecision(3);

        fout << "Test Case " << t << ":\n";
        cout << "Test Case " << t << ":\n";

        if (result == -1) {
            fout << "Error: Number of intervals (n) must be a multiple of 3.\n\n";
            cout << "Error: Number of intervals (n) must be a multiple of 3.\n\n";
        } else {
            fout << "Approximate integral from " << a << " to " << b << " is = " << result << "\n\n";
            cout << "Approximate integral from " << a << " to " << b << " is = " << result << "\n\n";
        }
    }

    fin.close();
    fout.close();

    return 0;
}
