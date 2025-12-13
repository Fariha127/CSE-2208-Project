#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

#define EP 0.01

ofstream fout("output.txt");
int degree;
vector<double> coefficients;

double func(double x) {
    double result = 0.0;
    for (int i = 0; i <= degree; i++) {
        result += coefficients[i] * pow(x, i);
    }
    return result;
}

void bisection(double a, double b) {
    int itr = 0;
    if (func(a) * func(b) >= 0) {
        fout << "No root found in range [" << a << ", " << b << "]" << endl;
        return;
    }

    double c;

    while ((b - a) >= EP) {
        c = (a + b) / 2;
        itr++;

        if (func(c) == 0) {
            break;
        }
        else if (func(a) * func(c) < 0) {
            b = c;
        }
        else {
            a = c;
        }
    }

    fout << "No. of iteration: " << itr << endl;
    fout << "The root is: " << c << endl << endl;
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

    fin >> degree;

    coefficients.resize(degree + 1);

    for (int i = degree; i >=0; i--) {
        fin >> coefficients[i];
    }

    fout << "The polynomial function is: f(x) = ";
    for (int i = degree; i >= 0; i--) {
        if (coefficients[i] != 0) {
            if (i != degree && coefficients[i] > 0) fout << "+";
            fout << coefficients[i];
            if (i > 0) fout << "x^" << i << " ";
        }
    }
    fout << endl << endl;

    double x_max;
    fin >> x_max;

    vector<pair<double, double>> ranges;

    for (double i = -x_max; i < x_max; i += 0.1) {
        double a = func(i);
        double b = func(i + 0.1);
        
        if (a * b < 0) {
            ranges.push_back({i, i + 0.1});
            fout << "Range found: [" << i << ", " << i + 0.1 << "]" << endl;
        }
        
        else if (fabs(a) < EP) {
            fout << "Exact root found at x = " << i << endl;
        }
    }
    
    if (fabs(func(x_max)) < EP) {
        fout << "Exact root found at x = " << x_max << endl;
    }

    fout << endl;

    for (int i = 0; i < ranges.size(); i++) {
        fout << "Root " << i + 1 << ":" << endl;
        bisection(ranges[i].first, ranges[i].second);
    }

    if (ranges.empty()) {
        fout << "No roots found in the given range." << endl;
    }

    return 0;
}