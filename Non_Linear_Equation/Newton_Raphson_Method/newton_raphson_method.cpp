#include <bits/stdc++.h>

#define E 0.001

using namespace std;


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

double d_func(double x) {
    double result = 0.0;
    for (int i = 1; i <= degree; i++) {
        result += i * coefficients[i] * pow(x, i - 1);
    }
    return result;
}

void newton_raphson(double x, int root_num) {

    int itr = 0;
    double h;

    do {
        h = func(x) / d_func(x);
        x = x - h;
        itr++;
    } while (fabs(h) >= E);

    fout << "Root " << root_num << ":" << endl;
    fout << "The root is: " << x << endl;
    fout << "Iteration needed: " << itr << endl << endl;
}

int main()
{
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

    vector<double> initial_guesses(degree);
    
    for (int i = 0; i < degree; i++) {
        fin >> initial_guesses[i];
    }

    for (int i = 0; i < degree; i++) {
        newton_raphson(initial_guesses[i], i + 1);
    }

    return 0;
}
 