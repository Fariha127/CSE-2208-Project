#include<bits/stdc++.h>
#include<cmath>
using namespace std;
#define E 0.001

ofstream fout("output.txt");
int degree;
vector<double> coefficients;

double f(double x) {
    double result = 0.0;
    for (int i = 0; i <= degree; i++) {
        result += coefficients[i] * pow(x, i);
    }
    return result;
}

void secant (double x1, double x2, int root_num) {

    int itr = 0;
    double x0;

    do {
        x0 = ((x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1)));

        x1 = x2;
        x2 = x0;

        itr++;

        if(f(x0) == 0) break;

    } while (fabs(f(x0)) >= E);

    fout << "Root " << root_num << ":" << endl;
    fout << "The root is: " << x0 << endl;
    fout << "Iteration needed: " << itr << endl << endl;

}


int main () {

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

    vector<pair<double, double>> initial_guesses(degree);
    
    fout << "Initial guess pairs: ";
    for (int i = 0; i < degree; i++) {
        fin >> initial_guesses[i].first >> initial_guesses[i].second;
        fout << "(" << initial_guesses[i].first << ", " << initial_guesses[i].second << ")";
        if (i < degree - 1) fout << ", ";
    }
    fout << endl << endl;

    for (int i = 0; i < degree; i++) {
        secant(initial_guesses[i].first, initial_guesses[i].second, i + 1);
    }

    return 0;
}
