
#include <bits/stdc++.h>
using namespace std;

vector<double> gaussianElimination(vector<vector<double>>& mat, vector<double>& y) {
    int n = mat.size();
    for(int i = 0; i < n; i++) {
        int maxRow = i;
        for(int k = i+1; k < n; k++)
            if(fabs(mat[k][i]) > fabs(mat[maxRow][i])) maxRow = k;
        swap(mat[i], mat[maxRow]);
        swap(y[i], y[maxRow]);

        for(int k = i+1; k < n; k++) {
            double factor = mat[k][i]/mat[i][i];
            for(int j = i; j < n; j++)
                mat[k][j] -= factor*mat[i][j];
            y[k] -= factor*y[i];
        }
    }

    vector<double> x(n);
    for(int i = n-1; i >= 0; i--) {
        x[i] = y[i];
        for(int j = i+1; j < n; j++)
            x[i] -= mat[i][j]*x[j];
        x[i] /= mat[i][i];
    }
    return x;
}

vector<double> polynomialRegression(const vector<double>& x,
                                    const vector<double>& y,
                                    int degree) {
    int n = x.size();
    int d = degree;

    vector<vector<double>> X(n, vector<double>(d+1));
    for(int i = 0; i < n; i++) {
        X[i][0] = 1.0;
        for(int k = 1; k <= d; k++)
            X[i][k] = pow(x[i], k);
    }

    vector<vector<double>> XtX(d+1, vector<double>(d+1,0.0));
    vector<double> Xty(d+1,0.0);

    for(int i = 0; i <= d; i++) {
        for(int j = 0; j <= d; j++)
            for(int k = 0; k < n; k++)
                XtX[i][j] += X[k][i]*X[k][j];

        for(int k = 0; k < n; k++)
            Xty[i] += X[k][i]*y[k];
    }

    return gaussianElimination(XtX, Xty);
}

double predictPolynomial(const vector<double>& coeff, double x) {
    double y = 0, term = 1;
    for(double a : coeff) {
        y += a*term;
        term *= x;
    }
    return y;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin.is_open() || !fout.is_open()) {
        cout << "Error opening file\n";
        return 0;
    }

    int testCases;
    fin >> testCases;

    for(int t = 1; t <= testCases; t++) {
        int n, degree;
        fin >> n >> degree;

        vector<double> x(n), y(n);
        for(int i = 0; i < n; i++) fin >> x[i] >> y[i];

        double xInput;
        fin >> xInput;

        vector<double> coeff = polynomialRegression(x, y, degree);
        double yPred = predictPolynomial(coeff, xInput);

        stringstream ss;
        ss << fixed << setprecision(4);
        ss << "Test Case " << t << ":\n";
        ss << "Polynomial Degree: " << degree << "\n";
        ss << "Function: y = ";
        for(int i = 0; i <= degree; i++) {
            if(i > 0) ss << " + ";
            ss << coeff[i];
            if(i > 0) ss << "*x^" << i;
        }
        ss << "\nPredicted y for x = " << xInput << ": " << yPred << "\n\n";

        string output = ss.str();
        cout << output;
        fout << output;
    }

    fin.close();
    fout.close();
    return 0;
}
