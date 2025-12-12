#include<bits/stdc++.h>
using namespace std;

double predictExponential(double a, double b, double x) { return a * exp(b * x); }
double predictLogarithmic(double a, double b, double x) { return a + b * log(x); }
double predictPower(double a, double b, double x) { return a * pow(x, b); }

pair<double,double> exponentialRegression(const vector<double>& x, const vector<double>& y) {
    double sumX=0, sumY=0, sumXY=0, sumXX=0;
    int n = x.size();
    for(int i=0;i<n;i++) {
        double Y = log(y[i]);
        sumX += x[i]; sumY += Y; sumXY += x[i]*Y; sumXX += x[i]*x[i];
    }
    double b = (n*sumXY - sumX*sumY) / (n*sumXX - sumX*sumX);
    double a = exp((sumY - b*sumX)/n);
    return {a,b};
}

pair<double,double> logarithmicRegression(const vector<double>& x, const vector<double>& y) {
    double sumLnX=0, sumY=0, sumLnX_Y=0, sumLnX2=0;
    int n = x.size();
    for(int i=0;i<n;i++) {
        double L = log(x[i]);
        sumLnX += L; sumY += y[i]; sumLnX_Y += L*y[i]; sumLnX2 += L*L;
    }
    double b = (n*sumLnX_Y - sumLnX*sumY) / (n*sumLnX2 - sumLnX*sumLnX);
    double a = (sumY - b*sumLnX)/n;
    return {a,b};
}

pair<double,double> powerLawRegression(const vector<double>& x, const vector<double>& y) {
    double sumLnX=0, sumLnY=0, sumLnX_LnY=0, sumLnX2=0;
    int n = x.size();
    for(int i=0;i<n;i++) {
        double LX = log(x[i]); double LY = log(y[i]);
        sumLnX += LX; sumLnY += LY; sumLnX_LnY += LX*LY; sumLnX2 += LX*LX;
    }
    double b = (n*sumLnX_LnY - sumLnX*sumLnY) / (n*sumLnX2 - sumLnX*sumLnX);
    double a = exp((sumLnY - b*sumLnX)/n);
    return {a,b};
}

double computePrediction(int type, const vector<double>& x, const vector<double>& y, double xInput, double &a, double &b) {
    if(type == 1) {
        auto res = exponentialRegression(x,y);
        a = res.first; b = res.second;
        return predictExponential(a,b,xInput);
    } else if(type == 2) {
        auto res = logarithmicRegression(x,y);
        a = res.first; b = res.second;
        return predictLogarithmic(a,b,xInput);
    } else if(type == 3) {
        auto res = powerLawRegression(x,y);
        a = res.first; b = res.second;
        return predictPower(a,b,xInput);
    }
    return 0;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin.is_open() || !fout.is_open()) return 0;

    int testCases;
    fin >> testCases;

    for(int t=0;t<testCases;t++) {
        int type, n;
        fin >> type >> n;

        if(type < 1 || type > 3) {
            for(int i=0;i<2*n+1;i++) { double temp; fin >> temp; }
            continue;
        }

        vector<double> x(n), y(n);
        for(int i=0;i<n;i++) fin >> x[i];
        for(int i=0;i<n;i++) fin >> y[i];

        double xInput;
        fin >> xInput;

        double a=0, b=0;
        double yPred = computePrediction(type, x, y, xInput, a, b);


        string funcStr;
        if(type == 1) funcStr = "y = " + to_string(a) + " * e^(" + to_string(b) + " * x)";
        else if(type == 2) funcStr = "y = " + to_string(a) + " + " + to_string(b) + " * ln(x)";
        else funcStr = "y = " + to_string(a) + " * x^" + to_string(b);

        string typeName = (type==1 ? "Exponential Regression" : type==2 ? "Logarithmic Regression" : "Power-law Regression");

        string output = "Test Case " + to_string(t+1) + ":\n" + typeName + "\nFunction: " + funcStr + "\nPredicted y for x=" + to_string(xInput) + ": " + to_string(yPred) + "\n\n";


        cout << output;
        fout << output;
    }

    fin.close();
    fout.close();
    return 0;
}
