#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

ofstream fout("output.txt");
void printMatrix(double M[][20], int n, const string &name)
{
  fout << "\nMatrix " << name << ":"<<endl;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      fout << setw(10) << fixed << setprecision(4) << M[i][j] << " ";
    fout << endl;
  }
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

    int n;
    fin >> n;
  double A[20][20], L[20][20] = {0}, U[20][20] = {0}, b[20], y[20] = {0}, x[20] = {0};

  for (int i = 0; i <= n; i++)
    for (int j = 0; j <= n; j++)
      fin >> A[i][j];

  for (int i = 0; i < n; i++)
     b[i] = A[i][n];

  for (int i = 0; i < n; i++)
    L[i][i] = 1;

  for (int i = 0; i < n; i++)
  {
    for (int k = i; k < n; k++)
    {
      double sum = 0;
      for (int j = 0; j < i; j++)
        sum += L[i][j] * U[j][k];
      U[i][k] = A[i][k] - sum;
      fout << "\nAfter Computing U[" << i + 1 << "][" << k + 1 << "]" << ":";
      printMatrix(L, n, "L");
      printMatrix(U, n, "U");
    }

    for (int k = i + 1; k < n; k++)
    {
      double sum = 0;
      for (int j = 0; j < i; j++)
        sum += L[k][j] * U[j][i];
      L[k][i] = (A[k][i] - sum) / U[i][i];
      fout << "\nAfter Computing L[" << k + 1 << "][" << i + 1 << "]" << ":";
      printMatrix(L, n, "L");
      printMatrix(U, n, "U");
    }
  }

  fout << "\nFinal L and U matrices:";
  printMatrix(L, n, "L");
  printMatrix(U, n, "U");

  bool singular = false;
  for (int i = 0; i < n; i++)
  {
    if (fabs(U[i][i]) < 1e-9)
      singular = true;
  }

  fout << "\nForward Substitution (Ly = b)"<<endl;
  for (int i = 0; i < n; i++)
  {
    double sum = 0;
    for (int j = 0; j < i; j++)
      sum += L[i][j] * y[j];
    y[i] = b[i] - sum;
    fout << "y" << i + 1 << " = " << fixed << setprecision(6) << y[i] << endl;
  }

  bool noSolution = false, infiniteSolution = false;

  for (int i = 0; i < n; i++)
  {
    bool rowZero = true;
    for (int j = 0; j < n; j++)
    {
      if (fabs(U[i][j]) > 1e-9)
      {
        rowZero = false;
        break;
      }
    }

    if (rowZero && fabs(y[i]) > 1e-9)
      noSolution = true;
    else if (rowZero && fabs(y[i]) < 1e-9)
      infiniteSolution = true;
  }

  if (noSolution)
  {
    fout << "\nThe system has No Solution."<<endl;
    return 0;
  }
  if (infiniteSolution)
  {
    fout << "\nThe system has Infinite Solutions."<<endl;
    return 0;
  }

  fout << "\nBackward Substitution (Ux = y)"<<endl;
  for (int i = n - 1; i >= 0; i--)
  {
    double sum = 0;
    for (int j = i + 1; j < n; j++)
      sum += U[i][j] * x[j];
    x[i] = (y[i] - sum) / U[i][i];
    fout << "x" << i + 1 << " = " << fixed << setprecision(6) << x[i] << endl;
  }

  fout << "\nThe system has a Unique Solution."<<endl;
  return 0;
}