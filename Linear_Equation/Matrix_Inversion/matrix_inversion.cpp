#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;

const double EPS = 1e-9;

double determinant(vector<vector<double>> A) {
    int n = A.size();
    double det = 1;
    for (int i = 0; i < n; i++) {
        int pivot = i;
        for (int r = i; r < n; r++)
            if (fabs(A[r][i]) > fabs(A[pivot][i])) pivot = r;

        if (fabs(A[pivot][i]) < EPS) return 0;

        if (i != pivot) swap(A[i], A[pivot]), det *= -1;

        det *= A[i][i];

        double div = A[i][i];
        for (int j = i; j < n; j++) A[i][j] /= div;

        for (int r = i + 1; r < n; r++) {
            double factor = A[r][i];
            for (int c = i; c < n; c++)
                A[r][c] -= factor * A[i][c];
        }
    }
    return det;
}

void getCofactor(const vector<vector<double>>& A, vector<vector<double>>& temp, int p, int q) {
    int n = A.size();
    int row = 0, col = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != p && j != q) {
                temp[row][col++] = A[i][j];
                if (col == n - 1) col = 0, row++;
            }
}

vector<vector<double>> adjoint(const vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> adj(n, vector<double>(n));

    if (n == 1) {
        adj[0][0] = 1;
        return adj;
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            vector<vector<double>> temp(n - 1, vector<double>(n - 1));
            getCofactor(A, temp, i, j);
            double sign = ((i + j) % 2 == 0) ? 1 : -1;
            adj[j][i] = sign * determinant(temp);
        }

    return adj;
}

int rankMatrix(vector<vector<double>> A) {
    int n = A.size(), m = A[0].size();
    int rank = 0;

    for (int col = 0, row = 0; col < m && row < n; col++) {
        int sel = row;
        for (int i = row; i < n; i++)
            if (fabs(A[i][col]) > fabs(A[sel][col])) sel = i;

        if (fabs(A[sel][col]) < EPS) continue;

        swap(A[sel], A[row]);

        for (int j = col + 1; j < m; j++)
            A[row][j] /= A[row][col];
        A[row][col] = 1;

        for (int i = 0; i < n; i++)
            if (i != row) {
                double factor = A[i][col];
                for (int j = col; j < m; j++)
                    A[i][j] -= factor * A[row][j];
            }

        row++;
        rank++;
    }
    return rank;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

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

    vector<vector<double>> aug(n, vector<double>(n + 1));
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
            fin >> aug[i][j];

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) A[i][j] = aug[i][j];
        b[i] = aug[i][n];
    }

    double det = determinant(A);
    fout << "\nDeterminant = " << det << endl;

    if (fabs(det) > EPS) {
        fout << "\nUnique solution exists."<<endl;

        vector<vector<double>> adj = adjoint(A);
        vector<vector<double>> inv(n, vector<double>(n));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                inv[i][j] = adj[i][j] / det;

        fout << "\nInverse Matrix:"<<endl;
        for (auto &row : inv) {
            for (auto &v : row) fout << setw(10) << v << " ";
            fout << endl;
        }

        vector<double> x(n, 0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                x[i] += inv[i][j] * b[j];

        fout << "\nSolution:"<<endl;
        for (double v : x) fout << v << " ";
        fout << endl;
    }
    else {
        int rankA = rankMatrix(A);
        int rankAug = rankMatrix(aug);

        if (rankA == rankAug)
            fout << "\nInfinite solutions exist."<<endl;
        else
            fout << "\nNo solution exists."<<endl;
    }

    return 0;
}
