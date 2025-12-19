- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Newton-Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss-Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Differential Equation Solving](#differential-equation-solving)
  - [Runge-Kutta 4th Order Method](#runge-kutta-4th-order-method)
    - [Theory](#runge-kutta-theory)
    - [Code](#runge-kutta-code)
    - [Input](#runge-kutta-input)
    - [Output](#runge-kutta-output)

- [Interpolation Methods](#interpolation-methods)
  - [Newton's Forward Interpolation](#newton-forward-interpolation)
    - [Theory](#newton-forward-theory)
    - [Code](#newton-forward-code)
    - [Input](#newton-forward-input)
    - [Output](#newton-forward-output)
  - [Newton's Backward Interpolation](#newton-backward-interpolation)
    - [Theory](#newton-backward-theory)
    - [Code](#newton-backward-code)
    - [Input](#newton-backward-input)
    - [Output](#newton-backward-output)
  - [Newton's Divided Difference Interpolation](#newton-divided-difference-interpolation)
    - [Theory](#newton-divided-difference-theory)
    - [Code](#newton-divided-difference-code)
    - [Input](#newton-divided-difference-input)
    - [Output](#newton-divided-difference-output)

- [Numerical Differentiation](#numerical-differentiation)
  - [Differentiation by Forward Interpolation](#differentiation-forward-interpolation)
    - [Theory](#differentiation-forward-theory)
    - [Code](#differentiation-forward-code)
    - [Input](#differentiation-forward-input)
    - [Output](#differentiation-forward-output)
  - [Differentiation by Backward Interpolation](#differentiation-backward-interpolation)
    - [Theory](#differentiation-backward-theory)
    - [Code](#differentiation-backward-code)
    - [Input](#differentiation-backward-input)
    - [Output](#differentiation-backward-output)

- [Curve Fitting / Regression](#curve-fitting--regression)
  - [Linear Regression](#linear-regression)
    - [Theory](#linear-regression-theory)
    - [Code](#linear-regression-code)
    - [Input](#linear-regression-input)
    - [Output](#linear-regression-output)
  - [Polynomial Regression](#polynomial-regression)
    - [Theory](#polynomial-regression-theory)
    - [Code](#polynomial-regression-code)
    - [Input](#polynomial-regression-input)
    - [Output](#polynomial-regression-output)
  - [Transcendental Regression](#transcendental-regression)
    - [Theory](#transcendental-regression-theory)
    - [Code](#transcendental-regression-code)
    - [Input](#transcendental-regression-input)
    - [Output](#transcendental-regression-output)

- [Numerical Integration](#numerical-integration)
  - [Simpson's 1/3 Rule](#simpsons-13-rule)
    - [Theory](#simpson-13-theory)
    - [Code](#simpson-13-code)
    - [Input](#simpson-13-input)
    - [Output](#simpson-13-output)
  - [Simpson's 3/8 Rule](#simpsons-38-rule)
    - [Theory](#simpson-38-theory)
    - [Code](#simpson-38-code)
    - [Input](#simpson-38-input)
    - [Output](#simpson-38-output)

---

<a id="solution-of-non-linear-equations"></a>
## Solution of Non-Linear Equations

<a id="bisection-method"></a>
### Bisection Method

<a id="bisection-theory"></a>
#### Bisection Theory

```
The bisection method is a numerical technique used to find the root of a non-linear 
equation f(x)=0 by repeatedly dividing an interval [a,b] into two equal parts where f(a) 
and f(b) have opposite signs, and selecting the subinterval that contains the root. 

❖ If a function f(x) is real and continuous in the interval a<x<b, and f(a) and f(b) are of 
opposite sign, that is, f(a)f(b)<0, then there is at least one real root in the interval 
between a and b. 

❖ Let x1=a and x2=b. 

Define x0 to be the midpoint between a and b, that is, x0 = (x1+x2)/2. 

❖ Now there exist the following 3 conditions:​
(i) If f(x0)=0, then the root is x0​
(ii) If f(x0)⋅f(x1)<0, then the root is between x0 and x1 
​(iii) If f(x0)⋅f(x2)<0, then the root is between x0 and x2 

Algorithm:

Step 1: Choose 2 real numbers x1 and x2 such that f(x1)∗f(x2)<0 and stopping criterion 
E.  

Step 2: Define root x0=(x1+x2)/2.  

Step 3: Find f(x0)  

Step 4: If f(x0)=0, then the root is x0 → Stop.  

If f(x0)∗f(x1)<0, then x2=x0  

If f(x0)∗f(x2)<0, then x1=x0  

Return to Step 2 until finding abs(x2−x1)/x2 < E.

```

<a id="bisection-code"></a>
#### Bisection Code

```cpp
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
```

<a id="bisection-input"></a>
#### Bisection Input

```
4
1 -3 2 6 0
3
```

<a id="bisection-output"></a>
#### Bisection Output

```
The polynomial function is: f(x) = 1x^4 -3x^3 +2x^2 +6x^1 

Range found: [-1.1, -1]
Exact root found at x = -1
Range found: [-0.1, 1.52656e-15]
Exact root found at x = 1.52656e-15

Root 1:
No. of iteration: 4
The root is: -1.00312

Root 2:
No. of iteration: 4
The root is: -0.003125
```

---

<a id="false-position-method"></a>
### False Position Method

<a id="false-position-theory"></a>
#### False Position Theory

``` 

The false position method (Regular Falsi method in Latin) finds the root of a non-linear 
equation f(x)=0 by repeatedly approximating the root using a straight line joining the 
points (a,f(a)) and (b,f(b)), where f(a) and f(b) have opposite signs. The point where this 
line intersects the x-axis gives a better approximation of the root. 

❖ If a function f(x) is real and continuous in the interval a<x<b, and f(a) and f(b) are of 
opposite sign, that is, f(a)f(b)<0, then there is at least one real root in the interval 
between a and b. 

❖ Let x1=a and x2=b. Define x0 to be the midpoint between a and b, that is,  

x0 = x1 - f(x1)*(x2-x1)/(f(x2)-f(x1)). 

❖ Now there exist the following 3 conditions:​
(i) If f(x0)=0, then the root is x0​
(ii) If f(x0)⋅f(x1)<0, then the root is between x0 and x1 
​(iii) If f(x0)⋅f(x2)<0, then the root is between x0 and x2 

Algorithm:

Step 1: Choose 2 real numbers x1 and x2 such that f(x1)∗f(x2)<0 and stopping criterion 
E.  

Step 2: Define root x0 = x1 - f(x1)*(x2-x1)/(f(x2)-f(x1)).  

Step 3: Find f(x0)  

Step 4: If f(x0)=0, then the root is x0 → Stop.  

If f(x0)∗f(x1)<0, then x2=x0  

If f(x0)∗f(x2)<0, then x1=x0  

Return to Step 2 until finding abs(x2−x1)/x2 < E.

```

<a id="false-position-code"></a>
#### False Position Code

```cpp
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

void false_position(double a, double b) {
    int itr = 0;
    if (func(a) * func(b) >= 0) {
        fout << "No root found in range [" << a << ", " << b << "]" << endl;
        return;
    }

    double c;

    while ((b - a) >= EP) {
        c = (a * func(b) - b * func(a)) / (func(b) - func(a));
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
        false_position(ranges[i].first, ranges[i].second);
    }

    if (ranges.empty()) {
        fout << "No roots found in the given range." << endl;
    }

    return 0;
}
```

<a id="false-position-input"></a>
#### False Position Input

```
4
1 -3 2 6 0
50
```

<a id="false-position-output"></a>
#### False Position Output

```
The polynomial function is: f(x) = 1x^4 -3x^3 +2x^2 +6x^1 

Range found: [-1.1, -1]
Exact root found at x = -1
Range found: [-0.1, 4.41619e-13]
Exact root found at x = 4.41619e-13

Root 1:
No. of iteration: 5
The root is: -1

Root 2:
No. of iteration: 1
The root is: -1.76831e-14
```

---

<a id="newton-raphson-method"></a>
### Newton-Raphson Method

<a id="newton-raphson-theory"></a>
#### Newton-Raphson Theory

```
Theory Here
```

<a id="newton-raphson-code"></a>
#### Newton-Raphson Code

```cpp
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
```

<a id="newton-raphson-input"></a>
#### Newton-Raphson Input

```
2
1 0 -4
-1 3
```

<a id="newton-raphson-output"></a>
#### Newton-Raphson Output

```
The polynomial function is: f(x) = 1x^2 -4

Root 1:
The root is: -2
Iteration needed: 4

Root 2:
The root is: 2
Iteration needed: 4
```

---

<a id="secant-method"></a>
### Secant Method

<a id="secant-theory"></a>
#### Secant Theory

```
The Secant Method is a numerical technique used to find the root of a non-linear equation
f(x) = 0 by repeatedly approximating the root using a straight line passing through two
points on the curve.

❖ Unlike the Bisection and False Position methods, the Secant Method does not require
f(x₀) and f(x₁) to have opposite signs.

❖ Let x₀ and x₁ be two initial approximations of the root.

❖ The next approximation x₂ is obtained by:

x₂ = x₁ − f(x₁)(x₁ − x₀) / [f(x₁) − f(x₀)]

❖ The process is repeated until the desired accuracy is achieved.

Algorithm:
Step 1: Choose two initial approximations x₀ and x₁.  
Step 2: Evaluate f(x₀) and f(x₁).  
Step 3: Compute the next approximation using the secant formula.  
Step 4: Replace x₀ by x₁ and x₁ by the new approximation.  
Step 5: Repeat Steps 2–4 until |x₁ − x₀| is less than the prescribed tolerance.  
Step 6: The final value of x₁ is taken as the approximate root.

```

<a id="secant-code"></a>
#### Secant Code

```cpp
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
```

<a id="secant-input"></a>
#### Secant Input

```
2
1 0 -4
-3 -1 1 3
```

<a id="secant-output"></a>
#### Secant Output

```
The polynomial function is: f(x) = 1x^2 -4

Initial guess pairs: (-3, -1), (1, 3)

Root 1:
The root is: -1.99987
Iteration needed: 4

Root 2:
The root is: 1.99995
Iteration needed: 4
```

---

<a id="solution-of-linear-equations"></a>
## Solution of Linear Equations

<a id="gauss-elimination-method"></a>
### Gauss Elimination Method

<a id="gauss-elimination-theory"></a>
#### Gauss Elimination Theory

```
The Gauss Elimination method solves a system of linear equations by systematically 
eliminating variables to convert the coefficient matrix into an upper triangular form, 
and then finding the solution using backward substitution. 

Algorithm:

Step 1: Take the given system of linear equations and write it in matrix form: AX=B  

Step 2: Form the augmented matrix [A|B].  

Step 3: Use elementary row operations to eliminate the unknowns below the main 
diagonal and convert the augmented matrix into an upper triangular form.  

Step 4: Continue the elimination process until all elements below the main diagonal 
become zero.  

Step 5: After obtaining the upper triangular matrix, solve the resulting system using 
backward substitution.  

Step 6: Compute the values of the unknown variables.  

Step 7: Obtain the solution vector X.

```

<a id="gauss-elimination-code"></a>
#### Gauss Elimination Code

```cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

ofstream fout("output.txt");

void printMatrix(double a[][20], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            fout << setw(10) << fixed << setprecision(4) << a[i][j] << " ";
        }
        fout << endl;
    }
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

    int n;
    fin >> n;

    double a[20][20], x[20];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            fin >> a[i][j];
        }
    }

    fout << "\nInitial Matrix:"<<endl;
    printMatrix(a, n);

    for (int i = 0; i < n - 1; i++) {

        for (int k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[i][i])) {
                for (int j = 0; j <= n; j++) {
                    swap(a[i][j], a[k][j]);
                }
                fout << "\nAfter swapping row " << i + 1 << " and row " << k + 1 << ":"<<endl;
                printMatrix(a, n);
            }
        }

        if (fabs(a[i][i]) < 1e-9) continue;

        for (int k = i + 1; k < n; k++) {
            double factor = a[k][i] / a[i][i];
            fout << "\nR" << k + 1 << " = R" << k + 1
                 << " - (" << factor << ") * R" << i + 1 << endl;

            for (int j = i; j <= n; j++) {
                a[k][j] -= factor * a[i][j];
            }
            printMatrix(a, n);
        }
    }

    fout << "\nMatrix after Forward Elimination (Upper Triangular):"<< endl;
    printMatrix(a, n);

    bool noSolution = false, infinite = false;

    for (int i = 0; i < n; i++) {
        bool allZero = true;
        for (int j = 0; j < n; j++) {
            if (fabs(a[i][j]) > 1e-9) {
                allZero = false;
                break;
            }
        }

        if (allZero && fabs(a[i][n]) > 1e-9) {
            noSolution = true;
            break;
        }
        else if (allZero && fabs(a[i][n]) < 1e-9) {
            infinite = true;
        }
    }

    if (noSolution) {
        fout << "\nThe system has No Solution (Inconsistent system)." << endl;
        return 0;
    }
    else if (infinite) {
        fout << "\nThe system has Infinite Solutions (Dependent equations)." << endl;
        return 0;
    }

    fout << "\nBack Substitution Steps:"<< endl;

    for (int i = n - 1; i >= 0; i--) {
        double sum = a[i][n];
        fout << "x" << i + 1 << " = ( " << a[i][n];

        for (int j = i + 1; j < n; j++) {
            sum -= a[i][j] * x[j];
            fout << " - (" << a[i][j] << " * x" << j + 1 << ")";
        }

        x[i] = sum / a[i][i];
        fout << " ) / " << a[i][i] << endl;
        fout << "x" << i + 1 << " = " << fixed << setprecision(6) << x[i] << endl << endl;
    }

    fout << "\nUnique Solution:"<<endl;
    for (int i = 0; i < n; i++) {
        fout << "x" << i + 1 << " = " << fixed << setprecision(6) << x[i] << endl;
    }

    return 0;
}
```

<a id="gauss-elimination-input"></a>
#### Gauss Elimination Input

```
3
2 1 1 10
3 2 3 18
1 4 9 16
```

<a id="gauss-elimination-output"></a>
#### Gauss Elimination Output

```

Initial Matrix:
    2.0000     1.0000     1.0000    10.0000 
    3.0000     2.0000     3.0000    18.0000 
    1.0000     4.0000     9.0000    16.0000 

After swapping row 1 and row 2:
    3.0000     2.0000     3.0000    18.0000 
    2.0000     1.0000     1.0000    10.0000 
    1.0000     4.0000     9.0000    16.0000 

R2 = R2 - (0.6667) * R1
    3.0000     2.0000     3.0000    18.0000 
    0.0000    -0.3333    -1.0000    -2.0000 
    1.0000     4.0000     9.0000    16.0000 

R3 = R3 - (0.3333) * R1
    3.0000     2.0000     3.0000    18.0000 
    0.0000    -0.3333    -1.0000    -2.0000 
    0.0000     3.3333     8.0000    10.0000 

After swapping row 2 and row 3:
    3.0000     2.0000     3.0000    18.0000 
    0.0000     3.3333     8.0000    10.0000 
    0.0000    -0.3333    -1.0000    -2.0000 

R3 = R3 - (-0.1000) * R2
    3.0000     2.0000     3.0000    18.0000 
    0.0000     3.3333     8.0000    10.0000 
    0.0000     0.0000    -0.2000    -1.0000 

Matrix after Forward Elimination (Upper Triangular):
    3.0000     2.0000     3.0000    18.0000 
    0.0000     3.3333     8.0000    10.0000 
    0.0000     0.0000    -0.2000    -1.0000 

Back Substitution Steps:
x3 = ( -1.0000 ) / -0.2000
x3 = 5.000000

x2 = ( 10.000000 - (8.000000 * x3) ) / 3.333333
x2 = -9.000000

x1 = ( 18.000000 - (2.000000 * x2) - (3.000000 * x3) ) / 3.000000
x1 = 7.000000


Unique Solution:
x1 = 7.000000
x2 = -9.000000
x3 = 5.000000
```

---

<a id="gauss-jordan-elimination-method"></a>
### Gauss-Jordan Elimination Method

<a id="gauss-jordan-theory"></a>
#### Gauss-Jordan Theory

```

The Gauss–Jordan elimination method is a numerical technique used to solve a system 
of linear equations by transforming the augmented matrix into a reduced row echelon 
form (identity matrix) using elementary row operations.  

Algorithm:

Step 1: Take the given system of linear equations and write it in matrix form: AX=B  

Step 2: Form the augmented matrix [A|B].  

Step 3: Use elementary row operations to make the leading diagonal elements equal 
to 1.  

Step 4: Use row operations to make all other elements in each column zero, both above 
and below the leading diagonal.  

Step 5: Continue the process until the coefficient matrix becomes an identity matrix.  

Step 6: Read the solution directly from the augmented matrix.

```

<a id="gauss-jordan-code"></a>
#### Gauss-Jordan Code

```cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

ofstream fout("output.txt");

void printMatrix(double a[][20], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            fout << setw(10) << fixed << setprecision(4) << a[i][j] << " ";
        }
        fout << endl;
    }
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

    int n;
    fin >> n;

    double a[20][20];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            fin >> a[i][j];
        }
    }

    fout << "\nInitial Matrix:"<<endl;
    printMatrix(a, n);

    for (int i = 0; i < n; i++) {
        for (int k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[i][i])) {
                for (int j = 0; j <= n; j++) {
                    swap(a[i][j], a[k][j]);
                }
                fout << "\nAfter swapping row " << i + 1 << " and row " << k + 1 << ":"<<endl;
                printMatrix(a, n);
            }
        }

        double pivot = a[i][i];
        if (fabs(pivot) < 1e-9) continue;

        for (int j = 0; j <= n; j++) {
            a[i][j] /= pivot;
        }
        fout << "\nAfter dividing row " << i + 1 << " by pivot (" << pivot << "):"<<endl;
        printMatrix(a, n);

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = a[k][i];
                for (int j = 0; j <= n; j++) {
                    a[k][j] -= factor * a[i][j];
                }
            }
        }
        fout << "\nAfter eliminating column " << i + 1 << ":"<<endl;
        printMatrix(a, n);
    }

    fout << "\nFinal Matrix (Reduced Row Echelon Form):"<<endl;
    printMatrix(a, n);

    bool noSolution = false, infinite = false;
    for (int i = 0; i < n; i++) {
        bool allZero = true;
        for (int j = 0; j < n; j++) {
            if (fabs(a[i][j]) > 1e-9) {
                allZero = false;
                break;
            }
        }
        if (allZero && fabs(a[i][n]) > 1e-9) {
            noSolution = true;
            break;
        }
        else if (allZero && fabs(a[i][n]) < 1e-9) {
            infinite = true;
        }
    }

    if (noSolution) {
        fout << "\nThe system has No Solution (Inconsistent system)." << endl;
    }
    else if (infinite) {
        fout << "\nThe system has Infinite Solutions (Dependent equations)." << endl;
    }
    else {
        fout << "\nThe system has a Unique Solution:"<<endl;
        for (int i = 0; i < n; i++) {
            fout << "x" << i + 1 << " = " << fixed << setprecision(6) << a[i][n] << endl;
        }
    }

    return 0;
}
```

<a id="gauss-jordan-input"></a>
#### Gauss-Jordan Input

```
3
2 1 1 10
3 2 3 18
1 4 9 16
```

<a id="gauss-jordan-output"></a>
#### Gauss-Jordan Output

```

Initial Matrix:
    2.0000     1.0000     1.0000    10.0000 
    3.0000     2.0000     3.0000    18.0000 
    1.0000     4.0000     9.0000    16.0000 

After swapping row 1 and row 2:
    3.0000     2.0000     3.0000    18.0000 
    2.0000     1.0000     1.0000    10.0000 
    1.0000     4.0000     9.0000    16.0000 

After dividing row 1 by pivot (3.0000):
    1.0000     0.6667     1.0000     6.0000 
    2.0000     1.0000     1.0000    10.0000 
    1.0000     4.0000     9.0000    16.0000 

After eliminating column 1:
    1.0000     0.6667     1.0000     6.0000 
    0.0000    -0.3333    -1.0000    -2.0000 
    0.0000     3.3333     8.0000    10.0000 

After swapping row 2 and row 3:
    1.0000     0.6667     1.0000     6.0000 
    0.0000     3.3333     8.0000    10.0000 
    0.0000    -0.3333    -1.0000    -2.0000 

After dividing row 2 by pivot (3.3333):
    1.0000     0.6667     1.0000     6.0000 
    0.0000     1.0000     2.4000     3.0000 
    0.0000    -0.3333    -1.0000    -2.0000 

After eliminating column 2:
    1.0000     0.0000    -0.6000     4.0000 
    0.0000     1.0000     2.4000     3.0000 
    0.0000     0.0000    -0.2000    -1.0000 

After dividing row 3 by pivot (-0.2000):
    1.0000     0.0000    -0.6000     4.0000 
    0.0000     1.0000     2.4000     3.0000 
   -0.0000    -0.0000     1.0000     5.0000 

After eliminating column 3:
    1.0000     0.0000     0.0000     7.0000 
    0.0000     1.0000     0.0000    -9.0000 
   -0.0000    -0.0000     1.0000     5.0000 

Final Matrix (Reduced Row Echelon Form):
    1.0000     0.0000     0.0000     7.0000 
    0.0000     1.0000     0.0000    -9.0000 
   -0.0000    -0.0000     1.0000     5.0000 

The system has a Unique Solution:
x1 = 7.000000
x2 = -9.000000
x3 = 5.000000
```

---

<a id="lu-decomposition-method"></a>
### LU Decomposition Method

<a id="lu-decomposition-theory"></a>
#### LU Decomposition Theory

```
The LU decomposition method is a numerical technique used to solve a system of 
linear equations by factorizing a square matrix A into the product of a lower triangular 
matrix L and an upper triangular matrix U, such that A=LU. 

Algorithm:

Step 1: Take the given system of linear equations and write it in matrix form: AX=B  

Step 2: Form the coefficient matrix A.  

Step 3: Decompose matrix A into: A=LU  

where  

●​ L is a lower triangular matrix with diagonal elements equal to 1  

●​ U is an upper triangular matrix​  

Step 4: Compute:  

First row of U : u11 ,u12,u13  

First column of L : l21 ,l32  

Second row of U : u22,u23  

Second column of L : l23  

Third row of U : u33  

Step 5: After finding matrices L and U, solve: LY=B by forward substitution.  

Step 6: Solve: UX=Y by backward substitution.Obtain the solution vector XXX.  

Step 7: Obtain the solution vector X.

```

<a id="lu-decomposition-code"></a>
#### LU Decomposition Code

```cpp
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
```

<a id="lu-decomposition-input"></a>
#### LU Decomposition Input

```
3
1 5 1 14
2 1 3 13
3 1 4 17
```

<a id="lu-decomposition-output"></a>
#### LU Decomposition Output

```

After Computing U[1][1]:
Matrix L:
    1.0000     0.0000     0.0000 
    0.0000     1.0000     0.0000 
    0.0000     0.0000     1.0000 

Matrix U:
    1.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000 

After Computing U[1][2]:
Matrix L:
    1.0000     0.0000     0.0000 
    0.0000     1.0000     0.0000 
    0.0000     0.0000     1.0000 

Matrix U:
    1.0000     5.0000     0.0000 
    0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000 

After Computing U[1][3]:
Matrix L:
    1.0000     0.0000     0.0000 
    0.0000     1.0000     0.0000 
    0.0000     0.0000     1.0000 

Matrix U:
    1.0000     5.0000     1.0000 
    0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000 

After Computing L[2][1]:
Matrix L:
    1.0000     0.0000     0.0000 
    2.0000     1.0000     0.0000 
    0.0000     0.0000     1.0000 

Matrix U:
    1.0000     5.0000     1.0000 
    0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000 

After Computing L[3][1]:
Matrix L:
    1.0000     0.0000     0.0000 
    2.0000     1.0000     0.0000 
    3.0000     0.0000     1.0000 

Matrix U:
    1.0000     5.0000     1.0000 
    0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000 

After Computing U[2][2]:
Matrix L:
    1.0000     0.0000     0.0000 
    2.0000     1.0000     0.0000 
    3.0000     0.0000     1.0000 

Matrix U:
    1.0000     5.0000     1.0000 
    0.0000    -9.0000     0.0000 
    0.0000     0.0000     0.0000 

After Computing U[2][3]:
Matrix L:
    1.0000     0.0000     0.0000 
    2.0000     1.0000     0.0000 
    3.0000     0.0000     1.0000 

Matrix U:
    1.0000     5.0000     1.0000 
    0.0000    -9.0000     1.0000 
    0.0000     0.0000     0.0000 

After Computing L[3][2]:
Matrix L:
    1.0000     0.0000     0.0000 
    2.0000     1.0000     0.0000 
    3.0000     1.5556     1.0000 

Matrix U:
    1.0000     5.0000     1.0000 
    0.0000    -9.0000     1.0000 
    0.0000     0.0000     0.0000 

After Computing U[3][3]:
Matrix L:
    1.0000     0.0000     0.0000 
    2.0000     1.0000     0.0000 
    3.0000     1.5556     1.0000 

Matrix U:
    1.0000     5.0000     1.0000 
    0.0000    -9.0000     1.0000 
    0.0000     0.0000    -0.5556 

Final L and U matrices:
Matrix L:
    1.0000     0.0000     0.0000 
    2.0000     1.0000     0.0000 
    3.0000     1.5556     1.0000 

Matrix U:
    1.0000     5.0000     1.0000 
    0.0000    -9.0000     1.0000 
    0.0000     0.0000    -0.5556 

Forward Substitution (Ly = b)
y1 = 14.000000
y2 = -15.000000
y3 = -1.666667

Backward Substitution (Ux = y)
x3 = 3.000000
x2 = 2.000000
x1 = 1.000000

The system has a Unique Solution.
```

---

<a id="matrix-inversion"></a>
### Matrix Inversion

<a id="matrix-inversion-theory"></a>
#### Matrix Inversion Theory

```  

The matrix inversion method solves a system of linear equations AX=B by finding the 
inverse of the coefficient matrix A. If A is non-singular (det⁡(A)≠0), the solution is 
obtained as: X=(A^-1)B 

Algorithm:

Step 1: Take the given system of linear equations and write it in matrix form: AX=B.  

Step 2: Form the coefficient matrix A and the constant matrix B.  

Step 3: Check whether matrix A is invertible (i.e., det⁡(A)≠0).  

Step 4: Find the inverse of matrix A, denoted by (A^-1), using any suitable method (adjoint  
method or elementary row operations).  

Step 5: Multiply both sides of the equation AX=B by  (A^-1): X=(A^-1)B.  

Step 6: Perform the matrix multiplication (A^-1)B.  

Step 7: Obtain the solution vector X.

```

<a id="matrix-inversion-code"></a>
#### Matrix Inversion Code

```cpp
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
```

<a id="matrix-inversion-input"></a>
#### Matrix Inversion Input

```
3
1 1 1 6
3 3 4 20
2 1 3 13
```

<a id="matrix-inversion-output"></a>
#### Matrix Inversion Output

```

Determinant = 1

Unique solution exists.

Inverse Matrix:
         5         -2          1 
        -1          1         -1 
        -3          1          0 

Solution:
3 1 2 
```

---

<a id="differential-equation-solving"></a>
## Differential Equation Solving

<a id="runge-kutta-4th-order-method"></a>
### Runge-Kutta 4th Order Method

<a id="runge-kutta-theory"></a>
#### Runge-Kutta Theory

```
Theory Here
```

<a id="runge-kutta-code"></a>
#### Runge-Kutta Code

```cpp
#include <bits/stdc++.h>
using namespace std;

// Differential equation: dy/dx = x^2 - y
double f(double x, double y) {
    return x*x - y;
}

double rungeKutta(double x0, double y0, double h, double x) {
    int n = (x - x0) / h;
    double k1, k2, k3, k4;

    for (int i = 0; i < n; i++) {
        k1 = h * f(x0, y0);
        k2 = h * f(x0 + h/2, y0 + k1/2);
        k3 = h * f(x0 + h/2, y0 + k2/2);
        k4 = h * f(x0 + h, y0 + k3);

        y0 = y0 + (k1 + 2*k2 + 2*k3 + k4) / 6;
        x0 = x0 + h;
    }

    return y0;
}

int main() {
    ifstream inFile("input.txt");
    ofstream outFile("output.txt");

    if (!inFile) {
        cout << "Error: Could not open Input.txt" << endl;
        return 1;
    }
    int T;
    inFile >> T;
     for (int tc = 1; tc <= T; tc++) {

    double x0, y0, h, x;

    inFile >> x0 >> y0 >> h >> x;

    double result = rungeKutta(x0, y0, h, x);

    cout << fixed << setprecision(4);
    cout << "Result of y at x = " << x << "is : " << result << endl;

    outFile << fixed << setprecision(4);
    outFile << "Result of y at x = " << x << "is : " << result << endl;
     }

    inFile.close();
    outFile.close();

    return 0;
}
```

<a id="runge-kutta-input"></a>
#### Runge-Kutta Input

```
3

0
1
0.1
0.5

1
2
0.05
1.5

2
3
0.1
2.5
```

<a id="runge-kutta-output"></a>
#### Runge-Kutta Output

```
Result of y at x = 0.500000is : 0.689680
Result of y at x = 1.500000is : 1.840128
Result of y at x = 2.500000is : 3.630321
```

---

<a id="interpolation-methods"></a>
## Interpolation Methods

<a id="newton-forward-interpolation"></a>
### Newton’s Forward Interpolation

<a id="newton-forward-theory"></a>
#### Newton’s Forward Interpolation Theory

```
Newton’s Forward Interpolation is a numerical technique used to estimate the value of a function at a given point using equally spaced tabulated data, when the value lies near the beginning of the data table.

❖ The independent variable values must be equally spaced.  
❖ Let the given data points be  
(x₀, y₀), (x₁, y₁), … , (xₙ, yₙ) where xᵢ = x₀ + ih.  
❖ Forward difference operator:  
Δyᵢ = yᵢ₊₁ − yᵢ  

❖ Interpolation formula:  
f(x) = y₀ + uΔy₀ + u(u−1)/2! Δ²y₀ + u(u−1)(u−2)/3! Δ³y₀ + …

where  
u = (x − x₀) / h

Algorithm:

Step 1: Arrange the given data in tabular form.  
Step 2: Construct the forward difference table.  
Step 3: Identify x₀, the first value of x, and compute h.  
Step 4: Compute u = (x − x₀)/h.  
Step 5: Substitute the values of u and the forward differences into the interpolation formula.  
Step 6: Evaluate the required value of f(x).
```



<a id="newton-forward-code"></a>
#### Newton’s Forward Interpolation Code

```cpp
#include <bits/stdc++.h>
using namespace std;

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

    
    fout << "Forward Difference Table:" << endl << endl;
    for (int i = 0; i < n; i++) {
        fout << setw(8) << x[i] << "\t";
        for (int j = 0; j < n - i; j++)
            fout << setw(10) << y[i][j] << "\t";
        fout << endl;
    }

    
    float sum = y[0][0];
    float u = (value - x[0]) / (x[1] - x[0]);
    
    float u_term = u;
    for (int i = 1; i < n; i++) {
        sum = sum + (u_term * y[0][i]) / factorial(i);
        u_term = u_term * (u - i);
    }

    fout << "\nInterpolated value at X = " << value << " is Y = " << sum << endl;

    fin.close();
    fout.close();

    return 0;
}
```

<a id="newton-forward-input"></a>
#### Newton’s Forward Interpolation Input

```
5
1891 46
1901 66
1911 81
1921 93
1931 101
1895
```

<a id="newton-forward-output"></a>
#### Newton’s Forward Interpolation Output

```
Forward Difference Table:

    1891	        46	        20	        -5	         2	        -3	
    1901	        66	        15	        -3	        -1	
    1911	        81	        12	        -4	
    1921	        93	         8	
    1931	       101	

Interpolated value at X = 1895 is Y = 54.8528
```

---

<a id="newton-backward-interpolation"></a>
### Newton's Backward Interpolation

<a id="newton-backward-theory"></a>
#### Newton's Backward Interpolation Theory
```
Newton’s Backward Interpolation is used to estimate the value of a function when the required value lies near the end of the data table.

❖ The independent variable values must be equally spaced.  
❖ Let xₙ be the last value of x.  
❖ Backward difference operator:  
∇yₙ = yₙ − yₙ₋₁  

❖ Interpolation formula:  
f(x) = yₙ + v∇yₙ + v(v+1)/2! ∇²yₙ + v(v+1)(v+2)/3! ∇³yₙ + …

where  
v = (x − xₙ) / h

Algorithm: 

Step 1: Arrange the given data in tabular form.  
Step 2: Construct the backward difference table.  
Step 3: Identify xₙ and compute h.  
Step 4: Compute v = (x − xₙ)/h.  
Step 5: Substitute values into the backward interpolation formula.  
Step 6: Evaluate the required value of f(x).

```

<a id="newton-backward-code"></a>
#### Newton's Backward Interpolation Code

```cpp
#include <bits/stdc++.h>
using namespace std;

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

    fout << "Backward Difference Table:" << endl << endl;

    for (int i = 0; i < n; i++) {
        fout << setw(8) << x[i] << "\t";
        for (int j = 0; j <= i; j++)
            fout << setw(10) << y[i][j] << "\t";
        fout << endl;
    }

    double sum = y[n - 1][0];
    double v = (value - x[n - 1]) / (x[1] - x[0]);
    
    double v_term = v;
    for (int i = 1; i < n; i++) {
        sum = sum + (v_term * y[n - 1][i]) / factorial(i);
        v_term = v_term * (v + i);
    }

    fout << "\nInterpolated value at X = " << value << " is Y = " << sum << endl;

    fin.close();
    fout.close();

    return 0;
}
```

<a id="newton-backward-input"></a>
#### Newton's Backward Interpolation Input

```
5
1891 46
1901 66
1911 81
1921 93
1931 101
1925
```

<a id="newton-backward-output"></a>
#### Newton's Backward Interpolation Output

```
Forward Difference Table:

    1891	        46	        20	        -5	         2	        -3	
    1901	        66	        15	        -3	        -1	
    1911	        81	        12	        -4	
    1921	        93	         8	
    1931	       101	

Interpolated value at X = 1895 is Y = 54.8528
```

---

<a id="newton-divided-difference-interpolation"></a>
### Newton’s Divided Difference Interpolation

<a id="newton-divided-difference-theory"></a>
#### Newton’s Divided Difference Theory

```
Newton’s Divided Difference Interpolation is used to estimate function values when the data points are not equally spaced.

❖ First divided difference:  
f[xᵢ, xⱼ] = (f(xⱼ) − f(xᵢ)) / (xⱼ − xᵢ)

❖ Interpolation polynomial:  
f(x) = f(x₀) + (x − x₀)f[x₀, x₁] + (x − x₀)(x − x₁)f[x₀, x₁, x₂] + …


Algorithm:

Step 1: Arrange the data points in ascending order of x.  
Step 2: Construct the divided difference table.  
Step 3: Select the required divided differences.  
Step 4: Substitute values into the divided difference polynomial.  
Step 5: Evaluate f(x) at the desired point.

```
<a id="newton-divided-difference-code"></a>
#### Newton’s Divided Difference Code

```cpp
#include <bits/stdc++.h>
using namespace std;

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
        for (int j = 0; j < n - i; j++) {
            y[j][i] = (y[j + 1][i - 1] - y[j][i - 1]) / (x[j + i] - x[j]);
        }
    }

    fout << "Divided Difference Table:" << endl << endl;
    for (int i = 0; i < n; i++) {
        fout << setw(8) << x[i] << "\t";
        for (int j = 0; j < n - i; j++)
            fout << setw(10) << y[i][j] << "\t";
        fout << endl;
    }

    float sum = y[0][0];
    float product = 1;
    
    for (int i = 1; i < n; i++) {
        product = product * (value - x[i - 1]);
        sum = sum + (product * y[0][i]);
    }

    fout << "\nInterpolated value at X = " << value << " is Y = " << sum << endl;

    fin.close();
    fout.close();

    return 0;
}
```

<a id="newton-divided-difference-input"></a>
#### Newton’s Divided Difference Input

```
5
1891 46
1901 66
1911 81
1921 93
1931 101
1895
```

<a id="newton-divided-difference-output"></a>
#### Newton’s Divided Difference Output

```
Divided Difference Table:

    1891	        46	         2	    -0.025	0.000333333	 -1.25e-05	
    1901	        66	       1.5	    -0.015	-0.000166667	
    1911	        81	       1.2	     -0.02	
    1921	        93	       0.8	
    1931	       101	

Interpolated value at X = 1895 is Y = 54.8528
```

---

<a id="numerical-differentiation"></a>
## Numerical Differentiation

<a id="differentiation-forward-interpolation"></a>
### Differentiation by Forward Interpolation

<a id="differentiation-forward-theory"></a>
#### Differentiation by Forward Interpolation Theory

```
Newton’s Forward Differentiation is used to approximate the derivative of a function using equally spaced tabulated data, when the point of differentiation lies near the beginning of the data table.

❖ The independent variable values must be equally spaced.  
❖ Forward differences are used.  
❖ Let h be the common interval and  
u = (x − x₀) / h  

❖ The first derivative obtained from Newton’s forward interpolation polynomial is:

(dy/dx) = [ Δy₀ + (2u − 1)/2 · Δ²y₀ + (3u² − 6u + 2)/6 · Δ³y₀ + (4u³ − 18u² + 22u − 6)/24 · Δ⁴y₀ + … ] / h

Algorithm:

Step 1: Arrange the given data in tabular form.  
Step 2: Construct the forward difference table.  
Step 3: Identify the initial value x₀ and compute the step size h.  
Step 4: Compute u = (x − x₀)/h.  
Step 5: Substitute u and the forward differences into the differentiation formula.  
Step 6: Evaluate the approximate value of dy/dx.

```

<a id="differentiation-forward-code"></a>
#### Differentiation by Forward Interpolation Code

```cpp
#include <bits/stdc++.h>
using namespace std;

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

    float h = x[1] - x[0];

    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n - i; j++)
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
    }

    int i = 0;
    float u = (value - x[i]) / h;
    
    float derivative = y[i][1] / h;

    if (n >= 3) {
        derivative += ((2 * u - 1) / 2.0) * (y[i][2] / h);
    }

    if (n >= 4) {
        derivative += ((3 * u * u - 6 * u + 2) / 6.0) * (y[i][3] / h);
    }

    if (n >= 5) {
        derivative += ((4 * u * u * u - 18 * u * u + 22 * u - 6) / 24.0) * (y[i][4] / h);
    }

    fout << "\nAt x = " << value << ", Derivative = " << derivative << endl;

    fin.close();
    fout.close();

    return 0;
}
```

<a id="differentiation-forward-input"></a>
#### Differentiation Forward Input

```
5
0 1
1 2.7183
2 7.3891
3 20.0855
4 54.5982
1.2
```

<a id="differentiation-forward-output"></a>
#### Differentiation Forward Output

```
At x = 1.2, Derivative = 3.54662
```

---

<a id="differentiation-backward-interpolation"></a>
### Differentiation by Backward Interpolation

<a id="differentiation-backward-theory"></a>
#### Differentiation by Backward Interpolation Theory

```
Newton’s Backward Differentiation is used to approximate the derivative of a function using equally spaced tabulated data, when the point of differentiation lies near the end of the data table.

❖ The independent variable values must be equally spaced.  
❖ Backward differences are used.  
❖ Let h be the common interval and  
v = (x − xₙ) / h  

❖ The first derivative obtained from Newton’s backward interpolation polynomial is:

(dy/dx) = [ ∇yₙ + (2v + 1)/2 · ∇²yₙ + (3v² + 6v + 2)/6 · ∇³yₙ + (4v³ + 18v² + 22v + 6)/24 · ∇⁴yₙ + … ] / h

Algorithm:

Step 1: Arrange the given data in tabular form.  
Step 2: Construct the backward difference table.  
Step 3: Identify the final value xₙ and compute the step size h.  
Step 4: Compute v = (x − xₙ)/h.  
Step 5: Substitute v and the backward differences into the differentiation formula.  
Step 6: Evaluate the approximate value of dy/dx.

```

<a id="differentiation-backward-code"></a>
#### Differentiation by Backward Interpolation Code

```cpp
#include <bits/stdc++.h>
using namespace std;

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

    float h = x[1] - x[0];

    for (int i = 1; i < n; i++) {
        for (int j = n - 1; j >= i; j--)
            y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
    }

    int i = n - 1;
    float v = (value - x[i]) / h;

    float derivative = y[i][1] / h;
    
    if (n >= 3) {
        derivative += ((2 * v + 1) / 2.0) * (y[i][2] / h);
    }
    
    if (n >= 4) {
        derivative += ((3 * v * v + 6 * v + 2) / 6.0) * (y[i][3] / h);
    }
    
    if (n >= 5) {
        derivative += ((4 * v * v * v + 18 * v * v + 22 * v + 6) / 24.0) * (y[i][4] / h);
    }

    fout << "\nAt x = " << value << ", Derivative = " << derivative << endl;

    fin.close();
    fout.close();

    return 0;
}
```

<a id="differentiation-backward-input"></a>
#### Differentiation by Backward Interpolation Input

```
5
0 1
1 2.7183
2 7.3891
3 20.0855
4 54.5982
3.4
```

<a id="differentiation-backward-output"></a>
#### Differentiation by Backward Interpolation Output

```
At x = 3.4, Derivative = 30.5605
```

---

<a id="curve-fitting--regression"></a>
## Curve Fitting / Regression

<a id="linear-regression"></a>
### Linear Regression

<a id="linear-regression-theory"></a>
#### Linear Regression Theory

```
Curve fitting is a numerical method used to determine a mathematical relationship that
best represents a given set of experimental data points.

❖ In linear regression, the data is assumed to follow a straight-line relationship of
the form:

y = a + bx

❖ The constants a and b are determined using the Principle of Least Squares, which states
that the sum of the squares of the deviations between the observed and calculated values
is minimum.

❖ The normal equations obtained are:

Σy = na + bΣx  
Σxy = aΣx + bΣx²

❖ Solving the above equations gives the explicit expressions for a and b:

b = [ nΣxy − (Σx)(Σy) ] / [ nΣx² − (Σx)² ]  
a = ( Σy − bΣx ) / n

❖ Substituting the values of a and b gives the best-fit straight line.

Algorithm:
Step 1: Arrange the given data points (x, y) in tabular form.  
Step 2: Compute x² and xy for each data point.  
Step 3: Calculate Σx, Σy, Σx², and Σxy.  
Step 4: Form the normal equations.  
Step 5: Solve the equations to determine a and b.  
Step 6: Write the equation of the best-fit straight line y = a + bx.

```

<a id="linear-regression-code"></a>
#### Linear Regression Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double predict(double a, double b, double x) {
    return a + b * x;
}

bool computeRegression(
    const vector<double>& x,
    const vector<double>& y,
    double& a,
    double& b
) {
    int n = x.size();
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    double denominator = (n * sumX2) - (sumX * sumX);
    if (denominator == 0) return false;

    b = (n * sumXY - sumX * sumY) / denominator;
    a = (sumY - b * sumX) / n;

    return true;
}

int main() {
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
int n, testCase = 1;

while (true) {
    if (!(fin >> n)) break;

    vector<double> x(n), y(n);

    for (int i = 0; i < n; i++) {
        fin >> x[i] >> y[i];
    }

    double X;
    fin >> X;

    cout << "Test Case " << testCase << "\n";
    fout << "Test Case " << testCase << "\n";

    double a, b;
    if (!computeRegression(x, y, a, b)) {
        cout << "Error: Regression cannot be computed.\n\n";
        fout << "Error: Regression cannot be computed.\n\n";
        testCase++;
        continue;
    }

    cout << "y = " << a << " + " << b << "x\n";
    fout << "y = " << a << " + " << b << "x\n";

    double Y = predict(a, b, X);

    cout << "Predicted Y = " << Y << " for X = " << X << "\n\n";
    fout << "Predicted Y = " << Y << " for X = " << X << "\n\n";

    testCase++;
}


    fin.close();
    fout.close();

    return 0;
}
```

<a id="linear-regression-input"></a>
#### Linear Regression Input

```
5
1 2
2 3
3 5
4 4
5 6
10

6
1 3
2 5
3 7
4 9
5 11
6 13
8

4
5 1
5 2
5 3
5 4
10
```

<a id="linear-regression-output"></a>
#### Linear Regression Output

```
Test Case 1
y = 1.3 + 0.9x
Predicted Y = 10.3 for X = 10

Test Case 2
y = 1 + 2x
Predicted Y = 17 for X = 8

Test Case 3
Error: Regression cannot be computed.
```

---

<a id="polynomial-regression"></a>
### Polynomial Regression

<a id="polynomial-regression-theory"></a>
#### Polynomial Regression Theory

```
❖ Polynomial regression is a numerical method used to fit a curve of degree m to a given set of experimental data points.

❖ The data is assumed to follow a polynomial relationship of the form:
\[
y = a_0 + a_1 x + a_2 x^2 + \cdots + a_m x^m
\]

❖ The unknown constants \( a_0, a_1, \dots, a_m \) are determined using the **Principle of Least Squares**, which states that the sum of the squares of the deviations between the observed and calculated values is minimum.

❖ The normal equations obtained are:
\[
\sum y = na_0 + a_1 \sum x + a_2 \sum x^2 + \cdots + a_m \sum x^m
\]
\[
\sum xy = a_0 \sum x + a_1 \sum x^2 + a_2 \sum x^3 + \cdots + a_m \sum x^{m+1}
\]
\[
\vdots
\]
\[
\sum x^m y = a_0 \sum x^m + a_1 \sum x^{m+1} + \cdots + a_m \sum x^{2m}
\]

❖ Solving these equations gives the values of the constants.

❖ Substituting these constants gives the best-fit polynomial curve.

---

Algorithm: 

Step 1: Arrange the given data points \((x, y)\) in tabular form.  
Step 2: Compute the required powers of \( x \) and products \( x^k y \).  
Step 3: Calculate the required summations.  
Step 4: Form the normal equations.  
Step 5: Solve the system of equations to determine the constants.  
Step 6: Write the equation of the best-fit polynomial.

```

<a id="polynomial-regression-code"></a>
#### Polynomial Regression Code

```cpp

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
```

<a id="polynomial-regression-input"></a>
#### Polynomial Regression Input

```
2

5 2
1 1
2 4
3 9
4 16
5 25
6

4 3
1 1
2 8
3 27
4 64
5
```

<a id="polynomial-regression-output"></a>
#### Polynomial Regression Output

```
Test Case 1:
Polynomial Degree: 2
Function: y = 0.0000 + -0.0000*x^1 + 1.0000*x^2
Predicted y for x = 6.0000: 36.0000

Test Case 2:
Polynomial Degree: 3
Function: y = 0.0000 + -0.0000*x^1 + 0.0000*x^2 + 1.0000*x^3
Predicted y for x = 5.0000: 125.0000
```

---

<a id="transcendental-regression"></a>
### Transcendental Regression

<a id="transcendental-regression-theory"></a>
#### Transcendental Regression Theory

```
Theory here
```

<a id="transcendental-regression-code"></a>
#### Transcendental Regression Code

```cpp
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
```

<a id="transcendental-regression-input"></a>
#### Transcendental Regression Input

```
3
1 5
1 2 3 4 5
2.718 7.389 20.086 54.598 148.413
6
2 4
1 2 3 4
0 0.6931 1.0986 1.3863
5
3 4
1 2 3 4
2 4 9 16
5
```

<a id="transcendental-regression-output"></a>
#### Transcendental Regression Output

```
Test Case 1:
Exponential Regression
Function: y = 0.999919 * e^(1.000021 * x)
Predicted y for x=6.000000: 403.446792

Test Case 2:
Logarithmic Regression
Function: y = -0.000017 + 1.000004 * ln(x)
Predicted y for x=5.000000: 1.609428

Test Case 3:
Power-law Regression
Function: y = 1.780428 * x^1.492058
Predicted y for x=5.000000: 19.652962
```

---

<a id="numerical-integration"></a>
## Numerical Integration

<a id="simpsons-13-rule"></a>
### Simpson's 1/3 Rule

<a id="simpson-13-theory"></a>
#### Simpson 1/3 Theory

```
Simpson’s 1/3 Rule is a numerical integration technique used to approximate the definite
integral of a function f(x) over a finite interval by fitting a parabola through three
consecutive data points.

❖ If the function f(x) is continuous over the interval a ≤ x ≤ b, then the definite
integral ∫ₐᵇ f(x) dx can be approximated using Simpson’s 1/3 Rule.

❖ The interval [a, b] is divided into an even number of subintervals of equal width h, where  
h = (b − a)/n.

❖ Let the ordinates be  
y₀ = f(x₀), y₁ = f(x₁), y₂ = f(x₂), … , yₙ = f(xₙ).

❖ The formula is:

∫ₐᵇ f(x) dx = (h/3) [ y₀ + yₙ + 4(y₁ + y₃ + … + yₙ₋₁) + 2(y₂ + y₄ + … + yₙ₋₂) ]

Algorithm:
Step 1: Choose the limits of integration a and b.  
Step 2: Divide the interval [a, b] into an even number of subintervals n.  
Step 3: Compute the step size h = (b − a)/n.  
Step 4: Calculate the function values y₀, y₁, … , yₙ.  
Step 5: Substitute the values into Simpson’s 1/3 formula.  
Step 6: Evaluate the approximate value of the integral.

```

<a id="simpson-13-code"></a>
#### Simpson 1/3 Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x) {
    return x*x + 2*x + 1;
}

double simpson13(double a, double b, int n) {
    if (n % 2 != 0) return -1;

    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        double x = a + i*h;
        if (i % 2 == 0)
            sum += 2 * f(x);
        else
            sum += 4 * f(x);
    }

    return (h / 3) * sum;
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

        double result = simpson13(a, b, n);

        cout << fixed << setprecision(3);
        fout << fixed << setprecision(3);

        fout << "Test Case " << t << ":\n";
        cout << "Test Case " << t << ":\n";

        if (result == -1)
        {
            fout << "Error: Number of intervals (n) must be even.\n\n";
            cout << "Error: Number of intervals (n) must be even.\n\n";

        }

        else
        {
            fout << "Approximate integral from " << a << " to " << b << " is = " << result << "\n\n";
            cout << "Approximate integral from " << a << " to " << b << " is = " << result << "\n\n";

        }


    }

    fin.close();
    fout.close();

    return 0;
}
```

<a id="simpson-13-input"></a>
#### Simpson 1/3 Input

```
2

0 2 4

1 3 6
```

<a id="simpson-13-output"></a>
#### Simpson 1/3 Output

```
Test Case 1:
Approximate integral from 0.000 to 2.000 is = 8.667

Test Case 2:
Approximate integral from 1.000 to 3.000 is = 18.667
```

---

<a id="simpsons-38-rule"></a>
### Simpson's 3/8 Rule

<a id="simpson-38-theory"></a>
#### Simpson 3/8 Theory

```
Simpson’s 3/8 Rule is a numerical integration method used to approximate definite integrals
by fitting a cubic polynomial through four consecutive data points.

❖ If the function f(x) is continuous over the interval a ≤ x ≤ b, then the definite
integral ∫ₐᵇ f(x) dx can be approximated using Simpson’s 3/8 Rule.

❖ The interval [a, b] is divided into several equal subintervals such that the total
number of subintervals is a multiple of 3.

❖ Let h = (b − a)/n be the common interval width.

❖ The formula is:

∫ₐᵇ f(x) dx = (3h/8) [ y₀ + yₙ + 3(y₁ + y₂ + y₄ + y₅ + …) + 2(y₃ + y₆ + …) ]

Algorithm:
Step 1: Choose the limits of integration a and b.  
Step 2: Divide the interval [a, b] into n subintervals where n is a multiple of 3.  
Step 3: Compute the step size h = (b − a)/n.  
Step 4: Calculate the function values y₀, y₁, … , yₙ.  
Step 5: Substitute the values into Simpson’s 3/8 formula.  
Step 6: Evaluate the approximate value of the integral.

```

<a id="simpson-38-code"></a>
#### Simpson 3/8 Code

```cpp
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
```

<a id="simpson-38-input"></a>
#### Simpson 3/8 Input

```
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
```

<a id="simpson-38-output"></a>
#### Simpson 3/8 Output

```
Test Case 1:
Error: Number of intervals (n) must be a multiple of 3.

Test Case 2:
Approximate integral from 1.000 to 3.000 is = 18.667
```

---
