#include<iostream>
#include<Eigen/core>
#include"mars/mars.h"
int main(){
	mars_algorithm::MARS mars;
	mars.mars_input(6.0, 3.0, 512, 256, 2,
		44, 3, 0.1,
		0.00001, 8,
		1,40);
	mars.solve();
	
} 
//# include <cmath>
//# include <complex>
//# include <cstdlib>
//# include <ctime>
//# include <iomanip>
//# include <iostream>
//
//using namespace std;
//
//#include"mars/mars.h"
//#include"mars/fsolve.h"
//int main();
//void fsolve_test1();
//void f1(int n, double x[], double fvec[]);
//void fsolve_test2();
//void f2(int n, double x[], double fvec[]);
//void fsolve_test3();
//void f3(int n, double x[], double fvec[]);
//void fsolve_test4();
//void f4(int n, double x[], double fvec[]);
//void r8vec2_print(int n, double a1[], double a2[], string title);
//void timestamp();
//
////****************************************************************************80
//
//int main()
//
////****************************************************************************80
////
////  Purpose:
////
////    fsolve_test() tests fsolve().
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    08 April 2021
////
////  Author:
////
////    John Burkardt
////
//{
//    timestamp();
//    cout << "\n";
//    cout << "fsolve_test()\n";
//    cout << "  C++ version:\n";
//    cout << "  Test fsolve().\n";
//
//    fsolve_test1();
//    fsolve_test2();
//    fsolve_test3();
//    fsolve_test4();
//    //
//    //  Terminate.
//    //
//    cout << "\n";
//    cout << "fsolve_test():\n";
//    cout << "  Normal end of execution.\n";
//    cout << "\n";
//    timestamp();
//
//    return 0;
//}
////****************************************************************************80
//
//void fsolve_test1()
//
////****************************************************************************80
////
////  Purpose:
////
////    fsolve_test1() tests fsolve() with a system of 1 equation.
////
////  Discussion:
////
////    This is a version of Kepler's equation.
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    08 April 2021
////
////  Author:
////
////    John Burkardt
////
//{
//    double* fx;
//    int info;
//    int lwa;
//    int n = 1;
//    double tol = 0.00001;
//    double* wa;
//    double* x;
//
//    lwa = (n * (3 * n + 13)) / 2;
//    fx = new double[n];
//    wa = new double[lwa];
//    x = new double[n];
//
//    cout << "\n";
//    cout << "fsolve_test1():\n";
//    cout << "  fsolve() solves a nonlinear system of 1 equation.\n";
//
//    x[0] = 0.0;
//    f1(n, x, fx);
//
//    r8vec2_print(n, x, fx, "  Initial X, F(X)");
//
//    info = fsolve(f1, n, x, fx, tol, wa, lwa);
//
//    cout << "\n";
//    cout << "  Returned value of INFO = " << info << "\n";
//    r8vec2_print(n, x, fx, "  Final X, FX");
//    //
//    //  Free memory.
//    //
//    delete[] fx;
//    delete[] wa;
//    delete[] x;
//
//    return;
//}
////****************************************************************************80
//
//void f1(int n, double x[], double fx[])
//
////****************************************************************************80
////
////  Purpose:
////
////    f1() evaluates a system of 1 nonlinear equation.
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    08 April 2010
////
////  Author:
////
////    John Burkardt
////
////  Input:
////
////    int N, the number of variables.
////
////    double X[N], the variable values.
////
////  Output:
////
////    double FX[N], the function values.
////
//{
//    double e;
//    double m;
//
//    e = 0.8;
//    m = 5.0;
//
//    fx[0] = x[0] - m - e * sin(x[0]);
//
//    return;
//}
////****************************************************************************80
//
//void fsolve_test2()
//
////****************************************************************************80
////
////  Purpose:
////
////    fsolve_test2() tests fsolve() with a system of 2 equations.
////
////  Discussion:
////
////    This is an example of what your main program would look
////    like if you wanted to use MINPACK to solve N nonlinear equations
////    in N unknowns.  In this version, we avoid computing the jacobian
////    matrix, and request that MINPACK approximate it for us.
////
////    The set of nonlinear equations is:
////
////      x1 * x1 - 10 * x1 + x2 * x2 + 8 = 0
////      x1 * x2 * x2 + x1 - 10 * x2 + 8 = 0
////
////    with solution x1 = x2 = 1
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    08 April 2021
////
////  Author:
////
////    John Burkardt
////
//{
//    double* fx;
//    int info;
//    int lwa;
//    int n = 2;
//    double tol = 0.00001;
//    double* wa;
//    double* x;
//
//    lwa = (n * (3 * n + 13)) / 2;
//    fx = new double[n];
//    wa = new double[lwa];
//    x = new double[n];
//
//    cout << "\n";
//    cout << "fsolve_test2():\n";
//    cout << "  fsolve() solves a nonlinear system of 2 equations.\n";
//
//    x[0] = 3.0;
//    x[1] = 0.0;
//    f2(n, x, fx);
//
//    r8vec2_print(n, x, fx, "  Initial X, F(X)");
//
//    info = fsolve(f2, n, x, fx, tol, wa, lwa);
//
//    cout << "\n";
//    cout << "  Returned value of INFO = " << info << "\n";
//    r8vec2_print(n, x, fx, "  Final X, FX");
//    //
//    //  Free memory.
//    //
//    delete[] fx;
//    delete[] wa;
//    delete[] x;
//
//    return;
//}
////****************************************************************************80
//
//void f2(int n, double x[], double fx[])
//
////****************************************************************************80
////
////  Purpose:
////
////    f2() evaluates a system of 2 nonlinear equations.
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    08 April 2010
////
////  Author:
////
////    John Burkardt
////
////  Input:
////
////    int N, the number of variables.
////
////    double X[N], the variable values.
////
////  Output:
////
////    double FX[N], the function values.
////
//{
//    fx[0] = x[0] * x[0] - 10.0 * x[0] + x[1] * x[1] + 8.0;
//    fx[1] = x[0] * x[1] * x[1] + x[0] - 10.0 * x[1] + 8.0;
//
//    return;
//}
////****************************************************************************80
//
//void fsolve_test3()
//
////****************************************************************************80
////
////  Purpose:
////
////    fsolve_test3() tests fsolve() with a system of 4 equations.
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    08 April 2021
////
////  Author:
////
////    John Burkardt
////
//{
//    double* fx;
//    int i;
//    int info;
//    int lwa;
//    int n = 4;
//    double tol = 0.00001;
//    double* wa;
//    double* x;
//
//    lwa = (n * (3 * n + 13)) / 2;
//    fx = new double[n];
//    wa = new double[lwa];
//    x = new double[n];
//
//    cout << "\n";
//    cout << "fsolve_test3():\n";
//    cout << "  fsolve() solves a nonlinear system of 4 equations.\n";
//
//    for (i = 0; i < n; i++)
//    {
//        x[i] = 0.0;
//    }
//    f3(n, x, fx);
//
//    r8vec2_print(n, x, fx, "  Initial X, F(X)");
//
//    info = fsolve(f3, n, x, fx, tol, wa, lwa);
//
//    cout << "\n";
//    cout << "  Returned value of INFO = " << info << "\n";
//    r8vec2_print(n, x, fx, "  Final X, FX");
//    //
//    //  Free memory.
//    //
//    delete[] fx;
//    delete[] wa;
//    delete[] x;
//
//    return;
//}
////****************************************************************************80
//
//void f3(int n, double x[], double fx[])
//
////****************************************************************************80
////
////  Purpose:
////
////    f3() evaluates a system of 4 nonlinear equations.
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    08 April 2010
////
////  Author:
////
////    John Burkardt
////
////  Input:
////
////    int N, the number of variables.
////
////    double X[N], the variable values.
////
////  Output:
////
////    double FX[N], the function values.
////
//{
//    int i;
//
//    for (i = 0; i < n; i++)
//    {
//        fx[i] = pow(x[i] - (i + 1), 2);
//    }
//
//    return;
//}
////****************************************************************************80
//
//void fsolve_test4()
//
////****************************************************************************80
////
////  Purpose:
////
////    fsolve_test4() tests fsolve() with a system of 8 equations.
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    08 April 2021
////
////  Author:
////
////    John Burkardt
////
//{
//    double* fx;
//    int i;
//    int info;
//    int lwa;
//    int n = 500;
//    double tol = 0.00001;
//    double* wa;
//    double* x;
//
//    lwa = (n * (3 * n + 13)) / 2;
//    fx = new double[n];
//    wa = new double[lwa];
//    x = new double[n];
//
//    cout << "\n";
//    cout << "fsolve_test4():\n";
//    cout << "  fsolve() solves a nonlinear system of 8 equations.\n";
//
//    for (i = 0; i < n; i++)
//    {
//        x[i] = 0.0;
//    }
//    f4(n, x, fx);
//
//    r8vec2_print(n, x, fx, "  Initial X, F(X)");
//
//    info = fsolve(f4, n, x, fx, tol, wa, lwa);
//
//    cout << "\n";
//    cout << "  Returned value of INFO = " << info << "\n";
//    r8vec2_print(n, x, fx, "  Final X, FX");
//    //
//    //  Free memory.
//    //
//    delete[] fx;
//    delete[] wa;
//    delete[] x;
//
//    return;
//}
////****************************************************************************80
//
//void f4(int n, double x[], double fx[])
//
////****************************************************************************80
////
////  Purpose:
////
////    f4() evaluates a system of 8 nonlinear equations.
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    08 April 2010
////
////  Author:
////
////    John Burkardt
////
////  Input:
////
////    int N, the number of variables.
////
////    double X[N], the variable values.
////
////  Output:
////
////    double FX[N], the function values.
////
//{
//    int i;
//
//    for (i = 0; i < n; i++)
//    {
//        fx[i] = (3.0 - 2.0 * x[i]) * x[i] + 1.0;
//
//        if (0 < i)
//        {
//            fx[i] = fx[i] - x[i - 1];
//        }
//
//        if (i < n - 1)
//        {
//            fx[i] = fx[i] - 2.0 * x[i + 1];
//        }
//    }
//
//    return;
//}
////****************************************************************************80
//
//void r8vec2_print(int n, double a1[], double a2[], string title)
//
////****************************************************************************80
////
////  Purpose:
////
////    r8vec2_print() prints an R8VEC2.
////
////  Discussion:
////
////    An R8VEC2 is a dataset consisting of N pairs of real values, stored
////    as two separate vectors A1 and A2.
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    19 November 2002
////
////  Author:
////
////    John Burkardt
////
////  Input:
////
////    int N, the number of components of the vector.
////
////    double A1[N], double A2[N], the vectors to be printed.
////
////    string TITLE, a title.
////
//{
//    int i;
//
//    cout << "\n";
//    cout << title << "\n";
//    cout << "\n";
//    for (i = 0; i <= n - 1; i++)
//    {
//        cout << setw(6) << i
//            << ": " << setw(14) << a1[i]
//            << "  " << setw(14) << a2[i] << "\n";
//    }
//
//    return;
//}
////****************************************************************************80
//
//void timestamp()
//
////****************************************************************************80
////
////  Purpose:
////
////    timestamp() prints the current YMDHMS date as a time stamp.
////
////  Example:
////
////    31 May 2001 09:45:54 AM
////
////  Licensing:
////
////    This code is distributed under the GNU LGPL license.
////
////  Modified:
////
////    08 July 2009
////
////  Author:
////
////    John Burkardt
////
//{
//# define TIME_SIZE 40
//
//    static char time_buffer[TIME_SIZE];
//    const struct std::tm* tm_ptr;
//    std::time_t now;
//
//    now = std::time(NULL);
//    tm_ptr = std::localtime(&now);
//
//    std::strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr);
//
//    std::cout << time_buffer << "\n";
//
//    return;
//# undef TIME_SIZE
//}
