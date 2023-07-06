#pragma once
#include <iostream>
#include <cmath>

using namespace std;

typedef double (*Function)(double);         //для передачи различных функций в виде указателя
typedef double (*FunctionArr[])(double);    //для передачи массива указателей на функции в качестве параметров

//Functions function(double, double, double) = {f1, f2, f3, f4, f5, f6, f7, f8};

//
// ******************************** Формулы и производные ******************************
// 

//уравнение для Ньютона на одно уравнение
double equation1 (double x) {
    return x * x + x - 6;
}

double GetDifff(double (*equation)(double), double x, double dp = 1.0e-6)
{
    return 2 * x + 1;
}


//
//      p01---->p1---->p2<----p02
//              |      |
//              v      v
//      po3---->p3---->p4----->p4
// 
// 
//функции для исходных уравнений
//производные в системе из 4-ех узлов = 1, -1, -2 * k * Q
double f1(double p1, double Q01, double p01, double k)          //k можно будет настроить для каждого участка отдельно
{
    return p01 - p1 - Q01 * Q01 * k;
}

double f2(double p1, double p2, double Q12, double k)
{
    return p1 - p2 - Q12 * Q12 * k;
}

double f3(double p2, double p4, double Q24, double k)
{
    return p2 - p4 - Q24 * Q24 * k;
}

double f4(double p1, double p3, double Q13, double k)
{
    return p1 - p3 - Q13 * Q13 * k;
}

double f5(double p3, double p4, double Q34, double k)
{
    return p3 - p4 - Q34 * Q34 * k;
}

double f6(double p4, double p04, double Q04, double k)
{
    return p4 - p04 - Q04 * Q04 * k;
}

double f7(double p2, double Q02, double p02, double k)
{
    return p02 - p2 - Q02 * Q02 * k;
}

double f8(double p3, double Q03, double p03, double k)
{
    return p03 - p3 - Q03 * Q03 * k;
}

double f9(double p0, double Q03, double p03, double k) { return p0; };
double f10(double p0, double Q03, double p03, double k) { return p0; };
double f11(double p0, double Q03, double p03, double k) { return p0; };
double f12(double p0, double Q03, double p03, double k) { return p0; };

//производные для системы уравнений
double df1dp1(double Q01, double k)
{
    return -1;
}

double df1dQ01(double p1, double Q01, double k)
{
    return -2 * Q01 * k;
}

double df2dp1()
{
    return 1;
}

double df2dp2()
{
    return -1;
}

double df2dQ12(double Q12, double k)
{
    return -2 * Q12 * k;
}

double df3dp2()
{
    return 1;
}

double df3dp4()
{
    return -1;
}

double df3dQ24(double Q24, double k)
{
    return -2 * Q24 * k;
}

double df4dp1()
{
    return 1;
}

double df4dp3()
{
    return -1;
}

double df4dQ13(double Q13, double k)
{
    return -2 * Q13 * k;
}

double df5dp3()
{
    return 1;
}

double df5dp4()
{
    return -1;
}

double df5dQ34(double Q34, double k)
{
    return -2 * Q34 * k;
}

double df6dp4()
{
    return 1;
}

double df6dQ04(double Q04, double k)
{
    return -2 * Q04 * k;
}

double df7dp2()
{
    return -1;
}

double df7dQ02(double Q02, double k)
{
    return -2 * Q02 * k;
}

double df8dp3()
{
    return -1;
}

double df8dQ03(double Q03, double k)
{
    return -2 * Q03 * k;
}