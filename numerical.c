#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ACC 1e-12
int main (void){
    double a, x, xt, res;
    char buf[64];
    double f(double, double), dfdx(double);
    printf("input a positive number : ");
    if (fgets(buf, sizeof(buf), stdin) == NULL ||
        sscanf(buf, "%lf", &a) != 1){
        fputs("Input error\n", stderr);
        exit(-1);
    }
    x = a;
    do {
        xt = x;
        x = xt - f(x, a) / dfdx(x);
        res = fabs(x-xt);
        printf("x=%25.16f, |xt-x|=%25.16f\n", x, res);
    }while (res >= ACC);
    printf("sqrt(%f)=%25.16f\n", a, x);
    return 0;
}

    double f(double x, double a){
        return x*x - a;
    }
    double dfdx(double x){
        return 2*x;
    }
