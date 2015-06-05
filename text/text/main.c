//
//  main.c
//  text
//
//  Created by Zhaowei Zhang on 6/5/15.
//  Copyright (c) 2015 Zhaowei Zhang. All rights reserved.
//

/* other C header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

double m = 9.1093826e-28;
double c = 2.99792458e10;
double theta_e = 10.;
double e = 4.80320680e-10;
double B = 30.;
double n_e = 1.;
double n_e_nt = 1.;
double theta = (M_PI  / 3.);
int C = 1;
double n_max = 30.;
double p = 3.;
double gamma_min = 1.;
double gamma_max = 1000.;
double gamma_cutoff = 1000.;
double check(double nu);

int main(int argc, const char * argv[])
{
    double nu = 1.;
    double hi = check(nu);
    printf("%f", hi);
    return 0;
}


double n_peak_zhaowei(double nu)
{
    double nu_c = (e * B)/(2. * M_PI * m * c);
    double n_overpeak = 91. * pow(21., log10(nu/nu_c));
    double n_minus = (nu/nu_c) * fabs(sin(theta));
    
    
    int    x_value[(int)n_overpeak];
    double y_value[(int)n_overpeak];
    
    double x_value_max = 0;
    double y_value_max = 0;
    
    double x_variable = 0;
    
    int step_div;
    
    int i, j, k;
    int a=0;
    int x_threshold;
    
    if(n_overpeak < 1000)
    {
        step_div = (int)n_overpeak;
    }
    else
    {
        step_div = 1000;
    }
    
    int step = (int)(n_overpeak / step_div);
    
    
    for(i = (int)n_minus; i <= (int)n_overpeak; i = i + step)
    {
        x_value[a] = i;
        y_value[a] = gamma_integration_result(i, nu);
        a = a  + 1;
    }
    
    for(j = 0; j < step_div; j++)
    {
        if(y_value[j] > y_value_max)
        {
            x_value_max = x_value[j];
            y_value_max = y_value[j];
            x_variable = j;
        }
    }
    
    double y_threshold = (1./1000.) * y_value_max;
    
    for(k = x_variable; k < step_div; k++)
    {
        if(y_value[k] < y_threshold)
        {
            x_threshold = x_value[k];
            return x_threshold;
        }
    }
    return 0;
    
}

