//
//  main.nonthermal.c
//  
//
//  Created by Zhaowei Zhang on 6/5/15.
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>

/* other C header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

//global variables
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


//double n_peak(double nu);
double K_s(double gamma, double n, double nu);
double my_Bessel_J(double n, double x);
double my_Bessel_dJ(double n, double x);
//double MJ_f(double gamma);
double I(double gamma, double n, double nu);
double trapez_gamma(double min, double max, double n, double nu);
double trapez_n(double min, double max, double nu);
double trapez_norm(double min, double max);

double gamma_integrand(double gamma, double n, double nu);
double gamma_integration_result(double n, double nu);
double n_summation(double nu);
double n_integration(double n_minus, double nu);
double Power_law_no_norm(double gamma);
double Power_law_f(double gamma);
double n_peak_zhaowei(double nu);

int main(int argc, char *argv[])
{
    //define parameters of calculation
    double nu_c = (e * B)/(2. * M_PI * m * c);
    double i;
    for (i = 1.; i < 1e6; i = i * 5.)
    {
        n_summation(i * nu_c);
    }
//    double nu = 1. * nu_c;

//    printf("\n%e\n", trapez_norm(gamma_min, gamma_max));
//    double check = gamma_integration_result(100, nu);
//    printf("%f", check);
    return 0;
}


double K_s(double gamma, double n, double nu)
{
    double nu_c = (e * B)/(2. * M_PI * m * c);
    double beta = sqrt(1. - 1./(gamma*gamma));
    double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
    double M = (cos(theta) - beta * cos_xi)/sin(theta);
    double N = beta * sqrt(1 - (cos_xi*cos_xi));
    double z = (nu * gamma * beta * sin(theta) * sqrt(1. - cos_xi*cos_xi))/nu_c;
    double K_xx = M*M * pow(my_Bessel_J(n, z), 2.);
    double K_yy = N*N * pow(my_Bessel_dJ(n, z), 2.);
    double ans = K_xx + K_yy;
    return ans;
}

double Power_law_no_norm(double gamma)
{
    double prefactor = n_e_nt * (p - 1.) / (pow(gamma_min, 1. - p) - pow(gamma_max, 1. - p));
    double body = pow(gamma, -p) * exp(- gamma / gamma_cutoff);
    double ans = prefactor * body;
    return ans;
}

double trapez_norm(double min, double max)
{
    int i;
    float inteval, sum = 0., x;
    int divisions = 1000;
    
    inteval = (max-min) / (divisions-1);
    
    for (i=2; i<divisions; i++)
    {
        x =  min + inteval * (i - 1);
        sum = sum + Power_law_no_norm(x);
//        printf("\n%e\n", sum);
    }
    return(sum);
}



double Power_law_f(double gamma)
{
    double beta = sqrt(1. - 1./(gamma*gamma));
    double prefactor_47 = 1./(pow(m,3) * pow(c,3) * gamma * gamma * beta);
//    double norm_2 = 1. / trapez_norm(1, 1e6);
    double norm_1 = 1.00199514442;
//    printf("%f", norm_1);
//    printf("%f", norm_2);
    
    double power_f = (1./(4. * M_PI)) * norm_1 * prefactor_47 * Power_law_no_norm(gamma);
    
    return power_f;
}


double I(double gamma, double n, double nu)
{
    double nu_c = (e * B)/(2. * M_PI * m * c);
    double beta = sqrt(1. - 1./(gamma*gamma));
    double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
    double ans = (2. * M_PI * e*e * nu*nu)/c * (pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta * 2. * M_PI) * Power_law_f(gamma) * K_s(gamma, n, nu);
    return ans;
}

double gamma_integrand(double gamma, double n, double nu)
{
    double nu_c = (e * B)/(2. * M_PI * m * c);
    double beta = sqrt(1. - 1./(gamma*gamma));
    double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
    double prefactor = 1./(nu * beta * fabs(cos(theta)));
    double ans = prefactor * I(gamma, n, nu);
//    printf("\n%e\n", ans);
    return ans;
}

double gamma_integration_result(double n, double nu)
{
    double nu_c = (e * B)/(2. * M_PI * m * c);
    double gamma_minus = ((n*nu_c)/nu - fabs(cos(theta))*sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(theta), 2.)))/(pow(sin(theta), 2));
    double gamma_plus  = ((n*nu_c)/nu + fabs(cos(theta))*sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(theta), 2.)))/(pow(sin(theta), 2));
    double result = trapez_gamma(gamma_minus, gamma_plus, n, nu);
//    printf("\n%e\n", result);
    return result;
}

double trapez_gamma (double min, double max, double n, double nu)
{
    int i;
    float interval, sum=0., x;
    int divisions = 1000;
    
    interval = ((max-min) / (divisions-1));
    
    for (i=2; i<divisions; i++)
   	{
      		x    = min + interval * (i-1);
      		sum += gamma_integrand(x, n, nu)*interval;
//            printf("\n%e\n", sum);
   	}
    
    //sum += 0.5 *(gamma_integrand(min, n, nu) + gamma_integrand(max, n, nu)) * interval;
   	return (sum);
}

double trapez_n(double min, double max, double nu)
{
    int i;
    float interval, sum=0., x;
    int divisions = 1000;
    
    interval = ((max-min) / (divisions-1));
    
    for (i=2; i<divisions; i++)
   	{
      		x    = min + interval * (i-1);
      		sum += gamma_integration_result(x, nu)*interval;
//            printf("\n%e\n", gamma_integration_result(x, nu));
   	}
    
    //sum += 0.5 *(gamma_integration_result(min, nu) + gamma_integration_result(max, nu)) * interval;
   	return (sum);
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
    
    
}



double n_integration(double n_minus, double nu)
{
    if(n_max < n_minus)
    {
        n_max = n_minus;
    }
    
    int n_peak = 20;
    
    double ans = trapez_n(n_max, C * n_peak_zhaowei(nu), nu);
//    printf("\n%e\n", ans);
    return ans;
}

double n_summation(double nu)
{
    double j_nu = 0.;
    double nu_c = (e * B)/(2. * M_PI * m * c);
    double n_minus = (nu/nu_c) * fabs(sin(theta));
    printf("%f", n_minus);
    int x;
    for(x = (int)(n_minus+1.); x <= n_max + (int)n_minus ; x++)
    {
        j_nu = j_nu + gamma_integration_result(x, nu);
//        printf("\n%e\n", j_nu);
//        printf("%i\n",x);
    }
    
//    printf("\n%e\n", j_nu);
    j_nu = j_nu + n_integration(n_minus, nu);
    printf("\n%e\n", j_nu);
    return j_nu;
}
