#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>

#define MAX(A,B) ( (A) > (B) ? (A):(B) );

vector<vector<int> > Identity(int n){return Identity(n);}
namespace vba
{ 
    int Pricing_European_Options(double S0, double Maturity, double Strike, double r, double divid,double sigma, double* ptprice, 
                                double* ptrdelta, double* ptrgamma, double* ptrtheta, double* ptrrho, 
                                double* ptrvega, int nTime, int nSpace){return vba::Pricing_European_Options(S0,
                                Maturity,Strike, r, divid,sigma,ptprice,ptrdelta,ptrgamma,ptrtheta,ptrrho, 
                                ptrvega, nTime,nSpace);}
}