#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <complex>
#include <vba_project.hpp>
 
 using namespace std;

  vector<vector<int> > Identity(int n)
    { 
       vector<vector<int> > matrix;
       matrix = vector<vector<int> >(n, vector<int>(n,0));

        // Set the diagonal to be 1
        for(unsigned int t = 0; t < n; t++)
            matrix[t][t] = 1;
        return matrix;
    } 

  

 int vba::Pricing_European_Options(double S0,double Maturity, double Strike, double r, double divid,double sigma, double* ptprice, double* ptrdelta, double* ptrgamma, double* ptrtheta, double* ptrrho, double* ptrvega, int nTime, int nSpace)
 {
   double h = Maturity/nTime,
          b = r - divid - sigma*sigma / nTime;
   double l = abs(b)*Maturity + 3.89*sigma*sqrt(Maturity),
          k = 2*l/(nSpace),
          alpha = (sigma*sigma)/(2.0*k*k) - b/(2*k);
   double beta = (sigma*sigma)/(k*k)+r;
   double gamma = (sigma*sigma)/(2*k*k) + b/(2*k);

        std::vector<double> x(nSpace+1), Payoff(nSpace+1), Price(nSpace+1);
        for (int i = 0; i < nSpace; i++)
        {
            x[i] = log(S0) - l + i*k;
            Payoff[i] = MAX(exp(x[i])-Strike, 0);
            Price[i] = Payoff[i];
        }
        
  }
 