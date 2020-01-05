#include "pch.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <complex>
#include "Header.hpp"

using namespace std;

#define Call = 1;
#define Put = 2;

double Phi(int option, double x, double Strike)
{
	double phi;
	switch (option)
	{
	case 1:
		phi = MAX((exp(x) - Strike), 0);
		break;

	case 2:
		phi = MAX((-exp(x) + Strike), 0);
		break;
	}
	return phi;
}


std::vector<double>solve_tridiag(const  std::vector<double>& low,const  std::vector<double>& diag, 
	                             const  std::vector<double>& up, std::vector<double>& rightside_result)
{
	int N = diag.size();
	if ((low.size() != N) || (low.size() != N) || (rightside_result.size() != N))
	{
		std::cout << "les tailles des vecteurs sont incorrectes!" << std::endl;
		std::exit(1);
	}

	std::vector<double> A(N - 1), B(N - 1);

	// descente
	const std::vector<double> & rightside = rightside_result;
	A[0] = up[0] / diag[0];
	B[0] = rightside[0] / diag[0];
	double denom;
	for (int i = 1;i < N - 1;i++)
	{
		denom = diag[i] - low[i] * A[i - 1];
		A[i] = up[i] / denom;
		B[i] = (rightside[i] - low[i] * B[i - 1]) / denom;
	};

	// remontee
	std::vector<double> & result = rightside_result;
	result[N - 1] = (rightside[N - 1] - low[N - 1] * B[N - 2]) / (diag[N - 1] - low[N - 1] * A[N - 2]);
	for (int i = N - 2;i >= 0;i--)
	{
		result[i] = -A[i] * result[i + 1] + B[i];
	};

	return(result);
}


double finum::Price_option(int optiontype,double s0, double Maturity, double Strike, double r, double  divid, double sigma, int nTime, int nSpace)
{


	double h = Maturity / nTime,
		b = r - divid - SQR(sigma) / 2.0,
		l = abs(b) * Maturity + 4.8 * sigma * sqrt(Maturity), //3.5 proba for the gaussian process S to leave the definition domain +/-l
		k = 2.0 * l / nSpace ,
		alpha = (SQR(sigma) / (2.0*SQR(k))) - (b / (2 * k));
	double beta = -((SQR(sigma) / SQR(k)) + r),
		gamma = (SQR(sigma) / (2 * SQR(k))) + (b / (2 * k));
	double pM_up = -0.5 * h * gamma,   
		pM_down = -0.5 * h * alpha,
		pM_diag = 1 - 0.5 * h * beta,
		pN_down = 0.5*h*alpha,
		pN_diag = 1 + 0.5*h*beta,
		pN_up = 0.5*h*gamma;
	std::vector<double> x(nSpace), Payoff(nSpace), Price(nSpace), Price1(nSpace);
	double T = Maturity, final_price;

	//Price at time T
	for (int i = 1; i <= nSpace; i++)
	{
		x[i - 1] = log(s0) - l + i * k;
		Price[i - 1] = Phi(optiontype, x[i - 1], Strike);
	}


	// Dirichlet boundary condition, u(t_i,x+l)= K*exp(r*(T-t)), u(t_i,x-l)= 0
	Price[0] = Price[0] + 0;
	//at the beginning, t=T so  exp(r*(T-t))=1
	Price[nSpace - 1] = Price[nSpace - 1];  //a verifier


	// ces vecteurs representent la matrice tridiagonale du systeme:
	std::vector<double> low(nSpace);      // diagonale inferieure
	std::vector<double> diag(nSpace); // diagonale principale
	std::vector<double> up(nSpace);       // diagonale superieure

	 // conditions au bord: propagation des valeurs initiales:
	low[0] = 0;
	diag[0] = pM_diag;
	up[nSpace - 1] = 0;
	diag[nSpace - 1] = pM_diag;
	
	
	for (int i = 1; i < nSpace; i++)
	{
		up[i - 1] = pM_up;
		diag[i] = pM_diag;
		low[i] = pM_down;
	};
	

	for (int t = nTime; t > 0; t--)
	{
		Price1[0] = Price[0];
		Price1[nSpace - 1] = pN_down * Price[nSpace - 2] + pN_diag * Price[nSpace - 1];

		for (int i = 1; i < nSpace - 1; i++)
		{
			Price1[i] = pN_down * Price[i - 1] + pN_diag * Price[i] + pN_up * Price[i + 1];
		}
		
		Price = solve_tridiag(low, diag, up, Price1);
		// Dirichlet boundary condition, u(t_i,x+l)= K*exp(r*(T-t)), u(t_i,x-l)= 0
		Price[nSpace - 1] = Price[nSpace - 1] - (Strike*exp(r*(T - (t - 1)*h))*pM_up);
	}

	if (nSpace % 2 == 0)
	{
		final_price = Price[(nSpace ) / 2 - 1];
	}
	else
	{
		final_price = (Price[(nSpace - 1) / 2] + Price[(nSpace - 3) / 2]) / 2;
	}


	/*for (int i = 1; i < nSpace; i++)
	{

		cout <<  Price[i] ;

		cout << endl;
	}*/

	/*for (int i = 1; i < nSpace; i++)
	{

		ptdelta[i] = exp(-x[i])*Price

	};


	*ptdelta = Price[nSpace / 2];
	*ptgamma = Price[nSpace / 2];
	*/
	return (final_price);

}

double finum:: delta(int optiontype, double S, double Maturity, double Strike, double r,
	                 double  divid, double sigma, int nTime, int nSpace)

{
	double delta_S, price1, price2, S1, delta;

	delta_S = S * (Maturity / nSpace);
	S1 =  S + delta_S;
	price1 = finum::Price_option(optiontype, S, Maturity, Strike, r, divid, sigma, nTime, nSpace);
	price2 = finum::Price_option(optiontype, S1, Maturity, Strike, r, divid, sigma, nTime, nSpace);
	delta = (price2 - price1) / delta_S;

	return(delta);
}
int main()
{

	double price, s0, Maturity, Strike, eqr, eqdiv, sigma_ref, volimpl, impliedDelta, impliedGamma, rebate;
	s0 = 105.;
	Maturity = 1.;
	Strike = 100.;
	eqr = 0.2;
	eqdiv = 0.;
	sigma_ref = .2;


	int nTime = 500;
	int nSpace = 499;

	const int nspace_2 = 10;
	double delta;
	//delta = 0;
	price=finum::Price_option(1,s0, Maturity, Strike, eqr, eqdiv, sigma_ref, nTime, nSpace);
	cout << "The price of ImplPriceNone is " << price << " " << endl;
	delta= finum::delta(1, s0, Maturity, Strike, eqr, eqdiv, sigma_ref, nTime, nSpace);
	
	cout << "delta is " << delta;



	return 0;
}
