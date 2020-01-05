#include <cmath>
#include <ctime>
#include <vector>
#include <iterator>

#ifndef _MAFIN_H
#define _MAFIN_H

typedef unsigned char boolean;

#define MAX(A,B) ( (A) > (B) ? (A):(B) )
#define MIN(A,B) ( (A) < (B) ? (A):(B) )
#define SQR(X) ((X)*(X))

/**
	@brief Thomas algorithm for tridiagonal systems Ax=b

	@param low	lower diagonal of A
	@param diag	 diagonal of A
	@param up	upper diagonal of A
 */
std::vector<double> solve_tridiag(const std::vector<double>& low, const std::vector<double>& diag, const
	std::vector<double>& up, std::vector<double>& rightside_result);


namespace finum
{

	double  Price_option(int type, double s0, double Maturity, double Strike, double r, double  divid, double sigma, int nTime, int nSpace);
	double delta(int type, double s0, double Maturity, double Strike, double r, double  divid, double sigma, int nTime, int nSpace, double delta_S);
	double gamma(int optiontype, double S, double Maturity, double Strike, double r, double  divid, double sigma, int nTime, int nSpace, double delta_S);
	double theta(int optiontype, double S, double Maturity, double Strike, double r, double  divid, double sigma, int nTime, int nSpace, double delta_t);
	double vega(int optiontype, double S, double Maturity, double Strike, double r, double  divid, double sigma, int nTime, int nSpace, double delta_v);
	double rho(int optiontype, double S, double Maturity, double Strike, double r, double  divid, double sigma, int nTime, int nSpace, double delta_r);
}
/*#define int_dll  __declspec(dllexport) int __stdcall
#define void_dll  __declspec(dllexport) void  __stdcall
#define double_dll __declspec(dllexport) double __stdcall
extern "C" {


	int_dll price_option(int type,double s0, double Maturity,
		double Strike, double r, double  divid, double sigma,
		double *ptprice, double *ptdelta, double *ptgamma, int nTime, int nSpace) {
		return finum::Price_option(type,s0, Maturity, Strike, r,
			divid, sigma, ptprice, ptdelta,
			ptgamma, nTime, nSpace);
	}


	//int_dll main(){return finum::main() ;} 
}*/

#endif