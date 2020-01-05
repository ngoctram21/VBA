#include <cmath>
#include <ctime>
#include <vector>
#include <iterator>

#define MAX(A,B) ( (A) > (B) ? (A):(B) );


vector<double> Phi (int option, vector<double> x, double Strike){return Phi(option, x,Strike);};
vector<vector<double>> M_diag(int n, double x){return M_diag(n,x);};
vector<vector<double>> M_low(int n, double x){return M_low(n,x);};
vector<vector<double>> M_up(int n, double x){return M_up(n,x);};
vector<vector<double>> M_sum(vector<vector<double>> A , vector<vector<double>> B){return M_sum(A, B);};
vector<vector<double>> M_theta_hA(int n, double a, double b, double c,int h){return M_theta_hA(n,a,b,c,h);};
vector<double> M_mult_v(vector<vector<double>> A, vector<double> u){return M_mult_v(A,u);};
vector<vector<double>> M_mult_M(vector<vector<double>> A, vector<vector<double>> B){return M_mult_M(A,B);};
vector<double> Decom_LU(vector<vector<double>> A, vector<double> b){return Decom_LU(A,b);};

namespace vba
{ 
    double Price(int option_type, double S0, double Maturity, double Strike, double r, double sigma, double divid, int nTime, int nSpace){return vba::Price(option_type, S0,Maturity,Strike,r,sigma,divid,nTime,nSpace);};
}