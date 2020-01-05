#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <complex>
#include <vector>
 
 using namespace std;

 #define MAX(A,B) ( (A) > (B) ? (A):(B) );
 #define Call = 1;
 #define Put = 2;

double Phi (int option, double x, double Strike)
{
    double phi;
    switch (option)
    {
    case 1:
        phi = MAX((exp(x)-Strike), 0);
        break;
    
    case 2:
        phi = MAX((-exp(x)+Strike), 0);
        break;
    }
    return phi;
}

vector<vector<double>> M_diag(vector<double> x)
{
    int n = x.size();
    vector<vector<double> > matrix_diag;
    matrix_diag = vector<vector<double>>(n, vector<double>(n,0));
    for (int i = 0; i < n; i++)
    {
        matrix_diag[i][i] = x[i];
    }
    return matrix_diag;
}

vector<vector<double>> M_low(vector<double> x)
{
    int n = x.size()+1;
    vector<vector<double> > matrix_low;
    matrix_low = vector<vector<double>>(n, vector<double>(n,0));
    for (int i = 0; i < n-1; i++)
    {
        matrix_low[i+1][i] = x[i];
    }
    return matrix_low;
}

vector<vector<double>> M_up(vector<double> x)
{
    int n = x.size()+1;
    vector<vector<double> > matrix_up;
    matrix_up = vector<vector<double>>(n, vector<double>(n,0));
    for (int i = 0; i < n-1; i++)
    {
        matrix_up[i][i+1] = x[i];
    }
    return matrix_up;
}
 
vector<vector<double>> M_sum(vector<vector<double>> A , vector<vector<double>> B)
{
    int n = A.size();
    vector<vector<double>> matrix_sum;
    matrix_sum = vector<vector<double>>(n, vector<double>(n,0));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            matrix_sum[i][j] = A[i][j] + B[i][j];
        }   
    }
    return matrix_sum;
}

vector<vector<double>> M_theta_hA(vector<double> mlow, vector<double> mdiag, vector<double> mup, double sign)
{
    int n = mdiag.size();
    vector<vector<double>> diag, low, up, sum;
    diag = vector<vector<double>>(n, vector<double>(n,0));
    low = vector<vector<double>>(n, vector<double>(n,0));
    up = vector<vector<double>>(n, vector<double>(n,0));                    
    sum = vector<vector<double>>(n, vector<double>(n,0));
    
    for (int i = 0; i < n-1; i++)
    {
        mlow[i] *= sign;
        mdiag[i] *= sign;
        mup[i] *= sign;
    }
    mdiag[n-1] *= sign;
    
    diag = M_diag(mdiag);
    low = M_low(mlow);
    up = M_up(mup);

    sum = M_sum(diag,low);
    sum = M_sum(sum,up);

    return sum;
}

vector<double> M_mult_v(vector<vector<double>> A, vector<double> u)
{
    int n = u.size();
    vector<double> mul(n);
    for (int i = 0; i < n; i++)
    {
        for(int j = 0; j<n; j++)
        {
            mul[i] += A[i][j]*u[j]; 
        }
    }
    return mul;
}

vector<vector<double>> M_mult_M(vector<vector<double>> A, vector<vector<double>> B)
{
    int n = A.size();
    vector<vector<double>> m_mult;
    m_mult = vector<vector<double>>(n, vector<double>(n,0));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            for(int k = 0; k < n; ++k)
            {
                m_mult[i][j] += A[i][k] * B[k][j];
            }
    return m_mult;
}

vector<double> Decom_LU(vector<vector<double>> A, vector<double> b)
{
    int n = b.size();
    vector<double> y(n);
    vector<vector<double>> L, U;
    L =  vector<vector<double>>(n, vector<double>(n,0));
    U =  vector<vector<double>>(n, vector<double>(n,0));
    U[0][0] = A[0][0];
    y[0] = b[0]; 
    for (int i = 0; i <= n-1; i++)
    {
        L[i][i] = 1;
        if(i<n-1)
            U[i][i+1] = A[i][i+1]; 
        if(i>0)
            U[i][i] = A[i][i] - L[i][i-1] * A[i-1][i];
        if(i<n-1)
            L[i+1][i] = A[i+1][i] / U[i][i];
        if(i>0)
            y[i] = b[i] - L[i][i-1] * y[i-1];
    }
    vector<double> x(n,0);
    x[n-1] = y[n-1] / U[n-1][n-1];
    for (int j = n-2; j >= 0; j--)
    {
        x[j] = (y[j] - A[j][j+1] * x[j+1])/U[j][j];
    }
    
    return x;
}

vector<vector<double> > M_FD_Price(int option_type, double S0, double Maturity, double Strike, double r, double sigma, double divid, int nTime, int nSpace)
{
    double h = Maturity/nTime,
        b = r - divid - (sigma*sigma)/2.0,
        l = abs(b)*Maturity + 3.5*sigma*sqrt(Maturity),
        k = 2.0*l/nSpace;
    vector<double> diag(nSpace-1), up(nSpace-2), low(nSpace-2);
    vector<vector<double> > Price;
    Price = vector<vector<double> >(nSpace-1, vector<double>(nTime+1,0));
    
    for (int i =0; i < nSpace-2; i++)
    {
        low[i] = 1./4. * h * (i+2) * ((i+2)*sigma*sigma - r + divid);
        up[i] = 1./4. * h * (i+1) * ((i+1)*sigma*sigma + r - divid);
        diag[i] = -1./2. * h * ( (i+1)*(i+1)*sigma*sigma + r - divid);
    }
    diag[nSpace-2] = -1./2. * h * ( (nSpace-1)*(nSpace-1)*sigma*sigma + r - divid);
    
    // Neumann boundary condition
    diag[0] += 2. * 1./4. * h * (sigma*sigma - r + divid);
    diag[nSpace-2] += 2. * 1./4. * h * (nSpace-1) * ((nSpace-1)*sigma*sigma + r - divid);
    
    up[0] -= 1./4. * h * (sigma*sigma - r + divid);
    low[nSpace-3] -= 1./4. * h * (nSpace-1) * ((nSpace-1)*sigma*sigma + r - divid);
    
    vector<double> x(nSpace-1), Payoff(nSpace-1), result(nSpace-1) , ulte_Price(nSpace-1), id_vec(nSpace-1);
    
    for (int i = 0; i <= nSpace-2; i++)
    {
        x[i] = log(S0) - l + (i+1)*k;
        Price[i][nTime] = Phi(option_type,x[i],Strike);
        id_vec[i] = 1.0;
    }    
    vector<vector<double> > A,B1;
    A = vector<vector<double> >(nSpace-1, vector<double>(nSpace-1,0));
    B1 = vector<vector<double> >(nSpace-1, vector<double>(nSpace-1,0));
    A =  M_sum(M_diag(id_vec),M_theta_hA(low, diag, up, -1.0));
    B1 = M_sum(M_diag(id_vec),M_theta_hA(low, diag, up, 1.0));
    
    for(int i = nTime-1; i >= 0; i--)
    {
        for (int j = 0; j <= nSpace-2; j++)
        {
            ulte_Price[j] = Price[j][i+1];
        }
        
        result = Decom_LU(A,M_mult_v(B1,ulte_Price)); 
        
        for (int k = 0; k <= nSpace-2; k++)
        {
            Price[k][i] = result[k];
        }
    }
    
    return Price;
}

double FD_Price(int option_type, double S0, double Maturity, double Strike, double r, double sigma, double divid, int nTime, int nSpace)
{
    vector<vector<double>> P;
    P = M_FD_Price( option_type,  S0,  Maturity,  Strike,  r,  sigma,  divid,  nTime,  nSpace);
    if(nSpace % 2 == 0)
    {
        return P[(nSpace)/2 - 1][0];
    }
    else
    {
        return (P[(nSpace-1)/2][0] + P[(nSpace-3)/2][0])/2;
    }
    
}

int main()
{
    double price,s0, Maturity,  Strike, eqr,  eqdiv, sigma_ref, volimpl,impliedDelta,  impliedGamma, barrierup, barrierdown, rebate;
    price=.2;
    s0=5630.;
    Maturity=1.;
    Strike=5620.;
    rebate=0.;
    eqr=0.05;
    eqdiv=0.;
    sigma_ref=0.3;
    volimpl=0.3;
    barrierup=1.1;
    barrierdown=0.9;

    double nSpace, nTime;
    nSpace = 10;
    nTime = 12;
    vector<vector<double>> Matrix_FDPrice;
    Matrix_FDPrice = M_FD_Price(1, s0, Maturity, Strike, eqr, sigma_ref, eqdiv, nTime, nSpace);

    for(int i = 0; i<= nSpace-2; i++)
    {
        for(int j = 0; j<= nTime; j++)
        {
            cout << Matrix_FDPrice[i][j] << " ";
        }
        cout << " \n ";
    }
    
    cout << " \n ";
    
    cout << FD_Price(1, s0, Maturity, Strike, eqr, sigma_ref, eqdiv, nTime, nSpace);
    
    cout << " \n ";
    
    
}

