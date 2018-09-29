#include <iostream>
#include <armadillo>
#include <cmath>
#include <math.h>


using namespace std;
using namespace arma;

mat samurai_jacobi(mat A, double tol);
mat toeplitz(int N);
void check(int Nstart,int Nstop, float tol);
void potpot(int N);

int main(int argc, char *argv[]){
    //check(10,10,1e-6);
    potpot(100);
    return 0;
}

mat samurai_jacobi(mat A, double tol){
    //preforms jacobismethod to find eigenvalues
    double hugest,tau,t,t2,s,c;
    int count,n,k,l;
    mat S;
    count = 0;
    n = A.n_cols;
    k = 0;
    l = 1;
    hugest = abs(A(k,l)); //biggest off-diagonal element
    for(int i = 0;i<n;i++){
        for(int j= i+1;j<n;j++){
            if(abs(A(i,j))>hugest){
                hugest = abs(A(i,j));
                k = i;
                l = j;
            }
        }
    }
    while(abs(A(k,l))>tol){
        tau = (A(l,l)-A(k,k))/(2*A(k,l));
        cout<<count<<endl;
        t   = -tau+sqrt(1+tau*tau);
        t2  = -tau-sqrt(1+tau*tau);
        if(abs(t2)<abs(t)){
            t = t2;
        }
        c   = 1/sqrt(1+t*t);
        s   = t*c;
        S = eye(size(A));
        S(k,k) = c;
        S(l,l) = c;
        S(k,l) = s;
        S(l,k) = -s;
        A = S.t()*A*S;
        hugest = 0;
        for(int i = 0;i<n;i++){
            for(int j= i+1;j<n;j++){
                if(abs(A(i,j))>hugest){
                    hugest = abs(A(i,j));
                    k = i;
                    l = j;
                }
            }
        }
        count +=1;
    }
    return A;
}

mat toeplitz(int N){
    //making a tridiagonal matrix with -1 ones on the upper
    //and lower diagonal and 2 on the diagonal. size NxN
    mat A = zeros(N,N);
    A.diag().ones();
    A.diag(1).ones();
    A.diag(-1).ones();
    A.diag()*=2;
    A.diag(1)*=-1;
    A.diag(-1)*=-1;
    return A;
}

void check(int Nstart,int Nstop,float tol){
    // checking the number of iterations needed to
    // calculate the eigenvetors of sises between Nstart x Nstart
    // and NstopxNstop
    //also checks the error between eig_syms egenvalues and the ones produesed by JACOBI
    vec err, N,count,eigval;
    mat start,approx,eigvec;
    double h2;
    for(int i = Nstart;i <=Nstop;i+=10){
        h2 = 1./(i*i);
        start = toeplitz(i);
        eig_sym(eigval,eigvec,toeplitz(i)/h2);
        approx = samurai_jacobi(start,tol)/h2;
        cout<<approx<<endl;
        cout<<eigval<<endl;
        cout<<2/h2-(2/h2)*cos(M_PI*5/(i+1))<<endl; //analytisk egenverdi nr x --> pi*x
        cout<<((sort(abs(approx.diag()))-sort(abs(eigval))));///sort(abs(eigval)))<<endl;
    }
}

void potpot(int N){

    double rho_0,rho_max,h2;
    double hbar = 1.0545718e-34;
    double m = 9.10938356e-31;
    double omega = 5;
    double k =5 ;

    vec eigval;
    mat eigvec;
    rho_0 = 0;
    rho_max = 10;
    vec rho = linspace(rho_0,rho_max,N);
    h2 = (rho_max/N-1)*(rho_max/N-1);
    vec V = rho%rho;
    mat A = toeplitz(N)/h2;
    cout<<A<<endl;
    A.diag() += V;
    //cout<<A<<endl;
    //mat s = samurai_jacobi(A,1e-6);
    eig_sym(eigval,eigvec,toeplitz(N)/h2);
    cout<<eigval<<endl;
    //cout<<sort(s.diag());
}
