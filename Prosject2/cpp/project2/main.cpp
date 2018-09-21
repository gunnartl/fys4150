#include <iostream>
#include <armadillo>
#include <cmath>


using namespace std;
using namespace arma;

mat samurai_jacobi(mat A, double tol);
mat toeplitz(int N);

int main(int argc, char *argv[]){
    int N = 100;
    mat A = toeplitz(N);
    mat s = samurai_jacobi(A,1e-6);
    cout<<s;
    return 0;
}

mat samurai_jacobi(mat A, double tol){
    double hugest,tau,t,t2,s,c;
    int count,n,k,l;
    mat S;
    count = 0;
    n = A.n_cols;
    k = 0;
    l = 1;
    hugest = abs(A(k,l));
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
        if(t2<t){
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
    mat A = zeros(N,N);
    A.diag().ones();
    A.diag(1).ones();
    A.diag(-1).ones();
    A.diag()*=2;
    A.diag(1)*=-1;
    A.diag(-1)*=-1;
    return A;
}
