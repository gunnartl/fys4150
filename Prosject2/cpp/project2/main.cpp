#include <iostream>
#include <armadillo>
#include <cmath>
#include <math.h>


using namespace std;
using namespace arma;

mat samurai_jacobi(double *A, double tol);
mat toeplitz(int N);
void check(int Nstart,int Nstop, float tol);
void potpot(int N);
void test_toeplitz();
void test_jacobi();

int main(int argc, char *argv[]){
    test_toeplitz();
    test_jacobi();
    //check(100,100,1e-6);
    //potpot(100);
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
    mat R = eye(n,n);
    while(abs(A(k,l))>tol){
        tau = (A(l,l)-A(k,k))/(2*A(k,l));
        //cout<<count<<endl;
        t   = -tau+sqrt(1+tau*tau);
        t2  = -tau-sqrt(1+tau*tau);
        if(abs(t2)<abs(t)){
            t = t2;
        }
        c   = 1/sqrt(1+t*t);
        s   = t*c;


        //S = eye(size(A));
        //S(k,k) = c;
        //S(l,l) = c;
        //S(k,l) = s;
        //S(l,k) = -s;
        //A = S.t()*A*S;

        double a_kk, a_ll, a_ik, a_il, r_ik, r_il;   //metode fra slides
        a_kk = A(k,k);
        a_ll = A(l,l);
        A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
        A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
        A(k,l) = 0.0;
        A(l,k) = 0.0;
        for ( int i = 0; i < n; i++ ) {
            if ( i != k && i != l ) {
                a_ik = A(i,k);
                a_il = A(i,l);
                A(i,k) = c*a_ik - s*a_il;
                A(k,i) = A(i,k);
                A(i,l) = c*a_il + s*a_ik;
                A(l,i) = A(i,l);
            }
            r_ik = R(i,k);
            r_il = R(i,l);
            R(i,k) = c*r_ik - s*r_il;
            R(i,l) = c*r_il + s*r_ik;
        }

        

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
    A = join_horiz(A,R);
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
    // calculate the eigenvetors matrises of sises between Nstart x Nstart
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
        //cout<<approx<<endl;
        cout<<eigval<<endl;
        cout<<2/h2-(2/h2)*cos(M_PI*1/(i+1))<<endl; //analytisk egenverdi nr x --> pi*x
        cout<<((sort((approx.diag()))));///sort(abs(eigval)))<<endl;
    }
}

void potpot(int N){

    double rho_0,rho_max,h;
    //double hbar = 1.0545718e-34;
    //double m = 9.10938356e-31;
    //double omega = 5;
    //double k =5 ;

    vec eigval;
    mat eigvec;
    rho_0 = 0;
    rho_max = 8.9;
    vec V = zeros(N);//linspace(rho_0,rho_max,N);
    h = (rho_max/N);
    for(int i = 0;i<N;i++){
        V[i] = ((i+1)*h)*((i+1)*h);
    }
    //vec V = rho%rho;


    mat A = toeplitz(N)/(h*h);
    cout<<V<<endl;
    A.diag() += V;
    cout<<A<<endl;
    mat jacobi_vals, jacobi_vecs;
    jacobi_vals = samurai_jacobi(A,1e-6); //matrix of eigenvalues and vectors contactenated
    jacobi_vecs = jacobi_vals;
    jacobi_vecs.shed_cols(0,N-1);   //sheds the eigenvalues from the vectors
    jacobi_vals.shed_cols(N,N+N-1); //sheds the eigenvectors from thevalues

    eig_sym(eigval,eigvec,A);
    cout<<sort(eigval)<<endl;
    //cout<<r<<endl;
    cout<<jacobi_vals.n_cols<<endl;
}

void test_toeplitz(){
    // Unit test for the toeplitz function
    // testing that toeplits returns Matrix A
    // for N = 4
    mat A,B;
    A= {{ 2,-1, 0, 0},
        {-1, 2,-1, 0},
        { 0,-1, 2,-1},
        { 0, 0,-1, 2}};
    B = toeplitz(4);
    for(int i = 0;i<4;i++){
        for(int j = 0;j<4;j++){
            if (A(i,j) != B(i,j)){
                cout<<"Error in function: toeplitz. Program aborted"<<endl;
                exit(1);
            }
        }
    }
}

void test_jacobi(){
    // Unit test for the jacobifunction.
    // checking that the eigenvalues and vectors returned
    //from the jacobi function for matrix A are the
    // same as our precalculated vectors and values within a tolerance = tol
    mat A;
    A = {{ 4, 1, 1, 0,-1},
         { 1, 4, 0, 1,-1},
         { 1, 0, 3, 0, 0},
         { 0, 1, 0, 3, 0},
         {-1,-1, 0, 0, 3}};
    vec analytical_eigval;
    analytical_eigval = {2,2,3,4,6};

    mat s = samurai_jacobi(A,1e-6);
    vec jacobi_eigval = sort(s.diag());
    for(int i= 0;i<5;i++){
        if (abs(jacobi_eigval(i)-analytical_eigval(i))>1e-10){
            cout<<"Error in function: JACOBI. Program aborted"<<endl;
            exit(1);
        }
    }
}
