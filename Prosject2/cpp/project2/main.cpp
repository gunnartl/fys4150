#include <iostream>
#include <armadillo>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <iomanip>


using namespace std;
using namespace arma;

mat samurai_jacobi(double *A, double tol, int& output);
mat toeplitz(int N);
void check(int Nstart,int Nstop, float tol);
void potpot(int N);
void test_toeplitz();
void test_jacobi();

int main(){
    test_toeplitz();
    test_jacobi();
    //check(1,1,1e-6);
    potpot(400);
    return 0;
}

mat samurai_jacobi(mat A, double tol, int& count){
    //preforms jacobismethod to find eigenvalues
    double hugest,tau,t,t2,s,c;
    int n,k,l;
    mat S;
    count = 0; // count of iterations needed before biggest off-diagonal-elemt < til
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
    mat R = eye(n,n); // soon to bee eigenmatrix
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


        //S = eye(size(A)); // dustemetode
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
                a_ik = A(i,k); // lagrer ellementene så de ikke blir overskrevet
                a_il = A(i,l);
                A(i,k) = c*a_ik - s*a_il;
                A(k,i) = A(i,k); // symetri
                A(i,l) = c*a_il + s*a_ik;
                A(l,i) = A(i,l); // symetri
            }
            r_ik = R(i,k);  // oppdaterer egenverdier
            r_il = R(i,l);
            R(i,k) = c*r_ik - s*r_il;
            R(i,l) = c*r_il + s*r_ik;
        }

        

        hugest = 0; // new biggest element
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
    A = join_horiz(A,R); // returns eigenvalue- and eigenvector matrix contactanated
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
    // calculate the eigenvetors matrises of sises between 10*Nstart x 10*Nstart
    // and 10*Nstopx10*Nstop
    //also checks the error between eig_syms egenvalues and the ones produesed by JACOBI

    // program fails when looping over more than two values of N, task solved in python
    //still searching for the reason
    vec err, N,count,eigval;
    mat start,jacobi_vals, jacobi_vecs,eigvec;
    count = zeros(Nstart-Nstop+1);
    N = zeros(Nstart-Nstop+1);
    double h2;
    int counter;
    int n;
    for(int i = Nstart;i <=Nstop;i++){
        n = i*10;
        N[i] = n;
        h2 = 1./(n*n);
        start = toeplitz(n);
        eig_sym(eigval,eigvec,toeplitz(n)/h2);
        counter = 0;
        jacobi_vals = samurai_jacobi(start/h2,tol,counter);
        count[i] = counter;
        jacobi_vecs = jacobi_vals;
        jacobi_vecs.shed_cols(0,n-1);   //sheds the eigenvalues from the vectors
        jacobi_vals.shed_cols(n,n+n-1); //sheds the eigenvectors from thevalues
        jacobi_vals = sort(jacobi_vals.diag());
        cout<<jacobi_vals[0]<<endl;
        cout<<eigval[0]<<endl;
        cout<<abs((jacobi_vals[0]-eigval[0])/eigval[0])<<endl;
    }
}

void potpot(int N){

    double rho_0,rho_max,h;
    vec eigval;
    mat eigvec;
    rho_0 = 0;
    rho_max = 5;
    vec V = zeros(N);//linspace(rho_0,rho_max,N);
    h = (rho_max/N);
    for(int i = 0;i<N;i++){
        V[i] = ((i+1)*h)*((i+1)*h);
    }
    //vec V = rho%rho;


    mat A = toeplitz(N)/(h*h);
    cout<<V<<endl;
    A.diag() += V;
    //cout<<A<<endl;
    int counter;
    mat jacobi_vals, jacobi_vecs;
    jacobi_vals = samurai_jacobi(A,1e-6,counter); //matrix of eigenvalues and vectors contactenated
    jacobi_vecs = jacobi_vals;
    jacobi_vecs.shed_cols(0,N-1);   //sheds the eigenvalues from the vectors
    jacobi_vals.shed_cols(N,N+N-1); //sheds the eigenvectors from thevalues
    jacobi_vals = jacobi_vals.diag();
    cout<<jacobi_vals<<1234<<endl;
    uvec sortindex = sort_index(jacobi_vals);
    eig_sym(eigval,eigvec,A);
    mat egenvector0,egenvector1,egenvector2;
    mat egenval0,egenval1,egenval2;
    egenvector0 = jacobi_vecs.col(sortindex(0));
    egenvector1 = jacobi_vecs.col(sortindex(1));
    egenvector2 = jacobi_vecs.col(sortindex(2));
    mat eigenvectors;
    eigenvectors = join_horiz(egenvector0,join_horiz(egenvector1,egenvector2));
    egenval0 = jacobi_vals(sortindex(0));
    egenval1 = jacobi_vals(sortindex(1));
    egenval2 = jacobi_vals(sortindex(2));

    eigenvectors.save("eigenvectors0-2.txt",raw_ascii);
    cout<<counter<<setprecision(17)<<egenval0<<egenval1<<egenval2<<endl;



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
    int counter;
    mat s = samurai_jacobi(A,1e-6,counter);
    vec jacobi_eigval = sort(s.diag());
    for(int i= 0;i<5;i++){
        if (abs(jacobi_eigval(i)-analytical_eigval(i))>1e-10){
            cout<<"Error in function: JACOBI. Program aborted"<<endl;
            exit(1);
        }
    }
}
