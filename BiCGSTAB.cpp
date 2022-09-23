#include<iostream>
#include<cmath>
#include<cstring>
#include<fstream>
#include<ctime>
#include"functions.h"


using namespace std;

# define getName(var)  #var 




double** BiCGSTAB(double **A, double **b,double **x,int row, int col , double tolerance , int max_iter)
{
    /* Algorithm taken from SIAM Paper :
    BI-CGSTAB: A FAST AND SMOOTHLY CONVERGING VARIANT OF BI-CG FOR THE SOLUTION OF NONSYMMETRIC LINEAR SYSTEMS* H. A. VAN DER VORSTt 
    2d matrrix details
    A = Input matrix - Square and Real
    b = Result matrix - suare and Real
    x =  intial guess matrix
    tolerance =  error tolerance for convergence check
    max_iter = maxximum Iteration;
    
    Algorithm used:
        r0 = b − Ax0
        Choose an arbitrary vector r̂0 such that (r̂0, r0) ≠ 0, e.g., r̂0 = r0 . Note that notation (x,y) applies for scalar product of vectors (x,y) = <x,y> = x·y = x ' y
        ρ0 = α = ω0 = 1
        v0 = p0 = 0
        For i = 1, 2, 3, …
        ρi = (r̂0, ri−1)
        β = (ρi/ρi−1)(α/ωi−1)
        pi = ri−1 + β(pi−1 − ωi−1vi−1)
        vi = Api
        α = ρi/(r̂0, vi)
        h = xi−1 + αpi
        If h is accurate enough, then set xi = h and quit
        s = ri−1 − αvi
        t = As
        ωi = (t, s)/(t, t)
        xi = h + ωis
        If xi is accurate enough, then quit
        ri = s − ωit
    */
    
   
    
    
    int i,j;
    int i_deviation,j_deviation;
    clock_t BiCGSTAB_start, BiCGSTAB_end;
    // Declaration 
    double **r0,**r0_hat;
    double **p0,**p0_new,**v,**v_new, **H,**x_new, **S, **T, **Ax;
    double rho,alpha,omega, rho_new,omega_new,error,beta;
    
    memset_2d(r0,row,col);memset_2d(r0_hat,row,col);
    memset_2d(p0,row,col);memset_2d(p0_new,row,col);memset_2d(v,row,col);memset_2d(v_new,row,col);
    memset_2d(H,row,col);memset_2d(x_new,row,col);memset_2d(S,row,col);memset_2d(T,row,col);memset_2d(Ax,row,col);
    
    double B_norm = norm_2(b,row);
    
    
    // Stp -1 ; Calulate R0 , Assign it to R0_new
    // Calculate Ax
    Ax = matrix_multiply(A,x,row,col,col,0,0);

    for (i = 0; i< row ; i++){
        for(j=0;j< col ; j++){
            r0[i][j] = b[i][j] - Ax[i][j];
            r0_hat[i][j] = r0[i][j];
        }
    }
    
    // decalatre values of Constants - Alpha , omega 
    alpha = rho = omega = rho= 1.0;
    
    // Set values of v0 and p0 as 0 Matrix, this would have been  done during declaration itself. 
    BiCGSTAB_start = clock();
    int iteration_count = 1;
    double L2_error = 1000;
    while ( (fabs(L2_error) > tolerance) && (iteration_count  <= max_iter))
    {
        //ρi = (r̂0, ri−1)
        rho_new = 0.0;
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                rho_new += r0_hat[i][j] * r0[i][j];
        
         if( iteration_count  > 1 )
         {
            // β = (ρi/ρi−1)(α/ωi−1)
            beta = (rho_new/rho) * (alpha/omega);
            //Pi = ri−1 + β(Pi−1 − ωi−1Vi−1)
            for (i = 0;i<row;i++)
                for (j=0;j<col;j++)
                    p0_new[i][j] = r0[i][j] + beta*(p0[i][j] -  (omega*v[i][j])) ;
         }
         else
             p0_new  = r0;
        
        // vi = APi
        v_new = matrix_multiply(A,p0_new,row,col,row,0,0);
               
        //α = ρi/(r̂0, vi)clc
        double temp  = 0.0;
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                temp += r0_hat[i][j] * v_new[i][j];
        alpha = rho_new/temp;
        
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                H[i][j] = x[i][j] + alpha*p0_new[i][j];
        
        // Early Exit Criteria , if the Norm of (alpha * p0 ) is less than tolerance , then exit the loop 
        L2_error = 0.0;
        L2_error = (alpha*norm_2(p0_new,row)/B_norm);
        if( (iteration_count > 1) && ( fabs(L2_error) < tolerance))
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has converged ********** " << endl;
            cout<<"Time Elapsed : "<< double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC<<"  ms "<< endl;
            cout << "Iteration :  " << iteration_count << ".5 " <<endl<<"Error :" <<  L2_error << endl;
            cout<<"Tolerance Provided : " << tolerance <<endl;
            //print_matrix(H,row,col,getName(h));
            return H;
        }

        // s = ri−1 − αvi
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                S[i][j] = r0[i][j] -  alpha*v_new[i][j];
        //print_matrix(S,row,col,getName(S));
        // t = As
        T = matrix_multiply(A,S,row,col,col,0,0);
        //print_matrix(T,row,col,getName(T));
        
        //ωi = (t, s)/(t, t)
        double temp2 = 0.0;
        double temp1 = 0.0;
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++){
                temp1+=T[i][j] * S[i][j];
                temp2+=T[i][j] * T[i][j];
        }
        
        omega = temp1/temp2;
        
        // xi = h + ωis
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                x_new[i][j] = H[i][j] + omega*S[i][j];
        
        //print_matrix(x_new,row,col,getName(x_new));
        //ri = s − ωit
         for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                r0[i][j] = S[i][j] - omega*T[i][j];
        
        //print_matrix(r0,row,col,getName(r0));
        L2_error = 0.0;
        L2_error = (norm_2(r0,row)/B_norm);
        if( (iteration_count > 1) && fabs(L2_error) < tolerance)
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has converged **********" << endl;
            cout<<"Time Elapsed : "<< double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC<<"  ms "<< endl;
            cout << "Iteration :  " << iteration_count   <<endl<<" Error :  "<< L2_error << endl;
            cout<<"Tolerance Provided : " << tolerance <<endl;
            return x_new;
        }
        
        // Error Exception handling
        if( omega == 0.0 || rho == 0.0 )
            cout <<" BICGSTAB has not converged , it may be due to singularity in the given 'A' matrix " << endl;
        
        // to print the outout for each 100 Iterations 
        if(iteration_count % 100 == 0 )
            cout<<" BICGSTAB - at Iteration " << iteration_count << " Error is : "<<L2_error<< endl;
        
        // assign x to x-new and rho_new to rho
        rho = rho_new;
        x = x_new;
        iteration_count++;
        //cout<<"   ----- Iteration ended ---- " << endl;
    
    }
    BiCGSTAB_end = clock();
    cout<< " BICGSTAB Iteration has ****NOT**** converged to the Solution "<< endl;
    cout<<" Time Elapsed : "<< double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC<<"  ms "<< endl;
    cout << " Iteration :  " << iteration_count   <<endl<< " Error :  "<< L2_error << endl;
    
    return x_new;
    
    // delete the variables
    delete r0,r0_hat, p0,p0_new,v,v_new, H,x_new, S, T, Ax;
//     double rho,alpha,omega, rho_new,omega_new,error,beta;
    
    
}


int main()
{
    // Set the precision flags for cout
    cout.flags( ios::dec | ios::scientific );
    cout.precision(5);
    
    int Nx = 32;
    int Ny = 32;
    int i,j;
    
    double **A, **b , **x , **c, **d;
    memset_2d(A,Nx,Ny);memset_2d(b,Nx,Ny);memset_2d(x,Nx,Ny);memset_2d(c,Nx,Ny);
    memset_2d(d,Nx,Ny);
    
    for(i=0;i<Nx;i++)
        for(j=0;j<Ny;j++){
            A[i][j] = rand()%200;
            b[i][j] = (j+1)/(i+1);
            x[i][j] = 0.0;
        
        }

    

    c = BiCGSTAB(A,b,x,Nx,Ny,1e-12,50000);
    
    d = matrix_multiply(A,c,Nx,Nx,Ny,0,0);

    //print_matrix(c,Nx,Ny,getName(c));
    //print_matrix(d,Nx,Ny,getName(d));
    bool flag = true;
    
    for(i=0;i<Nx;i++)
        for(j=0;j<Ny;j++){
            if( fabs(d[i][j] - b[i][j] ) > 0.1 )
                flag = false;
    }
        
    if( flag )
        cout<< " PASSED -- "<< endl;
    else
        cout <<" - FAILED "<<endl;
    //print_matrix(x,2,2,getName(x));
    return 0;
}

