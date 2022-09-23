#include<cstring>
#include<cmath>
#include<cstdlib>
#include <fstream>
#include<iostream>

// Set the value of pi
double const PI = 4.0*atan(1.0);

using namespace std;

/*---------- End of user input parameters ---------------------------------*/

void memset_2d( double**& arr, int row , int col)
{
    arr = new double*[row] ;
    for (int i = 0 ; i < row ; i ++ )
    {
        arr[i] = new double[col];
        for (int j = 0 ; j< col ; j++ )
            arr[i][j] = 0.0;
    }
}



void memset_1d( double*& arr, int row)
{
    arr = new double[row] ;
    for (int i = 0 ; i < row ; i ++ )
    {
        arr[i] = 0.0;
    }
}

double vector_inner_dot(double* vec1 , double* vec2,int Nx)
{
    double ans = 0.0;
    for(int i = 0 ; i< Nx ; i++)
        ans += vec1[i]*vec2[i];
    return ans;
}

void matrix_vector_multiply(double*& ans_arr,double** arr, double* vec, int size)
{
    for ( int i =0 ; i < size ; i++){
        ans_arr[i] = 0;
        for( int j =0 ; j< size ; j++)
            ans_arr[i] += arr[i][j] * vec[j] ; 
    }
}


void print_matrix(double** arr1,int row , int col , string var_name)
{
    
    std::cout<< var_name << " : " ;
    
    for(int i=0;i<row;i++)
    {
        cout<<endl;
        for (int j=0;j<col;j++){
            cout<<arr1[i][j]<<"\t";
        }
    }
    cout<<endl<<endl;;
}

void print_array(double* arr1,int row, string var_name)
{
    std::cout<< var_name << " : " ;
    for(int i=0;i<row;i++)
            cout<<arr1[i]<<"\t";
    cout<<endl<<endl;;
}

// void matrix_multiply( double**& ans_arr, double** arr_a , double** arr_b, int row1, int col1,int col2, int i_dev, int j_dev)
// {
//     
//     int i,j,k;
//     for( i=0 + i_dev ;i < row1-i_dev;i++)
//         for ( k=0 + j_dev ; k<col1 - j_dev ;k++)
//             for ( j=0 + j_dev ;j < col2 - j_dev ;j++)    
//                 ans_arr[i-i_dev][j-j_dev]  += (arr_a[i][k] * arr_b[k][j]);
// }

void matrix_multiply( double**& ans_arr, double** arr_a , double**
arr_b, int row1, int col1,int col2,int i_dev,int j_dev)
{
omp_set_num_threads (4);
int block_size = 16;
    int i,j,k;
    for (i = 0 ; i < (row1) ; i += block_size)
        {
            for (j = 0 ; j < (col2 ); j += block_size)
            {
                #pragma omp parallel for collapse(2)
                for (int x = 0; x < block_size; ++x)
                {
                    for (int y = 0; y < block_size; ++y)
                    {
                        for (int k = 0; k < (col1); ++k)
                        {
                         #pragma omp critical
                            ans_arr[(i) + x][(j) + y] += arr_a[i + x][k] * arr_b[k][j + y];
                        }
                    }
                }
            }
        }
}

void helmholtz_solver(double**& arr,double** RHS,double* Eig_val,double** Eig_vec,double** Eig_vec_trans,int Nx,int Ny,double alpha)
{
    int i,j,k;
    double hx = 1/ double(Nx);
    double hy = 1/ double(Ny);
    double **temp1 , **F_Tilde , **temp2, **U_Tilde,**temp3;
    memset_2d(F_Tilde,Nx-1,Ny-1);memset_2d(U_Tilde,Nx-1,Ny-1);
    memset_2d(temp1,Nx-1,Ny-1); memset_2d(temp2,Nx-1,Ny-1);memset_2d(temp3,Nx-1,Ny-1);
    
    matrix_multiply(temp1,Eig_vec_trans,RHS,Nx-1,Ny-1,Nx-1,0,0);
    matrix_multiply(F_Tilde,temp1,Eig_vec,Nx-1, Ny-1,Nx-1,0,0);
    
    for(i=0;i<Nx-1;i++)
        for(j=0;j<Ny-1;j++)
            U_Tilde[i][j]  = F_Tilde[i][j] / (Eig_val[i] + Eig_val[j] - ((alpha) *hx*hy)) ;
    
    matrix_multiply(temp2,Eig_vec,U_Tilde,Nx-1,Ny-1,Nx-1,0,0);
    
    for(i = 1 ; i < Nx ; i++)
        for (j = 1; j < Ny ; j++){
            arr[i][j] = 0.0;
            for (int k = 1; k < Nx ; k++)
                arr[i][j] += (temp2[i-1][k-1] * Eig_vec_trans[k-1][j-1]);   
        }

    // Free the memory
    //Free the array of pointers
    
    for(int i = 0; i < Nx-1; i++) {
        delete [] temp1[i];  
        delete [] F_Tilde[i];
        delete [] temp2[i];
        delete [] U_Tilde[i];
        delete [] temp3[i];
        //delete [] Eig_vec_trans[i];
        
    }
    delete [] temp1;delete [] F_Tilde; delete [] temp2; delete [] U_Tilde;delete [] temp3;//delete [] Eig_vec_trans;
    
    
}

void sherman_morrison_solver(double*& ans, double* yy, double* zz,double* n,int Nx)
{
    double *rhs;
    double nty = 0.0, ntz = 0.0; 
    memset_1d(rhs,Nx+1); 
    
    nty = vector_inner_dot(n,yy,Nx);
    ntz = vector_inner_dot(n,zz,Nx);
    ntz += 1.0;
    for(int i = 0 ; i<Nx+1 ; i++){
        rhs[i] = 0.0;
        rhs[i] = (nty/ntz) * zz[i];
        ans[i] = yy[i] - rhs[i];
    }
    delete rhs;
}


void print_to_output(ofstream& outfile,double** u, double** uExact, double** v,double** vExact,double**p, double* x, double* y, double Nx, double Ny,int time_step_no , double time, double& error_norm)
{
    double Maximum_error_u = 0.0,Maximum_error_v = 0.0;
    double  L2_error_u = 0.0, L2_error_v = 0.0;
   
    outfile<< " VARIABLES = \"x\", \"y\",\"U\",\"V\",\"P\" " <<endl;
    outfile << "Zone T = psi I = " << Ny+1 << " J = " << Nx+1 << endl ;
    outfile << " StrandID =1, SolutionTime = "<<time_step_no/200 <<endl;
    for(int i = 0 ; i <= Nx; i++)
        for(int j = 0 ; j <= Ny ; j++)
        {
            if( fabs(u[i][j] - uExact[i][j])  > Maximum_error_u)
                Maximum_error_u = fabs(u[i][j] - uExact[i][j]) ;
            
            if( fabs(v[i][j] - vExact[i][j])  > Maximum_error_v)
                Maximum_error_v = fabs(v[i][j] - vExact[i][j]) ;

            L2_error_u = L2_error_u + ( ((u[i][j] - uExact[i][j])*(u[i][j] - uExact[i][j])) / ((Nx+1)*(Ny+1)) );
            
            L2_error_v = L2_error_v + ( ((v[i][j] - vExact[i][j])*(v[i][j] - vExact[i][j])) / ((Nx+1)*(Ny+1)) );
            
            outfile <<x[i] <<"\t"<<y[j] <<"\t"<< u[i][j] <<"\t"<<v[i][j]<<"\t"<<p[i][j]<<endl;
        }
    outfile<<endl<<endl;
    L2_error_u = sqrt(L2_error_u);
    L2_error_v = sqrt(L2_error_v);
    
    error_norm = 0.0;
    error_norm = sqrt( pow(L2_error_u,2) + pow(L2_error_v,2) ) ;
    
    cout <<"Iteration : " << time_step_no<<"   Time : "<<time<<endl;
    cout<<"   L2_error(U) : " << L2_error_u <<"     L-inf error(U) : " <<Maximum_error_u << endl;
    cout<<"   L2_error(V) : " << L2_error_v <<"     L-inf error(V) : " <<Maximum_error_v << endl;
    cout<<"   L(2) Error Norm : " << error_norm <<endl;
}

double norm_2_matrix(double** b , int n)
{
    double norm = 0.0;
    for (int i =0 ; i< n ; i ++ )
        for(int j = 0 ; j < n ; j++)
            norm += (b[i][j]*b[i][j])/(n*n);
        
    return sqrt(norm);
    
}

double norm_2_vector(double* b , int n)
{
    double norm = 0.0;
    for (int i =0 ; i< n ; i ++ )
            norm += (b[i]*b[i])/(n);
        
    return sqrt(norm);
    
}

double** BiCGSTAB_Matrix(double **A, double **b,double **x,int row, int col , double tolerance , int max_iter)
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
    
    double B_norm = norm_2_matrix(b,row);
    
    
    // Stp -1 ; Calulate R0 , Assign it to R0_new
    // Calculate Ax
     matrix_multiply(Ax,A,x,row,col,col,0,0);


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
        matrix_multiply(v_new,A,p0_new,row,col,row,0,0);
               
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
        L2_error = (alpha*norm_2_matrix(p0_new,row)/B_norm);
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
        matrix_multiply(T,A,S,row,col,col,0,0);
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
        L2_error = (norm_2_matrix(r0,row)/B_norm);
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

void BiCGSTAB_Vector(double*& x_new,double **A, double *b,double *x,int row, int col , double tolerance , int max_iter)
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
    double *r0,*r0_hat;
    double *p0,*p0_new,*v,*v_new, *H, *S, *T, *Ax;
    double rho,alpha,omega, rho_new,omega_new,error,beta;
    
    memset_1d(r0,col);memset_1d(r0_hat,col);
    memset_1d(p0,col);memset_1d(p0_new,col);memset_1d(v,col);memset_1d(v_new,col);
    memset_1d(H,col);memset_1d(x_new,col);memset_1d(S,col);memset_1d(T,col);memset_1d(Ax,col);
    
    double B_norm = norm_2_vector(b,row);
    
    
    // Stp -1 ; Calulate R0 , Assign it to R0_new
    // Calculate Ax
    matrix_vector_multiply(Ax,A,x,col);

         
    for (i = 0; i< col ; i++){
            r0[i] = b[i] - Ax[i];
            r0_hat[i] = r0[i];
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
        rho_new = vector_inner_dot(r0_hat,r0,col);
        
         if( iteration_count  > 1 )
         {
            // β = (ρi/ρi−1)(α/ωi−1)
            beta = (rho_new/rho) * (alpha/omega);
            
            //Pi = ri−1 + β(Pi−1 − ωi−1Vi−1)
            for (i = 0;i<row;i++)
                    p0_new[i] = r0[i] + beta*(p0[i] -  (omega*v[i])) ;
         }
         else
             p0_new  = r0;
        
        // vi = APi
        matrix_vector_multiply(v_new,A,p0_new,col);
        
       
        //α = ρi/(r̂0, vi)clc
        double temp  = 0.0;
        temp =  vector_inner_dot(r0_hat,v_new,col);
        
        alpha = rho_new/temp;
        
        for (i = 0;i<col;i++)
                H[i] = x[i] + alpha*p0_new[i];
        
                // s = ri−1 − αvi
        for (i = 0;i<col;i++)
                S[i] = r0[i]-  alpha*v_new[i];
        
        // Early Exit Criteria , if the Norm of (S ) is less than tolerance , then exit the loop 
        L2_error = 0.0;
        L2_error = (norm_2_vector(S,col));
        if( ( fabs(L2_error) < tolerance))
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has converged in "<< iteration_count << " in " << double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC << "time ********** " << endl;
            //print_matrix(H,row,col,getName(h));
            return ;
        }
        if( ( fabs(L2_error) > 1e20))
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has NOT converged , Error blown up , Matrix might be Singular "<< endl;
            //print_matrix(H,row,col,getName(h));
            return ;
        }

        //print_matrix(S,row,col,getName(S));
        // t = As
         matrix_vector_multiply(T,A,S,col);
        //print_matrix(T,row,col,getName(T));

        //ωi = (t, s)/(t, t)
        double temp2 = 0.0;
        double temp1 = 0.0;
        temp1 = vector_inner_dot(T,S,col);
        temp2 = vector_inner_dot(T,T,col);
       
        omega = temp1/temp2;
        
        // xi = h + ωis
        for (i = 0;i<row;i++)
                x_new[i] = H[i] + omega*S[i];


        //print_matrix(x_new,row,col,getName(x_new));
        //ri = s − ωit
         for (i = 0;i<row;i++)
                r0[i] = S[i] - omega*T[i];
        
        //print_matrix(r0,row,col,getName(r0));
        L2_error = 0.0;
        L2_error = (norm_2_vector(r0,col)/B_norm);
        
        if( (iteration_count > 1) && fabs(L2_error) < tolerance)
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has converged in "<< iteration_count << " in " << double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC << "time ********** " << endl;
            return ;
        }
        
        // Error Exception handling
        if( omega == 0.0 || rho == 0.0 )
            cout <<" BICGSTAB has not converged , it may be due to singularity in the given 'A' matrix " << endl;
        
        // to print the outout for each 100 Iterations 
        if(iteration_count % 10000 == 0 )
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
    
    return ;
    
    // delete the variables
    delete r0,r0_hat, p0,p0_new,v,v_new, H, S, T, Ax;
}

double** FD_2D_Second_derivative_matrix(int Nx , int Ny)
{
    double **arr;
    double** diagonal;double** identity;double** zero;
    memset_2d(arr,Nx*Ny,Nx*Ny);
    memset_2d(diagonal,Nx,Ny);memset_2d(identity,Nx,Ny);memset_2d(zero,Nx,Ny);
    
    int block_x_id;int block_y_id;int block_x;int block_y;
       
    // fills the block matrices Individually ( A - Diagonal , B - Identity , C - Zero matrix )
    for (int i=0;i<Nx;i++)
        for(int j=0;j<Nx;j++){
            if(i==j){
                diagonal[i][j] = -4.0;
                identity[i][j] = 1.0;
                zero[i][j] = 0.0;
            }
            else if(abs(i-j) == 1){
                diagonal[i][j] = 1.0;
                identity[i][j] = 0.0;
                zero[i][j] = 0.0;
            }
            else{
                diagonal[i][j] = 0.0;
                identity[i][j] = 0.0;
                zero[i][j] = 0.0;
            }
        }
        
    // Now fils the Blocks of the overall  matrix  based on the block id's of the matrix 
    for ( int block_x_id = 0; block_x_id < Nx ; block_x_id++)
    {
        for(int block_y_id = 0 ; block_y_id<Ny; block_y_id++)
        {
            if(block_x_id == block_y_id)
            {
                // Diagonal Block - diagonal matrix
                for(int i = 0 ; i < Nx ; i++)
                    for(int j=0;j<Ny;j++)
                        arr[(block_y_id*Ny) + i][(block_x_id*Nx) + j] = diagonal[i][j];
                    
            }
            else if(abs(block_x_id - block_y_id) == 1)
            {
                 // Sub Diagonal - Identity
                for(int i = 0 ; i < Nx ; i++)
                    for(int j=0;j<Ny;j++){ 
                        
                        arr[(block_y_id*Ny) + i][(block_x_id*Nx) + j] = identity[i][j];
                    }
            }
            
            else
            {
                for(int i = 0 ; i < Nx ; i++)
                    for(int j=0;j<Ny;j++)
                        arr[(block_y_id*Ny) + i][(block_x_id*Nx) + j] = zero[i][j];
            }
        }
    }
    
    return arr;
    delete arr,diagonal,identity,zero;
}


// Data Structures //

class eigen_parameters
{
public:
    double *x ,*y, *Eig_val,*norm;
    double **Eig_vec,**Eig_vec_trans;
    int Nx,Ny;
    eigen_parameters(int NX,int NY)
    {
        Nx = NX;Ny = NY;
        int i,j,k;
        memset_2d(Eig_vec,(Nx-1),(Ny-1) ); 
        memset_2d(Eig_vec_trans,(Nx-1),(Ny-1) ); 
        x = new double[Nx+1]; 
        y = new double[Ny+1];
        Eig_val = new double[Nx-1];
        norm = new double[Nx-1];
        
        for ( int i = 0; i <= Nx ; i++)  x[i] = i/double(Nx);   
        for ( int i = 0; i <= Ny ; i++)  y[i] = i/double(Ny);
        for( i = 1 ; i <Nx ; i++ )
        {
            Eig_val[i-1] = -4.0*sin(PI*0.5*x[i])*sin(PI*0.5*x[i]) ;
            norm[i-1] = 0.0;
            for ( j = 1 ; j < Ny ; j++)
            {
                Eig_vec[i-1][j-1] = sin(i*PI*x[j]);  
                norm[i-1] = norm[i-1]  + (Eig_vec[i-1][j-1]*Eig_vec[i-1][j-1]) ;     
            }
        }
        //Normalising Eigen Vectors 
        for( i = 1 ; i<Nx ; i++)
            for(j = 1;j<Ny;j++)
                Eig_vec[i-1][j-1] = Eig_vec[i-1][j-1]/(sqrt(norm[i-1]));
        // Eigen Vector Inverse
        for( i = 0 ; i<Nx-1 ; i++)
            for(j = 0;j<Ny-1;j++)
            Eig_vec_trans[i][j] = Eig_vec[j][i]  ;      
    }
   
   
   
   ~eigen_parameters() {
           for(int i = 0; i < Nx-1; i++) {
        delete [] Eig_vec[i];  
        delete [] Eig_vec_trans[i];
        //delete [] Eig_vec_trans[i];
        
    }
    delete [] Eig_vec;delete [] Eig_vec_trans;
};
};


