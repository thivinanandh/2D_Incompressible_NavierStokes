#include<cstring>
#include<cmath>
#include<cstdlib>
#include <fstream>
#include<iostream>
#include"functions.h"

#include </home/thivin1/Softwares/gsl-1.0/gsl/gsl_blas.h>
#include </home/thivin1/Softwares/gsl-1.0/gsl/gsl_matrix.h>
#include </home/thivin1/Softwares/gsl-1.0/gsl/gsl_linalg.h>

using namespace std;



# define getName(var)  #var 

//////////////////////////////    INCOMPRESSIBLE NAVIER STOKES SOLVER          ////////////////////////////////////////
///////////// Coded for uniform mesh on Square domain with Dirichlet boundary conditions /////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*---------- Start of user input parameters ---------------------------------*/
#define time_step 0.001
#define grid_size 4
#define V_lid 1.0
#define RE_NR 100
#define max_iteration 200
#define EPSILON_ITERATION 1e-5
#define ITERATION_PER_RESULT 1
//#define PRINT 1       // YES to print console output for Debugging , NO to suppress output
double dirichlet_boundary_u  (int c )
{
    double a;
    switch(c)
    { 
        case 1: a = 0.0;       // Bottom parallel to X axis  i.e (y=0 ;  u(x,0))
            break;  
        case 2: a =  0.0;               //parallel to Y axis - right i.e (x = 1 ; u(1,y))
            break;
        case 3: a = V_lid;      // top parallel to X axis  i.e (y=1 ;  u(x,1))
            break;
        case 4: a =  0.0;              //parallel to Y axis - left i.e (x = 0 ; u(0,y))
            break;
        default: a =  0.0;
        break;
    }
    return a;
}


 double dirichlet_boundary_v  (int c)
{
    double a;
    switch(c)
    { 
        case 1: a = 0.0;       // Bottom parallel to X axis  i.e (y=0 ;  u(x,0))
            break;  
        case 2: a =  0.0;               //parallel to Y axis - right i.e (x = 1 ; u(1,y))
            break;
        case 3: a = 0.0;      // top parallel to X axis  i.e (y=1 ;  u(x,1))
            break;
        case 4: a =  0.0;              //parallel to Y axis - left i.e (x = 0 ; u(0,y))
            break;
        default: a =  0.0;
        break;
    }
    return a;
}

void set_boundary_u(double**& arr, int Nx , int Ny )
{
    
    // Boundary Conditions being appended to the uS Matrix ( Do not change Values here )
    for(int i = 0 ; i <Nx+1 ; i++ )
    {        // top
        arr[i][Ny+1] = arr[i][Ny+1] + dirichlet_boundary_u(3) ; 
        arr[i][Ny] = arr[i][Ny] + dirichlet_boundary_u(3) ;
        arr[i][0]    = arr[i][0] + dirichlet_boundary_u(1) ;
    }
    
    for(int j = 0 ; j <Nx+2 ; j++ )
    {        // top
        arr[0][j] = arr[0][j] + dirichlet_boundary_u(4) ; 
        arr[0][Nx] = arr[0][Nx] + dirichlet_boundary_u(2) ;
    }
}

void set_boundary_p(double**& arr, int Nx , int Ny )
{
    for (int i=0;i<Nx+2;i++)
        for(int j=0;j<Ny+2;j++)
            arr[i][j] = 1.0;
}

void set_boundary_v(double**& arr, int Nx , int Ny )
{
    // Boundary Conditions being appended to the uS Matrix ( Do not change Values here )
    for(int i = 0 ; i <Nx+1 ; i++ )
    {        // top
        arr[i][Ny+1] = arr[i][Ny+1] + dirichlet_boundary_v(3) ; 
        arr[i][Ny] = arr[i][Ny] + dirichlet_boundary_v(3) ;
        arr[i][0]    = arr[i][0] + dirichlet_boundary_v(1) ;
    }
    
    for(int j = 0 ; j <Nx+2 ; j++ )
    {        // top
        arr[0][j] = arr[0][j] + dirichlet_boundary_v(4) ; 
        arr[0][Nx] = arr[0][Nx] + dirichlet_boundary_v(2) ;
    }
}


void RHS_U_double_derivative(double**& arr, double** uS, int Nx , int Ny )
{
    // generates the RHS Adjustment that gors into the RHS of the system for it to be solved.
    for(int i = 1 ; i <Nx ; i++ )
        for (int j = 1 ; j < Ny ; j++)
        {
            if(j==Ny-1)  { arr[i-1][j-1] = arr[i-1][j-1] + uS[i][j+1] ; } // UXX - Top
            if(j==1)     { arr[i-1][j-1] = arr[i-1][j-1] + uS[i][j-1] ; } // UXX - Bottom
            if(i==1)     { arr[i-1][j-1] = arr[i-1][j-1] + uS[i-1][j] ; } // UYY - left
            if(i==Nx-1)  { arr[i-1][j-1] = arr[i-1][j-1] + uS[i+1][j] ; } // UYY - Right
        }
}

void RHS_V_double_derivative(double**& arr, double** vS, int Nx , int Ny )
{
    // generates the RHS Adjustment that gors into the RHS of the system for it to be solved.
    for(int i = 1 ; i <Nx ; i++ )
        for (int j = 1 ; j < Ny ; j++)
        {
            if(j==Ny-1)  { arr[i-1][j-1] = arr[i-1][j-1] + vS[i][j+1] ; } // UXX - Top
            if(j==1)     { arr[i-1][j-1] = arr[i-1][j-1] + vS[i][j-1] ; } // UXX - Bottom
            if(i==1)     { arr[i-1][j-1] = arr[i-1][j-1] + vS[i-1][j] ; } // UYY - left
            if(i==Nx-1)  { arr[i-1][j-1] = arr[i-1][j-1] + vS[i+1][j] ; } // UYY - Right
        }
}

void boundary_correction(double**& p , double**& uSnew , double**& vSnew,int Nx,int Ny)
{
    int i,j,k;   
     for(j=0;j<Nx+1;j++){
         p[0][j] = p[1][j];
         p[Nx+1][j] = p[Nx][j];
    }
//      p[0][1] = p[1][1] ; p[0][2] = p[1][2] ;
for(i=0;i<Nx+2;i++){
        p[i][0] = p[i][1];
        p[i][Nx+1] = p[i][Nx];
    }

// Correct the U Velocity
 
    for(i=0;i<Nx+1;i++)
        {
           uSnew[i][0] = -uSnew[i][1];
           uSnew[i][Nx+1] = (V_lid*2) - uSnew[i][Nx] ;   
        }

     for(j=0;j<Ny+1;j++)
        {
            vSnew[1][j] = -vSnew[0][j];
            vSnew[Nx+1][j] = -vSnew[Nx][j];
        }
}

void RHS_Pressure_poisson(double**& arr, int Nx, int Ny)
{
    for (int i = 1 ; i<Nx+1 ; i++ )
        for (int j = 1 ; j<Nx+1 ; j++ ){
            if(j==1 || j == Nx )
                arr[((i-1)*Nx) + (j-1)][((i-1)*Nx) + (j-1)] =  arr[((i-1)*Nx) + (j-1)][((i-1)*Nx) + (j-1)] + 1.0;
            
            if(i==1 || i == Nx)
                 arr[((i-1)*Nx) + (j-1)][((i-1)*Nx) + (j-1)] =  arr[((i-1)*Nx) + (j-1)][((i-1)*Nx) + (j-1)] + 1.0;
        }
}



    
int main()
{
    
    int Nx, Ny , i , j;
    // Calculating u.delu for the equation
    Nx = Ny = grid_size;
    double hx,hy,time = 0.0, L2_error_u = 0.0,L2_error_v = 0.0,error_norm = 0.0;
    int time_step_no = 1;
    double **u,**u_new, **uS ,**uSnew,**vSnew, **v ,**v_new, **vS , **p , **pS;
    memset_2d(u,Nx+1,Ny+1);memset_2d(uS,Nx+1,Ny+2);memset_2d(v,Nx+1,Ny+1);memset_2d(vS,Nx+2,Ny+1);
    memset_2d(p,Nx+1,Ny+1);memset_2d(v_new,Nx+1,Ny+1);memset_2d(u_new,Nx+1,Ny+1);
    memset_2d(uSnew,Nx+1,Ny+2);memset_2d(vSnew,Nx+2,Ny+1);memset_2d(pS,Nx+2,Ny+2);
    hx = 1/ double(Nx);
    hy = 1/ double(Ny);
    
    cout.flags( ios::dec | ios::scientific );
    cout.precision(4);
    
    //Set uS and  vS as per the initial boundary conditions
    //print_matrix(uSnew,Nx+1,Ny+2,getName(uSnew));
    /* 1, set boundary values for u and v*/
    set_boundary_u(uS,Nx,Ny);set_boundary_v(vS,Nx,Ny);set_boundary_p(pS,Nx,Ny);
    set_boundary_u(uSnew,Nx,Ny);set_boundary_v(vSnew,Nx,Ny);
    ///--------------------//
    
     //print_matrix(vS,Nx+2,Ny+1,getName(vS));
    
    double **uDelu , **vDelu , **RHS_D2U , **RHS_D2V , **RHS,**RHS2;
    memset_2d(uDelu,Nx-1,Ny-1);memset_2d(vDelu,Nx-1,Ny-1);memset_2d(RHS_D2U,Nx-1,Nx-1);
    memset_2d(RHS_D2V,Nx-1,Nx-1);memset_2d(RHS,Nx-1,Ny-1);memset_2d(RHS2,Nx-1,Ny-1);
//     double  **gradUstar , *gradUstar_vector, *intial_guess_x, **RHS_pressure, *p_vector;
//     memset_2d(gradUstar,Nx,Ny);memset_1d(gradUstar_vector,Nx*Ny);memset_2d(RHS_pressure,Nx*Ny,Nx*Ny);memset_1d(p_vector,Nx*Ny);memset_1d(intial_guess_x,Nx*Ny);
    
    //Calculate RHS value for U and V vectors which will be used in the calculations
    RHS_U_double_derivative(RHS_D2U,uS,Nx,Ny);
    RHS_V_double_derivative(RHS_D2V,vS,Nx,Ny);
    
    //Calculate the RHS value for Pressure poisson,. 
    // Since direct solution is not applicable for pressure poisson , we use BiCGSTAB to solve the system 
    //for which pressure is considered in a vector format.
    RHS_pressure = FD_2D_Second_derivative_matrix(Nx,Ny);
    RHS_Pressure_poisson(RHS_pressure,Nx,Ny);
    
    /////////////////////// ----------- GSL BLASSS  ----------------------- //////////////////////////////
    gsl_matrix *RHS_pressure_blas = gsl_matrix_alloc(Nx*Nx,Nx*Nx);
    gsl_vector *gradUstar_vector1 = gsl_vector_alloc(Nx*Nx);
    gsl_vector *p_vector1 = gsl_vector_alloc(Nx*Nx);
    gsl_vector *tau = gsl_vector_alloc(Nx*Nx);
        
    for(i=0;i<Nx*Nx;i++)
        for(j=0;j<Nx*Nx;j++)
            gsl_matrix_set(RHS_pressure_blas, i, j,RHS_pressure[i][j] );
    
    
    /////////////////////// ----------- GSL BLASSS  ----------------------- //////////////////////////////
    
  
    // Output Stream 
    ofstream outfile("Output001.dat", ios::out) ;
    outfile << "TITLE = Temp_distribution" <<endl;
    if(!outfile) {cerr<< "Error: Output file couldnot be opened.\n";}
    outfile.flags( ios::dec | ios::scientific);
    outfile.precision(4);
    
    // Declaring objects for classes to call eigen_parameters calss.
    eigen_parameters eig(Nx,Ny); 
    eigen_parameters eig1(Nx+1,Ny+1);
    
    double  alpha = RE_NR / time_step ;
    
    //------------------------------- Looping Starts ---------------------------------------------------------//
    error_norm = 1000;
    
    while((time_step_no-1) < max_iteration && error_norm > EPSILON_ITERATION)
    {
        
    time = time_step + time;
     // Helmholtz Constant for 1st equation
       
// Declaring object for eigen_parameters class.
    
    // calculate the Convective term for U term ,.//
    for( i=1 ; i < Nx ; i++ )
        for( j=1 ; j<Ny ; j++ ){
            uDelu[i-1][j-1]  = ( pow(uS[i+1][j],2) - pow(uS[i-1][j],2) ) / (2*hx) ;
            uDelu[i-1][j-1] += ( 0.5*0.5*( uS[i][j+1] + uS[i][j] ) *(vS[i+1][j] + vS[i][j]) ) / (hy);
            uDelu[i-1][j-1]  = uDelu[i-1][j-1] - ( ( 0.5*0.5* (uS[i][j] + uS[i][j-1] ) * ( vS[i+1][j-1] + vS[i][j-1] ) ) / hy);
        }
    
    // calculate the Convective term for V term ,.//
    for( i=1 ; i < Nx ; i++ )
        for( j=1 ; j<Ny ; j++ ){
            vDelu[i-1][j-1] = ( pow(vS[i][j+1],2) - pow(vS[i][j-1],2) ) / (2*hy) ;
            vDelu[i-1][j-1] += ( 0.5*0.5*( uS[i][j+1] + uS[i][j] ) *(vS[i+1][j] + vS[i][j]) ) / (hy);
            vDelu[i-1][j-1] = vDelu[i-1][j-1] - (0.5*0.5* (uS[i-1][j] + uS[i-1][j+1] ) * ( vS[i][j] + vS[i-1][j] ) /hx );
        }
#ifdef PRINT
    print_matrix(uDelu,Nx-1,Ny-1,getName(uDelu));
    print_matrix(vDelu,Nx-1,Ny-2,getName(vDelu));
    print_matrix(uS,Nx+1,Ny+2,getName(uS));
#endif
    
     for(i=0;i<Nx-1;i++)
        for(j=0;j<Nx-1;j++)
            RHS[i][j] = ( (hx*hy) * (alpha*uS[i+1][j+1] + (RE_NR *uDelu[i][j]) ) ) -  RHS_D2U[i][j]  ;
//      
    
    
    ///helmholtz_solver(double**& ans_arr,double** RHS,double* Eig_val,double** Eig_vec,double** Eig_vec_trans,int Nx,int Ny,double alpha)
    helmholtz_solver(uS,RHS,eig.Eig_val,eig.Eig_vec,eig.Eig_vec_trans,Nx,Ny,alpha);
    
    for(i=0;i<Nx-1;i++)
        for(j=0;j<Nx-1;j++)
            RHS2[i][j] = (hx*hy) * ( alpha*vS[i+1][j+1] + (RE_NR * vDelu[i][j]) ) -  RHS_D2V[i][j]  ;


    
    helmholtz_solver(vS,RHS2,eig.Eig_val,eig.Eig_vec,eig.Eig_vec_trans,Nx,Ny,alpha);
    
#ifdef PRINT
    print_matrix(RHS,Nx-1,Ny-1,getName(RHS));
    print_matrix(uS,Nx+1,Ny+2,getName(uS));
    print_matrix(RHS2,Nx-1,Ny-1,getName(RHS2));
    print_matrix(vS,Nx+1,Ny+2,getName(vS));
#endif
    
    // Calculate the Pressure Poisson equation. 

    
    // Calculate grad.u* which forms RHS of pressure poisson equation // 
    // Multiply - (h^2)/del(t) with it so that it directly forms the RHS of Helmholtz solver //

    for(i=1;i<Nx+1;i++)
        for(j=1;j<Ny+1;j++)
        {
         gradUstar[i-1][j-1] = ((uS[i][j] - uS[i-1][j])/hx) + ((vS[i][j] - vS[i][j-1])/(hy));
         gradUstar[i-1][j-1] *=  (hx*hy) / (time_step);
         gsl_vector_set(gradUstar_vector1, ((i-1)*Nx) + (j-1), ( ((uS[i][j] - uS[i-1][j])/hx) + ((vS[i][j] - vS[i][j-1])/(hy)))  * ( (hx*hy) / (time_step) )   );
        }

    cout<<endl;
    
    gsl_linalg_QR_decomp (RHS_pressure_blas,tau);
    gsl_linalg_QR_solve (RHS_pressure_blas,tau,gradUstar_vector1,p_vector1);
//    SYNTAX -  double* BiCGSTAB_Vector(double **A, double *b,double *x,int row, int col , double tolerance , int max_iter)
    
    for(i=1;i<Nx+1;i++)
        for(j=1;j<Ny+1;j++){
            pS[i][j] = gsl_vector_get(p_vector1, ( ((i-1)*Nx) + j-1 ) ) ;
            p_vector[((i-1)*Nx) + (j-1)] = 0.0;
        }
    
    // Solve Laplacian of P = (1/time_step)*div.U
    //helmholtz_solver(pS,gradUstar,eig1.Eig_val,eig1.Eig_vec,eig1.Eig_vec_trans,Nx+1,Ny+1,0);
    
    //print_matrix(pS,Nx+2,Ny+2,getName(pS));

#ifdef PRINT
    print_matrix(vS,Nx+2,Ny+1,getName(vS));
    print_matrix(gradUstar,Nx,Ny,getName(gradUstar));
    print_matrix(pS,Nx,Ny,getName(pS));
#endif
    
    //STEP 3 , Calulate final Velocity by including the change in pressure term.
    
    for(i=1;i<Nx;i++)
        for(j=1;j<Ny;j++)
        {
            uSnew[i][j] = uS[i][j] - (time_step/hx)*(pS[i+1][j] - pS[i][j]) ;
            vSnew[i][j] = vS[i][j] - (time_step/hy)*(pS[i][j+1] - pS[i][j]) ;
        }
    // Correct the pressure and Velocity values for boundary cells. 
    boundary_correction(pS,uSnew,vSnew,Nx,Ny);
#ifdef PRINT
    print_matrix(pS,Nx+2,Ny+2,getName(pS));
    print_matrix(vSnew,Nx+2,Ny+1,getName(vSnew));
    print_matrix(uSnew,Nx+1,Ny+2,getName(uSnew));
#endif
    
    // Calculate Normal Pressure and velocity value
    for(i=0;i<Nx+1;i++)
        for(j=0;j<Ny+1;j++){
            u_new[i][j] = (uSnew[i][j] + uSnew[i][j+1])/2;
        }
    for(i=0;i<Nx+1;i++)
        for(j=0;j<Ny+1;j++){
            v_new[i][j] = (vSnew[i][j] + vSnew[i+1][j])/2;
        }
    
    for(i=0;i<Nx+1;i++)
        for(j=0;j<Ny+1;j++){
            p[i][j] = (pS[i][j] + pS[i+1][j] + pS[i+1][j+1] + pS[i][j+1])/4;
        }
    
#ifdef PRINT
    print_matrix(u_new,Nx+1,Ny+1,getName(u_new));
    print_matrix(v_new,Nx+1,Ny+1,getName(v_new));
#endif
    // Print to print_to_output
    //void print_to_output(ofstream& outfile,double** u, double** uExact, double** v,double** vExact,double**p, double* x, double* y, double Nx, double Ny,int time_step_no , double time)
    if(time_step_no % ITERATION_PER_RESULT == 0 || time_step_no == 1)
    {
     print_to_output(outfile,u_new,u,v_new,v,p,eig.x,eig.y,Nx,Ny,time_step_no,time,error_norm);
    }
//     print_matrix(u,Nx+1,Ny+1,getName(u));
//     print_matrix(u_new,Nx+1,Ny+1,getName(u_new));
//     print_matrix(v_new,Nx+1,Ny+1,getName(v_new));
//     print_matrix(v,Nx+1,Ny+1,getName(v));
//     
    
    // Asssign the value to the old variable
    for(i=0;i<Nx+1;i++)
        for(j=0;j<Ny+2;j++)
            uS[i][j] = uSnew[i][j];
        
    for(i=0;i<Nx+2;i++)
        for(j=0;j<Ny+1;j++)
            vS[i][j] = vSnew[i][j];
        
    for(i=0;i<Nx+1;i++)
        for(j=0;j<Ny+1;j++){
            v[i][j] = v_new[i][j];
            u[i][j] = u_new[i][j];
        }
    //increement time_step_no
    time_step_no ++;
    
    
    
    }
    outfile.close() ; 
    
    
    ofstream outfile3("Central.dat", ios::out) ;
    outfile3 << "ZONE F=POINT" <<endl;
    outfile3 << "I = "<<grid_size <<endl;
    if(!outfile3) {cerr<< "Error: Output file couldnot be opened.\n";}
    outfile3.flags( ios::dec | ios::scientific);
    outfile3.precision(4);
    
    for(int j = 0 ; j <= Ny ; j++){
        double y_val = j*hy;
        outfile3 << y_val<<"\t" << 0.5* (u[(Nx/2)+1][j] + u[(Nx/2)][j]) << endl;}
    cout<<"................. Process Completed ..................." <<endl;
    
    delete u,uS,v,vS,p,pS,v_new,u_new,uSnew,vSnew,**uDelu , **vDelu , **RHS_D2U , **RHS_D2V , **RHS;
    
}
