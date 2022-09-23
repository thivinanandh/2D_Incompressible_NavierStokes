#include<cstring>
#include<cmath>
#include<cstdlib>
#include <fstream>
#include<iostream>
#include"functions.h"
#include"omp.h"


using namespace std;


// To get the Variable name to be printed while debugging
# define getName(var)  #var 

//////////////////////////////    INCOMPRESSIBLE NAVIER STOKES SOLVER          ////////////////////////////////////////
///////////// Coded for uniform mesh on Square domain with Dirichlet boundary conditions         /////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*---------------- Start of User INput  Values ---------------------------------*/
// #define time_step 0.00005
// #define grid_size 16
 #define V_lid 1.0
// #define RE_NR 100
 //#define max_iteration 20000
// #define EPSILON_ITERATION 1e-8
// #define ITERATION_PER_RESULT 100
//#define PRINT 1       // YES to print console output for Debugging , NO to suppress output
//#define PARALLEL 1
//#define NUM_THREADS 4


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

//-------------------------- END OF USER INPUT VARRIABLES ---------------------------------------//


// Sets the FD matrix based in the dirichlet boundary conditions provided by the user in above code for U , V , P//
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
// *** END***//

// Calculates the RHS ( g(x) )  value for the Double Derivative Matrix  for including them in the RHS of the solving system //

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

// ***** END ****** //

// This function to append the U,V,P of the Boundary cells of the Staggered grid , based on the conditions , i.e
// in out case , the U at bottom 2 of setup has to be such that , their Sum is Zero. 
// for pressure , their difference has to be zero //.

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
           uSnew[i][Nx+1] = uSnew[i][Nx] = V_lid;   
        }
    
// Correct the V Velocity
     for(j=0;j<Ny+1;j++)
        {
            vSnew[1][j] = -vSnew[0][j];
            vSnew[Nx+1][j] = -vSnew[Nx][j];
        }
}

// ***** END ****** //


//// START OF MAIN FUNCTION //////

int main(int n, char **argv)
{
    
    double RE_NR = stod(argv[1]);
    double time_step = stod(argv[2]);
    int grid_size = stoi(argv[3]);
    int max_iteration = stoi(argv[4]);
    int ITERATION_PER_RESULT = stoi(argv[5]);
    double EPSILON_ITERATION = stod(argv[6]);
    
    cout <<"   Time Step : " << time_step<<endl;
    cout <<"   RE NO : " << RE_NR<<endl;
    cout <<"   Max_iter : " << max_iteration<<endl;
    cout<<"   Iteration per result : "<< ITERATION_PER_RESULT<<endl;
    cout<<"   grid Size : "<<grid_size<<endl;
    cout<<"   EPSILON_ITERATION : "<< EPSILON_ITERATION<<endl;
    #ifdef PARALLEL
    cout<< " ---- PARALLELISED --- Threads : " << NUM_THREADS <<endl;
    #endif
    cout <<" ------------------------------------------------------- " <<endl;
    #ifdef PARALLEL
    omp_set_num_threads (NUM_THREADS);
    #endif
    // Declaring the Variables of the process
    int Nx, Ny , i , j;
    Nx = grid_size;
     Ny = grid_size;
    double hx,hy,time = 0.0, L2_error_u = 0.0,L2_error_v = 0.0,error_norm = 0.0;
    int time_step_no = 1;
    double **u,**u_new, **uS ,**uSnew,**vSnew, **v ,**v_new, **vS , **p , **pS;
    double **uDelu , **vDelu , **RHS_D2U , **RHS_D2V , **RHS,**RHS2,**gradUstar ;
//     double *gradUstar_vector,*p_vector ,*x_int_guess;
    memset_2d(u,Nx+1,Ny+1);memset_2d(uS,Nx+1,Ny+2);memset_2d(v,Nx+1,Ny+1);memset_2d(vS,Nx+2,Ny+1);
    memset_2d(p,Nx+1,Ny+1);memset_2d(v_new,Nx+1,Ny+1);memset_2d(u_new,Nx+1,Ny+1);
    memset_2d(uSnew,Nx+1,Ny+2);memset_2d(vSnew,Nx+2,Ny+1);memset_2d(pS,Nx+2,Ny+2);
   // memset_1d(gradUstar_vector,Nx*Nx);memset_1d(p_vector,Nx*Nx);memset_1d(x_int_guess,Nx*Nx);
    memset_2d(uDelu,Nx-1,Ny-1);memset_2d(vDelu,Nx-1,Ny-1);memset_2d(RHS_D2U,Nx-1,Nx-1);
    memset_2d(RHS_D2V,Nx-1,Nx-1);memset_2d(RHS,Nx-1,Ny-1);memset_2d(gradUstar,Nx,Nx);
    hx = 1/ double(Nx);
    hy = 1/ double(Ny);

    //Setting up flags for COUT Stream for precision printing // 
    cout.flags( ios::dec | ios::scientific );
    cout.precision(4);
    

    
    //Set u Staggered, vStaggered, P Staggered as per the initial boundary conditions given by the user 
    /* 1, set boundary values for u and v*/
    set_boundary_u(uS,Nx,Ny);set_boundary_v(vS,Nx,Ny);set_boundary_p(pS,Nx,Ny);

    
    //Calculate RHS value for ∇2U and ∇2V Matrices which will be added to the RHS of the system
    RHS_U_double_derivative(RHS_D2U,uS,Nx,Ny);
    RHS_V_double_derivative(RHS_D2V,vS,Nx,Ny);
   
    // Declaring objects for classes to call eigen_parameters class .This Holds the values such as 
    // grid_structure , Eigen values ,eigen_vectors for a given grid size double derivative matirx , which will be used to 
    // solve the "Helmholtz part" of the Navier Stokes //
    eigen_parameters eig(Nx,Ny); 
    eigen_parameters eig1(Nx+1,Ny+1);
    
    
    // Output Stream parameters for prining to file
    ofstream outfile("Output001.dat", ios::out) ;
    outfile << "TITLE = Temp_distribution" <<endl;
    if(!outfile) {cerr<< "Error: Output file couldnot be opened.\n";}
    outfile.flags( ios::dec | ios::scientific);
    outfile.precision(4);
    
    // USed in fining U* and V*
    double  alpha = double(RE_NR) / double(time_step) ;
    
    //-----------------------------------Time Looping Starts ---------------------------------------------------------//
    error_norm = 1000;
    
    while( (time_step_no-1) < max_iteration  && error_norm > EPSILON_ITERATION)
    {
        
    time = time_step + time;
      
    // calculate the Convective term for U term ,.//
    // d(uu)/ dx + d(uv)/dy 
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif 
    for( i=1 ; i < Nx ; i++ )
        for(int j=1 ; j<Ny ; j++ ){
            uDelu[i-1][j-1]  = ( pow(uS[i+1][j],2) - pow(uS[i-1][j],2) ) / (2*hx) ;
            uDelu[i-1][j-1] += ( 0.5*0.5*( uS[i][j+1] + uS[i][j] ) *(vS[i+1][j] + vS[i][j]) ) / (hy);
            uDelu[i-1][j-1]  = uDelu[i-1][j-1] - ( ( 0.5*0.5* (uS[i][j] + uS[i][j-1] ) * ( vS[i+1][j-1] + vS[i][j-1] ) ) / hy);
        }
    
    // calculate the Convective term for V term ,.//
    // d(uv)/ dx + d(vv)/dy
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif 
    for( i=1 ; i < Nx ; i++ )
        for(int j=1 ; j<Ny ; j++ ){
            vDelu[i-1][j-1] = ( pow(vS[i][j+1],2) - pow(vS[i][j-1],2) ) / (2*hy) ;
            vDelu[i-1][j-1] += ( 0.5*0.5*( uS[i][j+1] + uS[i][j] ) *(vS[i+1][j] + vS[i][j]) ) / (hy);
            vDelu[i-1][j-1] = vDelu[i-1][j-1] - (0.5*0.5* (uS[i-1][j] + uS[i-1][j+1] ) * ( vS[i][j] + vS[i-1][j] ) /hx );
        }
    #ifdef PRINT  // Print for debugging purposes
    print_matrix(uDelu,Nx-1,Ny-1,getName(uDelu));
    print_matrix(vDelu,Nx-1,Ny-2,getName(vDelu));
    print_matrix(uS,Nx+1,Ny+2,getName(uS));
    #endif
    // --------------- Solve for U* at nth time step ---------------------------//
    // Solve for U* at the nth Time step by Ignoring the pressure gradient part. 
    // [ ∇^2 - alpha*(h^2) ]U*   = h^2 * ( Re* (d(uu)/ dx + d(uv)/dy)) -  RHS_D2U ( boundary terms from  ∇^2U* part )
    
    // Construct the RHS Matrix for "U" part //
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif 
    for( i=0;i<Nx-1;i++)
        for(int j=0;j<Nx-1;j++)
            RHS[i][j] = ( (hx*hy) * ((RE_NR *uDelu[i][j]) - alpha*uS[i+1][j+1]) ) -  RHS_D2U[i][j]  ;
        
    // For debugging
    #ifdef PRINT
    print_matrix(uS,Nx+1,Ny+2,getName(uS));
    print_matrix(RHS,Nx-1,Ny-1,getName(RHS2));  
    #endif
    
    // Solve for U* usng the custom written function Helmholtz matrix solver , which uses the Eigen based solver for double derivative matrix //
    // ** SYNTAX**helmholtz_solver( Answer array  ,RHS MAtrix , Eig Value vector ,Eig_vector matrix ,Eig_vec_transpose matrix ,Nx,Ny, alpha < Helmholtz constant>)
    //  All values must be sent as Pointers , and the ans array should be declated in main funciton and will be updated by the function //
    helmholtz_solver(uS,RHS,eig.Eig_val,eig.Eig_vec,eig.Eig_vec_trans,Nx,Ny,alpha);
    
    // Construct the RHS Matrix for "U" part //
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif 
    for( i=0;i<Nx-1;i++)
        for(int j=0;j<Nx-1;j++)
            RHS[i][j] = (hx*hy) * ((RE_NR * vDelu[i][j]) -  alpha*vS[i+1][j+1]) -  RHS_D2V[i][j]  ;

    // Solve for V* usng the custom written function Helmholtz matrix solver 
    helmholtz_solver(vS,RHS,eig.Eig_val,eig.Eig_vec,eig.Eig_vec_trans,Nx,Ny,alpha);
    
    // Debugging 
    #ifdef PRINT
    print_matrix(RHS,Nx-1,Ny-1,getName(RHS));
    print_matrix(vS,Nx+2,Ny+1,getName(vS));
    #endif
    
    // ------- END of alculation of U* and V* ---------------------//
    
    // ------------ Solve Pressure POisssion ----------------------//
    // Solve the Pressure Poisson equation, which calculates both the pressure at Nth step and satisfies Continuity equation //
    // (∇^2)P  = (1/∆t) (∇U*)
    
    // Calculate grad.u* which forms RHS while solving pressure poisson equation // 
    // Multiply - (h^2)/del(t) with it so that it directly forms the RHS of Helmholtz solver //
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif 
    for( i=1;i<Nx+1;i++)
        for(int j=1;j<Ny+1;j++)
         gradUstar[i-1][j-1] =  ( ((uS[i][j] - uS[i-1][j])/hx) + ((vS[i][j] - vS[i][j-1])/(hy)) )  / ( (hx*hy) / (time_step) )  ;        

    // Solve (∇^2)P = (1/time_step)(∇U*)
    helmholtz_solver(pS,gradUstar,eig1.Eig_val,eig1.Eig_vec,eig1.Eig_vec_trans,Nx+1,Ny+1,0);
    
    // Debugging part
    #ifdef PRINT
    //print_matrix(vS,Nx+2,Ny+1,getName(vS));
    //print_matrix(gradUstar,Nx,Ny,getName(gradUstar));
    print_matrix(pS,Nx,Ny,getName(pS));
    #endif
    
    // -------END of Pressure Poisson Solver -------//
    
    //STEP 3 , Calulate Un+1 and  Vn+1 by including the Gradient of Pressure Term //
    // U(n+1) = U* - ∆t*∇P(n)
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif 
    for( i=1;i<Nx;i++)
        for(int j=1;j<Ny;j++)
        {
            uSnew[i][j] = uS[i][j] - (time_step/hx)*(pS[i+1][j] - pS[i][j]) ;
            vSnew[i][j] = vS[i][j] - (time_step/hy)*(pS[i][j+1] - pS[i][j]) ;
        }
    // Correct the pressure and Velocity values for boundary cells so that the end conditions remains satisfied
    boundary_correction(pS,uSnew,vSnew,Nx,Ny);
    
    // Debugging
    #ifdef PRINT
    print_matrix(pS,Nx+2,Ny+2,getName(pS));
    print_matrix(vSnew,Nx+2,Ny+1,getName(vSnew));
    print_matrix(uSnew,Nx+1,Ny+2,getName(uSnew));
    #endif
    // ***** END of Staggered Grid Calculation ******//
    
    
    // Calculate Pressure and Velocity in the normal grid by Interpolating between the staggered grid_size
    // U = average of 2 in y direction ; V = aerage of 2 in X direction ; P = Average of $ pressure cells in all direction 
    // U INterpolation
    
    for(i=0;i<Nx+1;i++)
        for(j=0;j<Ny+1;j++){
            u_new[i][j] = (uSnew[i][j] + uSnew[i][j+1])/2;
        }
    // V INterpolation
    for(i=0;i<Nx+1;i++)
        for(j=0;j<Ny+1;j++){
            v_new[i][j] = (vSnew[i][j] + vSnew[i+1][j])/2;
        }
    // P INterpolation
    for(i=0;i<Nx+1;i++)
        for(j=0;j<Ny+1;j++){
            p[i][j] = (pS[i][j] + pS[i+1][j] + pS[i+1][j+1] + pS[i][j+1])/4;
        }
    
    #ifdef PRINT
    print_matrix(u_new,Nx+1,Ny+1,getName(u_new));
    print_matrix(v_new,Nx+1,Ny+1,getName(v_new));
    #endif
    
    //***** END Of U,V,P Interpolation ******//
    
    // Print to print_to_output and Console using custom function

    
    double Maximum_error_u = 0.0,Maximum_error_v = 0.0;
    double  L2_error_u = 0.0, L2_error_v = 0.0;
        for(i = 0 ; i <= Nx; i++)
        for(j = 0 ; j <= Ny ; j++)
        {
            if( fabs(u[i][j] - u_new[i][j])  > Maximum_error_u)
                Maximum_error_u = fabs(u[i][j] - u_new[i][j]) ;
            
            if( fabs(v[i][j] - v_new[i][j])  > Maximum_error_v)
                Maximum_error_v = fabs(v[i][j] - v_new[i][j]) ;

            L2_error_u = L2_error_u + ( ((u[i][j] - u_new[i][j])*(u[i][j] - u_new[i][j])) / ((Nx+1)*(Ny+1)) );
            
            L2_error_v = L2_error_v + ( ((v[i][j] - v_new[i][j])*(v[i][j] - v_new[i][j])) / ((Nx+1)*(Ny+1)) );
        }
    L2_error_u = sqrt(L2_error_u);
    L2_error_v = sqrt(L2_error_v);
    
    error_norm = 0.0;
    error_norm = sqrt( pow(L2_error_u,2) + pow(L2_error_v,2) ) ;

    if(time_step_no % ITERATION_PER_RESULT == 0 || time_step_no == 1)
    {
        print_to_output(outfile,u_new,u,v_new,v,p,eig.x,eig.y,Nx,Ny,time_step_no,time,error_norm);
    }
    
    // Allocate New values to the old values for error calculation//
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
    
    // *********************************** END of TIME LOOPING ****************************************//
    outfile.close() ; 
    
    // U vs Y plot generation for benmark Comparsion //
    ofstream outfile3("Central.plt", ios::out) ;
    outfile3<<"VARIABLES=\"U\",\"Y\" "<<endl;
    outfile3 << "ZONE F=POINT" <<endl;
    outfile3 << "I="<<grid_size <<endl;
    if(!outfile3) {cerr<< "Error: Output file couldnot be opened.\n";}
    outfile3.flags( ios::dec | ios::scientific);
    outfile3.precision(4);
    
    for(int j = 0 ; j <= Ny ; j++){
        double y_val = j*hy;
        outfile3 << 0.5 * (u[(Nx/2)+1][j] + u[(Nx/2)][j]) <<"\t" << y_val << endl;}

    cout<<"................. Process Completed ..................." <<endl;
    outfile3.close() ; 
    delete u,uS,v,vS,p,pS,v_new,u_new,uSnew,vSnew,**uDelu , **vDelu , **RHS_D2U , **RHS_D2V , **RHS , **gradUstar;
}
