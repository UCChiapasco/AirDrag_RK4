//
//  chiapasco_project2.cpp
//  Project2
//
//  Created by Umberto Carlo Chiapasco
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

#define DEBUG 0

const double g_vel_i = 10.0; // the total speed of the projectile is
// initially assumed equal (in module) to 10 m/s
const double g_xf = 5.0; // the target is the point (5.0m,1.0m)
const double g_yf = 0.0;
const double g_grav = 9.81;
// note: the bullet mass is assumed equal to 1 kg


void RK4Step(double, double *, void(*)(double, double *, double *),
  double, int);
void dYdt(double, double *, double *);
int bracketing(double(*)(double),double, double, int, int&,
  double*, double*);
double Residual(double);
int bisection_rf(double(*)(double), double, double, double,
  int&, double, double&);

int main(){
  double a = 0., b = M_PI*0.5, tol = 1e-7, xtol = 0.01;
  // note that:
  // tol is used in bisection as the width of the interval around the
  // root
  // xtol is the tolerance for the distance between the x-result at
  // the chosen theta and the actual target
  int N = 10, nroots, nbisec, op;
  double sr[64], sl[64], root[64];
  // part 1: identify the roots of Residual()
  op = bracketing(Residual, a, b, N, nroots, sl, sr);
  if(op != 0) { cout << "Error occurred on bracketing()\n"; }
  #if DEBUG
  cout << "Bracketing executed. Number of roots: " << nroots << endl;
  #endif
  for(int i=0; i < nroots; i++){
    op = bisection_rf(Residual,sl[i],sr[i],tol,nbisec,xtol,root[i]);
    if(op != 0) { cout << "Error occurred on bisection_rf()\n"; }
    #if DEBUG
    cout << "Bisection: root n. " << i+1 << ": " << root[i] << endl;
    #endif
  }

  // part 2: check solutions and save trajectory data
  double ytol = 0.001;
  // ytol is the tolerance on the y-axis, useful to determine
  // the proximity to the ground in the 2nd cycle
  int nstepsmax = 10000, neq = 4, n;
  double rb = 0.0, r = rb, dr = 0.01, E; //dr = 0.1*g_xf/g_vel_i;
  double v_iy, halfway;
  double Y[4];
  ofstream fdata;
  fdata.open("chiapasco_project2_data.dat");
  // loop over the roots of Residual() found in part 1:
  for(int i=0; i < nroots; i++){
    n = 0;
    r = rb;
    v_iy = g_vel_i*sin(root[i]);
    Y[0] = 0.0; Y[1] = 0.0; // initial point
    Y[2] = g_vel_i*cos(root[i]);   // initial x speed
    Y[3] = v_iy;                   // initial y speed
    // instead of using the condition on the velocity Y[3] one could
    // use a condition on the x (Y[0] < halfway), with:
    // halfway really as half way (on x) point:
    //// halfway = g_xf*0.5;
    // halfway point using kinematics & initial conditions:
    //// halfway=(Y[2]+sqrt(Y[2]*Y[2]-2*g_grav*g_yf))*0.5*Y[3]/g_grav;
    // halfway point as approximated top point of the trajectory:
    //// halfway = Y[2]*Y[3]/g_grav;
    #if DEBUG
    cout << "I cicle\n";
    #endif
    /////////////////////////////////////////////////////
    // Compute the result obtained with RK4 shooting:
    // The loop has been split in order to avoid requiring the
    // condition on Y in the first part of the trajectory.
    // 1st half of the trajectory:
    while(n < nstepsmax && Y[3] > ytol) {
      // save trajectory data on "chiapasco_project2_data.dat"
      E = 0.5*Y[2]*Y[2] + 0.5*Y[3]*Y[3] + g_grav*Y[1];
      fdata << r << "  " << Y[0] << "  " << Y[1] << "  " << Y[2];
      fdata << "  " << Y[3] << "  " << E << endl;
      RK4Step(r,Y,dYdt,dr,neq);
      r+=dr;
      #if DEBUG
      cout << "*";
      #endif
      n++;
    }
    /////////////////////////////////////////////////////
    // 2nd half of the trajectory:
    dr /= 5.; // increase precision when getting near to the target
    #if DEBUG
    cout << "\nII cicle\n";
    #endif
    while(n < nstepsmax && Y[1] - g_yf > ytol) {
      // the condition requires to be at the same height as
      // the target in order to stop the trajectory
      E = 0.5*Y[2]*Y[2] + 0.5*Y[3]*Y[3] + g_grav*Y[1];
      fdata << r << "  " << Y[0] << "  " << Y[1] << "  " << Y[2];
      fdata << "  " << Y[3] << "  " << E << endl;
      RK4Step(r,Y,dYdt,dr,neq);
      r+=dr;
      n++;
      #if DEBUG
      cout << "*";
      #endif
    }
    E = 0.5*Y[2]*Y[2] + 0.5*Y[3]*Y[3] + g_grav*Y[1];
    fdata << r << "  " << Y[0] << "  " << Y[1] << "  " << Y[2];
    fdata << "  " << Y[3] << "  " << E << endl;
    /////////////////////////////////////////////////////
    #if DEBUG
    cout << endl;
    if(n >= nstepsmax){ cout << "nstepsmax reached.\n"; }
    else if(Y[1] - g_yf < ytol){cout << "mass reached the ground.\n";}
    else { cout << "an error occurred.\n"; }
    #endif
    fdata << endl << endl;
    cout << "Theta value " << i + 1 << ": " << root[i];
    cout << ".\n  x result: " << Y[0] << ", target: " << g_xf;
    cout << ".\n  Difference: " << fabs(Y[0] - g_xf);
    if(fabs(Y[0] - g_xf) < xtol){
      cout << " has acceptable deviation\n";
    }
    else{
      cout << " is beyond the required tolerance on x: ";
      cout << xtol << "\n";
    }
  }
  fdata.close();
  return 0;
}

void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *,
  double *), double h, int Neq){
  double Y1[256], k1[256], k2[256], k3[256], k4[256];
  int i;
  RHSFunc(t,Y,k1);
  for (i=0;i<Neq;i++) Y1[i] = Y[i] + 0.5*h*k1[i];
  RHSFunc(t+0.5*h,Y1,k2);
  for (i=0;i<Neq;i++) Y1[i] = Y[i] + 0.5*h*k2[i];
  RHSFunc(t+0.5*h,Y1,k3);
  for (i=0;i<Neq;i++) Y1[i] = Y[i] + h*k3[i];
  RHSFunc(t+h,Y1,k4);
  for (i=0;i<Neq;i++) Y[i] += h/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
}

int bisection_rf(double (*F_fr)(double), double a, double b,
  double tol, int &n, double ytol, double &root){
  // on return, bisection_rf returns:
  // 0: no errors ocurred
  // 1: initial values problems
  // 2: result problems
  if(a>b){
    cout << "! bisection_rf: Error: b must be greater than a\n";
    return 1;
  }
  if(tol<0){
    cout << "! bisection_rf: Error: tol must be greater than 0\n";
    return 1;
  }
  double f_l, f_r;
  f_l = F_fr(a);
  f_r = F_fr(b);
  
  if(f_l*f_r >= 0){
    cout << "! bisection_rf: Error: interval must have positive";
    cout << " & negative extremes in order to use this function;\n";
    return 1;
  }
  double x_mid,f_mid;
  int i;
  n = log2(1./tol);
  #if DEBUG
  cout << "f_l: " << f_l << " || f_r: " << f_r << endl;
  #endif
  for(i=0; i<n; i++){
    x_mid = (a+b)/2.0;
    f_mid = F_fr(x_mid);
    #if DEBUG
    cout << "x_mid: " << x_mid << " || f_mid: " << f_mid << endl;
    #endif
    if(f_l*f_mid>=0){a = x_mid;}
    else{b = x_mid;}
  }
  root = a;
  if(F_fr(root) >= ytol){
    cout << "! bisection_rf: An error occurred, the function did ";
    cout << "not converge to zero. Try again with lower tolerance.\n";
    return 2;
  }
  return 0;
}

int bracketing(double (*F_fr)(double), double a, double b,
  int N_intervals, int &n_roots, double *xL_arr,double *xR_arr){
  // on return, bracketing returns:
  // 0: no errors ocurred
  // 1: initial values problems
  // 2: result problems

  if(a>b){
    cout << "! bracketing: Error: b must be greater than a\n";
    return 1;
  }
  int i;
  double len = (b - a)/N_intervals, f_l,f_r;
  n_roots = 0;
  f_r = F_fr(a);
  for(i=1;i<=N_intervals;i++){
    f_l = f_r;
    f_r = F_fr(a + i*len);
    if(f_l*f_r<0){ // => the function changes sign in this interval
      *xL_arr = a + (i-1.0)*len;
      *xR_arr = a + i*len;
      xL_arr++;
      xR_arr++;
      n_roots++;
    }
  }
  return 0;
}

void dYdt(double t, double *Y, double *R){
  double B = 4.e-5;
  R[0] = Y[2];
  R[1] = Y[3];
  R[2] = -B*Y[2]*sqrt(Y[2]*Y[2] + Y[3]*Y[3]);
  R[3] = -g_grav - B*Y[3]*sqrt(Y[2]*Y[2] + Y[3]*Y[3]);
}

double Residual(double theta){
  // Residual function returns the error committed on the
  // shooting at angle theta
  // Note that the only residual (and therefore tolerance)
  // we consider is that on the x result. The height is only
  // useful in order to calculate the resulting x range.
  // note that a lot of lines in this function are commented.
  // This is mainly useful if one wants to debug this function
  // separately from the rest of the code
  double ytol = 0.001;
  // ytol has the same meaning as in main()
  int nstepsmax = 10000, neq = 4, n = 0, E;
  double rb = 0.0, r = rb, dr = 0.01; //dr = 0.02*g_xf/g_vel_i;
  double v_iy = g_vel_i*sin(theta), halfway;
  double Y[4];
  Y[0] = 0.0; Y[1] = 0.0; // initial point
  Y[2] = g_vel_i*cos(theta);
  Y[3] = v_iy;
  // halfway can be chosen in different ways (see main()). We
  // choose it as approximed top point of the trajectory;
  // halfway = Y[2]*Y[3]/g_grav;
  // #if DEBUG
  // cout << "I cicle\n";
  // #endif
  // ofstream fdata;
  // fdata.open("chiapasco_project2_data.dat");
  /////////////////////////////////////////////////////
  // 1st half of the trajectory:
  while(n < nstepsmax && Y[3] > ytol) {
    // save data on file to show the trajectory:
    // E = 0.5*Y[2]*Y[2] + 0.5*Y[3]*Y[3] + g_grav*Y[1];
    // fdata << r << "  " << Y[0] << "  " << Y[1] << "  " << Y[2];
    // fdata << "  " << Y[3] << "  " << E << endl;
    RK4Step(r,Y,dYdt,dr,neq);
    r+=dr;
    // #if DEBUG
    // cout << "*";
    // #endif
    n++;
  }
  /////////////////////////////////////////////////////
  // 2nd half of the trajectory:
  dr /= 5.; // increase precision when getting near to the target
  // #if DEBUG
  // cout << "\nII cicle\n";
  // #endif
  while(n < nstepsmax && Y[1] - g_yf > ytol) {
    // save data on file to show the trajectory:
    // E = 0.5*Y[2]*Y[2] + 0.5*Y[3]*Y[3];
    // fdata << r << "  " << Y[0] << "  " << Y[1] << "  " << Y[2];
    // fdata << "  " << Y[3] << "  " << E << endl;
    RK4Step(r,Y,dYdt,dr,neq);
    r+=dr;
    n++;
    // #if DEBUG
    // cout << "*";
    // #endif
  }
  // #if DEBUG
  // cout << endl;
  // if(n >= nstepsmax){ cout << "nstepsmax reached.\n"; }
  // else if(Y[1] - g_yf < ytol){ cout << "mass reached the ground.\n"; }
  // else { cout << "an error occurred.\n"; }
  // #endif
  // fdata.close();
  return Y[0] - g_xf;
}
