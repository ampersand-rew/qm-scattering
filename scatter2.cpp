#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include "TText.h"
#include "TMultiGraph.h"
#include "TGraph.h"

#include <gsl/gsl_sf_bessel.h>

using namespace std;

const double hbar = 1.054571817e-34; //reduced planck's constant, J*s
const double m = 9.1093837015e-31; //electron mass, kg
const double q = 1.602176634e-19; //elementary charge, Coulombs
const double eV = q; //Electronvolt, J
const double MeV = 1e6 * eV; //Mega Electron Volt, J
const double meV = 1e-3 * eV; //Milli Electron Volt, J
const double rB = 5.29e-11; //Bohr radius, meters
const double A = 1e-10; //Angstrom, meters

const double epsilon = 5.9; //epsilon = 5.9 meV
const double rho = 3.57; //ro = 3.57 angstrom
const double alpha = 6.12;

//Harmonic oscillator effective potential.
double FHO(int l, double r, double E) {
  if (l > 0) {
    return r*r + l*(l+1)/(r*r) - E;
  }
  return r*r - E;
}

//Leonard-Jones potential.
double VLJ(double r) {
  return epsilon*( TMath::Power((rho/r),12) - 2*TMath::Power((rho/r),6) );
}

//Leonard-Jones effective potential.
double FLJ(int l, double r, double E) {
  return VLJ(r) + l*(l+1)/(r*r) - E;
  //(void)l;
  //return VLJ(r) - E;
}
	     
// h: step size
// fvals: input vector of energy funcion (V_eff - E) along trajectory
// psi1: starting value of wavefcn u(x0)
// psi2: 2nd starting value u(x0+h)
// traj: return values of trajectory u(x)
void Numerov(double h, const vector<double> &fvals,
             double psi1, double psi2, vector<double> &traj){
  traj = fvals;
  double hsq = h * h;
  double fac = hsq / 12.0;

  int ih = 1;   // step direction
  int startI = 0;
  int endI   = traj.size() - 1;
  if(h < 0) {
    ih = -1;
    startI = traj.size() - 1;
    endI   = 0;
  }
  // to do: handle singular condition
  double wPrev = (1.0 - fac * fvals[startI]) * psi1;
  traj[startI] = psi1;

  double psi = psi2;
  traj[startI + ih] = psi2;
  double w = (1 - fac * fvals[startI + ih]) * psi;

  for(int i = startI + ih; i != endI; i += ih) {
    double wNext = w * 2.0 - wPrev + hsq * psi * fvals[i];
    wPrev = w;
    w     = wNext;
    psi   = w / (1 - fac * fvals[i + ih]);
    traj[i + ih] = psi;
  }
}

//Approximation of u(r), the radial wavefunction, for small r
double u_small(double r) {
  double C = TMath::Sqrt(epsilon*alpha/25);
  return TMath::Exp(-1*C*TMath::Power(r, -5));
}

int main(int argc, char ** argv) {
  // ------------------Setting up ROOT Display----------------------------

  TApplication theApp("App", &argc, argv); // init ROOT App for display
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height 
  UInt_t dw = 1.1*dh;
  TCanvas *c1 = new TCanvas("c1","Solutions",dw,dh);
  c1->cd();
  //---------------------------------------------------------------------


  double h = 0.01;
  double r0 = 1;
  //Our initial guesses u0 and u1 are based on the small-value approximation.
  double u0 = u_small(r0);
  double u1 = u_small(r0+h);
  //range of r values to be searched
  double rmin = r0;
  double rmax = 5.0;
  double r;

  //Number of r values that will be calculated
  int nsteps = (int) (((rmax-rmin)/h)+0.5);
  
  //Range of energies.
  double Emin = 0.1;
  double Emax = 3.5;
  double E;
  //Number of different energy values that will be calculated.
  int Esteps = 1000;
  
  //Range of l, the quantum angular momentum number.
  int lmin = 0;
  int lmax = 6;
  
  //These arrays hold calculations for later use.
  double results[Esteps][lmax-lmin+1];
  double E_vals[Esteps];

  //The initial values of the effective potential, which are used for Numerov
  vector<double> fVals(nsteps);
  //u(r), the radial wavefunction, will have values stored here
  vector<double> traj(nsteps);
  //here we store what r value is at which index
  vector<double> r_values(nsteps);
  
  
  //Here begins the loop for calculating values of sigma_l for various energies.
  for (int e = 0; e < Esteps; e++) {
    //Looping over different values of energy
    E = Emin + ((Emax-Emin)/Esteps)*e;
    E_vals[e] = E;
    for (int l = lmin; l <= lmax; l++) {
      //For a certain value of energy, looping over
      //different values of l
      for (int i = 0; i <= nsteps; i++) {
	//Putting in the values of the effective potential into fVals
	r = rmin + h*i;
	fVals[i] = FLJ(l,r,E);
	if ((e==0) and (l == lmin)) { //Lets us remember the r values
	  r_values[i] = r;
	}
      } 
      Numerov(h, fVals, u0, u1, traj); //Solving with Numerov

      //-------The following code calculates sigma_l---------
      double k = TMath::Sqrt(alpha) * TMath::Sqrt(E);
      double lam = 2*TMath::Pi() / k;
      
      int index2 = nsteps-1;
      //int index_difference = (int)(k*h);
      int index_difference = 1;
      //cout << "index diff " << index_difference << endl;
      int index1 = index2 - index_difference;

      double u_1 = traj[index1];
      double u_2 = traj[index2];
      double r_1 = r_values[index1];
      double r_2 = r_values[index2];
      //cout << "u_1 = " << u_1 <<", u_2 = "<<u_2<<endl;
      //cout << "r_1 = " << r_1 <<", r_2 = "<<r_2<<endl;
      double K = (r_1*u_2) / (r_2*u_1);
      //cout << "K = " << K << endl;
      double j_1 = gsl_sf_bessel_jl(l,k*r_1);
      double j_2 = gsl_sf_bessel_jl(l,k*r_2);
      double n_1 = gsl_sf_bessel_yl(l,k*r_1);
      double n_2 = gsl_sf_bessel_yl(l,k*r_2);

      double tan_sigma = (K*j_1 - j_2) / (K*n_1 - n_2);
      //cout << "tan_sigma = " << tan_sigma << endl; 
      double sigma = TMath::ATan(tan_sigma);
      
      //Storing our result
      results[e][l] = sigma;
      
    }
  }
  

  //The cross section for a given energy is calculated by summing up
  //sigma_l for all values of l at that Energy. Here we do this to obtain
  //our final results for the cross section.
  double cross_sections[Esteps];
  for (int e = 0; e < Esteps; e++) {
    double sig_tot = 0;
    for (int l = lmin; l <= lmax; l++) {
      sig_tot += (2*l + 1)*TMath::Sin(results[e][l])*TMath::Sin(results[e][l]);
    }
    double k = TMath::Sqrt(alpha) * TMath::Sqrt(E_vals[e]);
    sig_tot *= 4*TMath::Pi() / (k*k);
    cross_sections[e] = sig_tot;
  }
  TGraph *cs = new TGraph(Esteps, E_vals, cross_sections);
  
  /*
  double u[nsteps];
  double x[nsteps];
  for (int i = 0; i < nsteps; i++) {
    u[i] = traj[i];
    x[i] = (double)i;
  }
  TGraph *cs = new TGraph(nsteps,x,u);
  */
  cs->Draw();


  /*
  TGraph *tgTraj = new TGraph();
  for (int i = 0; i <= nsteps; i++) {
    r = rmin + i*h;
    tgTraj->SetPoint(i,r,traj[i]);
  }
  */
  
  //---------Drawing display and running the application-----------------
  
  //tgTraj->Draw();
  c1->Draw();
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(300,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
  //-------------------------------------------------------------------
}
