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
const double ro = 3.57; //ro = 3.57 angstrom
const double alpha = 6.12;

double FHO(int l, double r, double E) {
  if (l > 0) {
    return r*r + l*(l+1)/(r*r) - E;
  }
  return r*r - E;
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

int main(int argc, char ** argv) {
  // ------------------Setting up ROOT Display----------------------------

  TApplication theApp("App", &argc, argv); // init ROOT App for display
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height 
  UInt_t dw = 1.1*dh;
  TCanvas *c1 = new TCanvas("c1","Solutions",dw,dh);
  c1->cd();
  //---------------------------------------------------------------------




  double h = 0.001;
  double u0 = 0;
  double u1 = h;
  int l = 0;
  double E = 3.0;
  
  double rmin = 0;
  double rmax = 5.0;
  double r;

  int nsteps = (int) (((rmax-rmin)/h)+0.5);

  vector<double> fVals;
  vector<double> traj;

  for (int i = 0; i <= nsteps; i++) {
    r = rmin + h*i;
    fVals.push_back( FHO(l,r,E) );
  }

  Numerov(h, fVals, u0, u1, traj);

  TGraph *tgTraj = new TGraph();
  for (int i = 0; i <= nsteps; i++) {
    r = rmin + i*h;
    tgTraj->SetPoint(i,r,traj[i]);
  }

  //---------Drawing display and running the application-----------------
  
  tgTraj->Draw();
  c1->Draw();
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(300,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
  //-------------------------------------------------------------------
}
