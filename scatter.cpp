%%cpp
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
