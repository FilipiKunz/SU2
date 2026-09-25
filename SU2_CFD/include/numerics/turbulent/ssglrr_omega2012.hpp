/*
 * SSG/LRR-RSM-omega 2012 local model terms.
 * Equations and coefficients: NASA TMR, rsm-ssglrr.html (2012 model).
 * R is the positive Favre covariance per unit density, not the negative stress.
 * gradU[i][j] is dU_i/dx_j. All source terms below are per unit mass.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SSGLRROmega2012 {

using Vector = std::array<double, 3>;
using Tensor = std::array<std::array<double, 3>, 3>;

struct State {
  Tensor R;
  double omega;
  double rho;
  double mu;
  double wallDistance;
  Vector gradK;
  Vector gradOmega;
  Tensor gradU;
};

struct Terms {
  double k;
  double F1;
  double epsilon;
  double omegaSource;
  double omegaProduction;
  double omegaDestruction;
  double omegaCrossDiffusion;
  double omegaDiffusivity;
  double stressDiffusivity[3][3];
  Tensor production;
  Tensor pressureStrain;
  Tensor stressSource;
};

inline double blend(double f, double inner, double outer) {
  return f*inner + (1.0-f)*outer;
}

inline Terms evaluate(const State& s) {
  // Positivity is a transport/BC invariant. Do not silently clip it here.
  Terms t{};
  t.k = 0.5*(s.R[0][0] + s.R[1][1] + s.R[2][2]);
  if (!(t.k > 0.0 && s.omega > 0.0 && s.rho > 0.0 && s.wallDistance > 0.0)) {
    throw std::domain_error("SSG/LRR requires positive k, omega, rho and wall distance");
  }

  double gradDot = 0.0;
  for (int i=0; i<3; ++i) gradDot += s.gradK[i]*s.gradOmega[i];
  gradDot = std::max(gradDot, 0.0);
  const double CD = 1.712*s.rho*gradDot/s.omega;
  const double d = s.wallDistance;
  const double zeta1 = std::max(std::sqrt(t.k)/(0.09*s.omega*d),
                               500.0*s.mu/(s.rho*s.omega*d*d));
  const double zeta2 = CD > 0.0 ? 4.0*0.856*s.rho*t.k/(CD*d*d)
                                : std::numeric_limits<double>::max();
  const double zeta = std::min(zeta1, zeta2);
  const double zeta2pow = zeta*zeta;
  t.F1 = std::tanh(zeta2pow*zeta2pow);

  const double C1 = blend(t.F1, 1.8, 1.7);
  const double C1star = blend(t.F1, 0.0, 0.9);
  const double C2 = blend(t.F1, 0.0, 1.05);
  const double C3 = 0.8;
  const double C3star = blend(t.F1, 0.0, 0.65);
  const double C4 = blend(t.F1, 0.5*(18.0*0.52+12.0)/11.0, 0.625);
  const double C5 = blend(t.F1, 0.5*(-14.0*0.52+20.0)/11.0, 0.2);
  const double D = blend(t.F1, 0.75*0.09, 0.22);
  const double alphaOmega = blend(t.F1, 0.5556, 0.44);
  const double betaOmega = blend(t.F1, 0.075, 0.0828);
  const double sigmaOmega = blend(t.F1, 0.5, 0.856);
  const double sigmaD = blend(t.F1, 0.0, 1.712);

  t.epsilon = 0.09*t.k*s.omega;
  t.omegaDiffusivity = s.mu + sigmaOmega*s.rho*t.k/s.omega;
  Tensor a{}, S{}, W{}, a2{};
  const double divU = s.gradU[0][0]+s.gradU[1][1]+s.gradU[2][2];
  double aa = 0.0, aS = 0.0, Pkk = 0.0;
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) {
    const double delta = i==j ? 1.0 : 0.0;
    a[i][j] = s.R[i][j]/t.k - (2.0/3.0)*delta;
    S[i][j] = 0.5*(s.gradU[i][j]+s.gradU[j][i]);
    W[i][j] = 0.5*(s.gradU[i][j]-s.gradU[j][i]);
    t.stressDiffusivity[i][j] = s.mu*delta + D*s.rho*s.R[i][j]/(0.09*s.omega);
    aa += a[i][j]*a[i][j];
    aS += a[i][j]*S[i][j];
    for (int k=0; k<3; ++k)
      t.production[i][j] -= s.R[i][k]*s.gradU[j][k]
                           + s.R[j][k]*s.gradU[i][k];
  }
  for (int i=0; i<3; ++i) Pkk += t.production[i][i];
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
    for (int k=0; k<3; ++k) a2[i][j] += a[i][k]*a[k][j];

  t.omegaProduction = alphaOmega*s.omega*Pkk/(2.0*t.k);
  t.omegaDestruction = -betaOmega*s.omega*s.omega;
  t.omegaCrossDiffusion = sigmaD*gradDot/s.omega;
  t.omegaSource = t.omegaProduction + t.omegaDestruction + t.omegaCrossDiffusion;

  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) {
    const double delta = i==j ? 1.0 : 0.0;
    double aSsum = 0.0, aWsum = 0.0;
    for (int k=0; k<3; ++k) {
      aSsum += a[i][k]*S[j][k] + a[j][k]*S[i][k];
      aWsum += a[i][k]*W[j][k] + a[j][k]*W[i][k];
    }
    t.pressureStrain[i][j] = -(C1*t.epsilon + 0.5*C1star*Pkk)*a[i][j]
      + C2*t.epsilon*(a2[i][j]-(aa/3.0)*delta)
      + (C3-C3star*std::sqrt(aa))*t.k*(S[i][j]-(divU/3.0)*delta)
      + C4*t.k*(aSsum-(2.0/3.0)*aS*delta)
      + C5*t.k*aWsum;
    t.stressSource[i][j] = t.production[i][j] + t.pressureStrain[i][j]
                            - (2.0/3.0)*t.epsilon*delta;
  }
  return t;
}

} // namespace SSGLRROmega2012
