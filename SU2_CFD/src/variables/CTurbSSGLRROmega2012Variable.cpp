#include "../../include/variables/CTurbSSGLRROmega2012Variable.hpp"

CTurbSSGLRROmega2012Variable::CTurbSSGLRROmega2012Variable(
    su2double kine, su2double omega, su2double mut, unsigned long npoint,
    unsigned long ndim, unsigned long nvar, CConfig *config)
  : CTurbVariable(npoint, ndim, nvar, config) {
  if (nvar != 7) SU2_MPI::Error("SSG/LRR requires six stresses and omega", CURRENT_FUNCTION);
  for (unsigned long iPoint=0; iPoint<nPoint; ++iPoint) {
    for (unsigned short iVar=0; iVar<6; ++iVar) Solution(iPoint,iVar) = 0.0;
    Solution(iPoint,0) = Solution(iPoint,1) = Solution(iPoint,2) = (2.0/3.0)*kine;
    Solution(iPoint,6) = omega;
  }
  Solution_Old = Solution;
  F1.resize(nPoint) = su2double(1.0);
  F2.resize(nPoint) = su2double(0.0);
  CDkw.resize(nPoint) = su2double(0.0);
  muT.resize(nPoint) = mut;
}

void CTurbSSGLRROmega2012Variable::SetBlendingFunc(
    unsigned long iPoint, su2double mu, su2double d, su2double rho,
    TURB_TRANS_MODEL) {
  // NASA TMR 2012 F1; grad(k) = 1/2 trace grad(R).
  const su2double k = 0.5*(Solution(iPoint,0)+Solution(iPoint,1)+Solution(iPoint,2));
  const su2double omega = Solution(iPoint,6);
  if (k <= 0.0 || omega <= 0.0 || rho <= 0.0 || d <= 0.0)
    SU2_MPI::Error("Nonpositive RSM state in F1", CURRENT_FUNCTION);
  su2double gradDot = 0.0;
  for (unsigned short j=0; j<nDim; ++j)
    gradDot += 0.5*(Gradient(iPoint,0,j)+Gradient(iPoint,1,j)+Gradient(iPoint,2,j))
                     *Gradient(iPoint,6,j);
  gradDot = max(gradDot, su2double(0.0));
  CDkw(iPoint) = 1.712*rho*gradDot/omega;
  const su2double zeta1 = max(sqrt(k)/(0.09*omega*d), 500.0*mu/(rho*omega*d*d));
  const su2double zeta2 = CDkw(iPoint) > 0.0 ? 4.0*0.856*rho*k/(CDkw(iPoint)*d*d) : 1e100;
  const su2double zeta = min(zeta1,zeta2);
  F1(iPoint) = tanh(pow(zeta,4.0));
}
