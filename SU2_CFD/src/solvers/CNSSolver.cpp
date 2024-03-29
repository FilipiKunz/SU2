/*!
 * \file CNSSolver.cpp
 * \brief Main subroutines for solving Finite-Volume Navier-Stokes flow problems.
 * \author F. Palacios, T. Economon
 * \version 7.5.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/solvers/CNSSolver.hpp"
#include "../../include/variables/CNSVariable.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CFVMFlowSolverBase.inl"

/*--- Explicit instantiation of the parent class of CEulerSolver,
 *    to spread the compilation over two cpp files. ---*/
template class CFVMFlowSolverBase<CEulerVariable, ENUM_REGIME::COMPRESSIBLE>;


CNSSolver::CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) :
           CEulerSolver(geometry, config, iMesh, true) {

  /*--- This constructor only allocates/inits what is extra to CEulerSolver. ---*/

  unsigned short iMarker, iDim;
  unsigned long iVertex;


  /*--- Allocates a 2D array with variable "outer" sizes and init to 0. ---*/

  auto Alloc2D = [](unsigned long M, const std::vector<unsigned long>& N, su2double**& X) {
      X = new su2double*[M];

      for (unsigned long i = 0; i < M; ++i) {
          X[i] = new su2double[N[i]]();
      }
  };

  /*--- Allocates a 3D array with variable "middle" sizes and init to 0. ---*/

  auto Alloc3D = [](unsigned long M, const std::vector<unsigned long>& N, unsigned short P, su2double***& X)->void {
      X = new su2double**[M];

      for (unsigned long i = 0; i < M; ++i) {
          X[i] = new su2double*[N[i]];

          for (unsigned long j = 0; j < N[i]; ++j) {
              X[i][j] = new su2double[P]();
          }
      }
  };


  /*--- Buffet sensor in all the markers and coefficients ---*/

  Buffet_Sensor.resize(nMarker);
  for (unsigned long i = 0; i< nMarker; ++i) Buffet_Sensor[i].resize(nVertex[i], 0.0);
  Buffet_Metric.resize(nMarker, 0.0);
  Surface_Buffet_Metric.resize(config->GetnMarker_Monitoring(), 0.0);
  /*--- Read farfield conditions from config ---*/

  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();
  Tke_Inf         = config->GetTke_FreeStreamND();


  /*--- Set the SGS model in case an LES simulation is carried out.
   Make a distinction between the SGS models used and set SGSModel and
  SGSModelUsed accordingly. ---*/

  SGSModel = nullptr;
  SGSModelUsed = false;

  switch( config->GetKind_SGS_Model() ) {
    /* No LES, so no SGS model needed.
     Set the pointer to NULL and the boolean to false. */
    case TURB_SGS_MODEL::NONE: case TURB_SGS_MODEL::IMPLICIT_LES:
      SGSModel = nullptr;
      SGSModelUsed = false;
      break;
    case TURB_SGS_MODEL::SMAGORINSKY:
      SGSModel     = new CSmagorinskyModel;
      SGSModelUsed = true;
      break;
    case TURB_SGS_MODEL::WALE:
      SGSModel     = new CWALEModel;
      SGSModelUsed = true;
      break;
    case TURB_SGS_MODEL::VREMAN:
      SGSModel     = new CVremanModel;
      SGSModelUsed = true;
      break;
    default:
      SU2_MPI::Error("Unknown SGS model encountered", CURRENT_FUNCTION);
  }

  /*--- Set the wall model to NULL ---*/

  WallModel = NULL;

  if (config->GetWall_Models()){

    /*--- Set the WMLES class  ---*/
    /*--- First allocate the auxiliary variables ---*/
    Alloc2D(nMarker, nVertex, TauWall_WMLES);
    Alloc2D(nMarker, nVertex, HeatFlux_WMLES);
    Alloc3D(nMarker, nVertex, nDim, FlowDirTan_WMLES);
    Alloc3D(nMarker, nVertex, nDim, VelTimeFilter_WMLES);

    /*--- Check if the Wall models or Wall functions are unique. ---*/
    /*--- OBS: All the markers must have the same wall model/function ---*/

    vector<unsigned short> WallFunctions_;
    vector<string> WallFunctionsMarker_;
    for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
      switch (config->GetMarker_All_KindBC(iMarker)) {
        case ISOTHERMAL:
        case HEAT_FLUX: {
          string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if(config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NO_WALL_FUNCTION){
            WallFunctions_.push_back(static_cast<unsigned short>(config->GetWallFunction_Treatment(Marker_Tag)));
            WallFunctionsMarker_.push_back(Marker_Tag);
          }
          break;
        }
        default:  /* Just to avoid a compiler warning. */
          break;
      }
    }

    if (!WallFunctions_.empty()){
      sort(WallFunctions_.begin(), WallFunctions_.end());
      vector<unsigned short>::iterator it = std::unique( WallFunctions_.begin(), WallFunctions_.end() );
      WallFunctions_.erase(it, WallFunctions_.end());

      if(WallFunctions_.size() == 1) {
        switch (config->GetWallFunction_Treatment(WallFunctionsMarker_[0])) {
          case WALL_FUNCTIONS::EQUILIBRIUM_WALL_MODEL:
            WallModel = new CWallModel1DEQ(config,WallFunctionsMarker_[0]);
            break;
          case WALL_FUNCTIONS::LOGARITHMIC_WALL_MODEL:
            WallModel = new CWallModelLogLaw(config,WallFunctionsMarker_[0]);
            break;
          case WALL_FUNCTIONS::ALGEBRAIC_WALL_MODEL:
            WallModel = new CWallModelAlgebraic(config,WallFunctionsMarker_[0]);
            break;
          case WALL_FUNCTIONS::APGLL_WALL_MODEL:
            WallModel = new CWallModelAPGLL(config,WallFunctionsMarker_[0]);
            break;
          case WALL_FUNCTIONS::TEMPLATE_WALL_MODEL:
            WallModel = new CWallModelTemplate(config,WallFunctionsMarker_[0]);
          break;

          default:
            break;
        }
      }
      else{
        SU2_MPI::Error("Wall function/model type must be the same in all wall BCs", CURRENT_FUNCTION);
      }
    }
  }

  /*--- Initialize the seed values for forward mode differentiation. ---*/

  switch(config->GetDirectDiff()) {
    case D_VISCOSITY:
      SU2_TYPE::SetDerivative(Viscosity_Inf, 1.0);
      break;
    default:
      /*--- Already done upstream. ---*/
      break;
  }

}

void CNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                              unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  const auto InnerIter = config->GetInnerIter();
  const bool muscl = config->GetMUSCL_Flow() && (iMesh == MESH_0);
  const bool center = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  const bool limiter = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (InnerIter <= config->GetLimiterIter());
  const bool van_albada = (config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE);
  const bool wall_functions = config->GetWall_Functions();
  bool wall_models          = config->GetWall_Models();

  /*--- Common preprocessing steps (implemented by CEulerSolver) ---*/

  CommonPreprocessing(geometry, solver_container, config, iMesh, iRKStep, RunTime_EqSystem, Output);

  /*--- Compute gradient for MUSCL reconstruction, for output (i.e. the
   turbulence solver, and post) only temperature and velocity are needed ---*/

  const auto nPrimVarGrad_bak = nPrimVarGrad;
  if (Output) ompMasterAssignBarrier(nPrimVarGrad, 1+nDim);

  if (config->GetReconstructionGradientRequired() && muscl && !center) {
    switch (config->GetKind_Gradient_Method_Recon()) {
      case GREEN_GAUSS:
        SetPrimitive_Gradient_GG(geometry, config, true); break;
      case LEAST_SQUARES:
      case WEIGHTED_LEAST_SQUARES:
        SetPrimitive_Gradient_LS(geometry, config, true); break;
      default: break;
    }
  }

  /*--- Compute gradient of the primitive variables ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
  }
  else if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
  }

  if (Output) ompMasterAssignBarrier(nPrimVarGrad, nPrimVarGrad_bak);

  /*--- Compute the limiters ---*/

  if (muscl && !center && limiter && !van_albada && !Output) {
    SetPrimitive_Limiter(geometry, config);
  }

  ComputeVorticityAndStrainMag(*config, geometry, iMesh);

  /*--- Compute the TauWall from the wall functions ---*/
  /*--- Calculate the eddy viscosity using a SGS model ---*/

  if (SGSModelUsed){
    SU2_OMP_MASTER
    Setmut_LES(geometry, solver_container, config);
    SU2_OMP_BARRIER
  }

  /*--- Compute the wall shear stress from the wall functions ---*/

  if (wall_functions) {
    SU2_OMP_MASTER
    SetTau_Wall_WF(geometry, solver_container, config);
    SetEddyViscFirstPoint(geometry, solver_container, config);
    SU2_OMP_BARRIER
  }

  /*--- Compute the wall shear stress from the wall model ---*/

  if (wall_models && (iRKStep==0) && (iMesh == MESH_0)){
    SU2_OMP_MASTER
    SetTauWallHeatFlux_WMLES1stPoint(geometry, solver_container, config, iRKStep);
    SU2_OMP_BARRIER
  }

}

unsigned long CNSSolver::SetPrimitive_Variables(CSolver **solver_container, const CConfig *config) {

  /*--- Number of non-physical points, local to the thread, needs
   *    further reduction if function is called in parallel ---*/
  unsigned long nonPhysicalPoints = 0;

  const TURB_MODEL turb_model = config->GetKind_Turb_Model();
  const bool tkeNeeded = (turb_model == TURB_MODEL::SST);

  AD::StartNoSharedReading();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Retrieve the value of the kinetic energy (if needed). ---*/

    su2double eddy_visc = 0.0, turb_ke = 0.0;

    if (turb_model != TURB_MODEL::NONE && solver_container[TURB_SOL] != nullptr) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);

      if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES) {
        su2double DES_LengthScale = solver_container[TURB_SOL]->GetNodes()->GetDES_LengthScale(iPoint);
        nodes->SetDES_LengthScale(iPoint, DES_LengthScale);
      }
    }

    /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/

    bool physical = static_cast<CNSVariable*>(nodes)->SetPrimVar(iPoint, eddy_visc, turb_ke, GetFluidModel());
    nodes->SetSecondaryVar(iPoint, GetFluidModel());

    /*--- Check for non-realizable states for reporting. ---*/

    nonPhysicalPoints += !physical;

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

  return nonPhysicalPoints;
}

void CNSSolver::Viscous_Residual(unsigned long iEdge, CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config) {

  Viscous_Residual_impl(iEdge, geometry, solver_container, numerics, config);
}

void CNSSolver::Buffet_Monitoring(const CGeometry *geometry, const CConfig *config) {

  unsigned long iVertex, iMarker;
  unsigned short iMarker_Monitoring;
  const su2double* Vel_FS = Velocity_Inf;
  const su2double k = config->GetBuffet_k(), lam = config->GetBuffet_lambda(), Sref = config->GetRefArea();

  const su2double VelMag_FS = GeometryToolbox::Norm(nDim, Vel_FS);

  /*-- Variables initialization ---*/

  Total_Buffet_Metric = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_Buffet_Metric[iMarker_Monitoring] = 0.0;
  }

  /*--- Loop over the Euler and Navier-Stokes markers ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Buffet_Metric[iMarker] = 0.0;

    const auto Monitoring = config->GetMarker_All_Monitoring(iMarker);

    if (config->GetViscous_Wall(iMarker)) {

      /*--- Loop over the vertices to compute the buffet sensor ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        /*--- Perform dot product of skin friction with freestream velocity ---*/

        const su2double SkinFrictionMag = GeometryToolbox::Norm(nDim, CSkinFriction[iMarker][iVertex]);
        su2double SkinFrictionDot = GeometryToolbox::DotProduct(nDim, CSkinFriction[iMarker][iVertex], Vel_FS);

        /*--- Normalize the dot product ---*/

        SkinFrictionDot /= SkinFrictionMag*VelMag_FS;

        /*--- Compute Heaviside function ---*/

        Buffet_Sensor[iMarker][iVertex] = 1./(1. + exp(2.*k*(SkinFrictionDot + lam)));

        /*--- Integrate buffet sensor ---*/

        if (Monitoring == YES){

          auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          su2double Area = GeometryToolbox::Norm(nDim, Normal);

          Buffet_Metric[iMarker] += Buffet_Sensor[iMarker][iVertex]*Area/Sref;

        }

      }

      if (Monitoring == YES){

        Total_Buffet_Metric += Buffet_Metric[iMarker];

        /*--- Per surface buffet metric ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          auto Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag)
            Surface_Buffet_Metric[iMarker_Monitoring] = Buffet_Metric[iMarker];
        }

      }

    }

  }

  /*--- Add buffet metric information using all the nodes ---*/

  su2double MyTotal_Buffet_Metric = Total_Buffet_Metric;
  SU2_MPI::Allreduce(&MyTotal_Buffet_Metric, &Total_Buffet_Metric, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  /*--- Add the buffet metric on the surfaces using all the nodes ---*/

  auto local_copy = Surface_Buffet_Metric;
  SU2_MPI::Allreduce(local_copy.data(), Surface_Buffet_Metric.data(), local_copy.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

}

void CNSSolver::Evaluate_ObjFunc(const CConfig *config, CSolver**) {

  unsigned short iMarker_Monitoring, Kind_ObjFunc;
  su2double Weight_ObjFunc;

  /*--- Evaluate objective functions common to Euler and NS solvers ---*/

  CEulerSolver::Evaluate_ObjFunc(config, nullptr);

  /*--- Evaluate objective functions specific to NS solver ---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {

    Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);
    Kind_ObjFunc = config->GetKind_ObjFunc(iMarker_Monitoring);

    switch(Kind_ObjFunc) {
      case BUFFET_SENSOR:
          Total_ComboObj +=Weight_ObjFunc*Surface_Buffet_Metric[iMarker_Monitoring];
          break;
      default:
          break;
    }
  }

}

void CNSSolver::SetRoe_Dissipation(CGeometry *geometry, CConfig *config){

  const unsigned short kind_roe_dissipation = config->GetKind_RoeLowDiss();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    if (kind_roe_dissipation == FD || kind_roe_dissipation == FD_DUCROS){

      su2double wall_distance = geometry->nodes->GetWall_Distance(iPoint);

      nodes->SetRoe_Dissipation_FD(iPoint, wall_distance);

    } else if (kind_roe_dissipation == NTS || kind_roe_dissipation == NTS_DUCROS) {

      const su2double delta = geometry->nodes->GetMaxLength(iPoint);
      assert(delta > 0 && "Delta must be initialized and non-negative");
      nodes->SetRoe_Dissipation_NTS(iPoint, delta, config->GetConst_DES());
    }
  }
  END_SU2_OMP_FOR

}

void CNSSolver::AddDynamicGridResidualContribution(unsigned long iPoint, unsigned long Point_Normal,
                                                   const CGeometry* geometry,  const su2double* UnitNormal,
                                                   su2double Area, const su2double* GridVel,
                                                   su2double** Jacobian_i, su2double& Res_Conv,
                                                   su2double& Res_Visc) const {

  su2double ProjGridVel = Area * GeometryToolbox::DotProduct(nDim, GridVel, UnitNormal);

  /*--- Retrieve other primitive quantities and viscosities ---*/

  su2double Density = nodes->GetDensity(iPoint);
  su2double Pressure = nodes->GetPressure(iPoint);
  su2double laminar_viscosity = nodes->GetLaminarViscosity(iPoint);
  su2double eddy_viscosity = nodes->GetEddyViscosity(iPoint);
  su2double total_viscosity = laminar_viscosity + eddy_viscosity;

  /*--- Compute the viscous stress tensor ---*/

  su2double tau[MAXNDIM][MAXNDIM] = {{0.0}};
  CNumerics::ComputeStressTensor(nDim, tau, nodes->GetVelocityGradient(iPoint), total_viscosity);

  /*--- Dot product of the stress tensor with the grid velocity ---*/

  su2double tau_vel[MAXNDIM] = {0.0};
  for (auto iDim = 0u; iDim < nDim; iDim++)
    tau_vel[iDim] = GeometryToolbox::DotProduct(nDim, tau[iDim], GridVel);

  /*--- Compute the convective and viscous residuals (energy eqn.) ---*/

  Res_Conv += Pressure*ProjGridVel;
  Res_Visc += GeometryToolbox::DotProduct(nDim, tau_vel, UnitNormal) * Area;

  /*--- Implicit Jacobian contributions due to moving walls ---*/

  if (Jacobian_i != nullptr) {

    /*--- Jacobian contribution related to the pressure term ---*/

    su2double GridVel2 = GeometryToolbox::SquaredNorm(nDim, GridVel);

    Jacobian_i[nDim+1][0] += 0.5*(Gamma-1.0)*GridVel2*ProjGridVel;

    for (auto jDim = 0u; jDim < nDim; jDim++)
      Jacobian_i[nDim+1][jDim+1] += -(Gamma-1.0)*GridVel[jDim]*ProjGridVel;

    Jacobian_i[nDim+1][nDim+1] += (Gamma-1.0)*ProjGridVel;

    /*--- Now the Jacobian contribution related to the shear stress ---*/

    /*--- Get coordinates of i & nearest normal and compute distance ---*/

    const auto Coord_i = geometry->nodes->GetCoord(iPoint);
    const auto Coord_j = geometry->nodes->GetCoord(Point_Normal);

    su2double dist_ij = GeometryToolbox::Distance(nDim, Coord_i, Coord_j);

    const su2double theta2 = 1.0;

    su2double factor = total_viscosity*Area/(Density*dist_ij);

    if (nDim == 2) {
      su2double thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
      su2double thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;

      su2double etaz = UnitNormal[0]*UnitNormal[1]/3.0;

      su2double pix = GridVel[0]*thetax + GridVel[1]*etaz;
      su2double piy = GridVel[0]*etaz   + GridVel[1]*thetay;

      Jacobian_i[nDim+1][0] += factor*(-pix*GridVel[0]+piy*GridVel[1]);
      Jacobian_i[nDim+1][1] += factor*pix;
      Jacobian_i[nDim+1][2] += factor*piy;
    }
    else {
      su2double thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
      su2double thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
      su2double thetaz = theta2 + UnitNormal[2]*UnitNormal[2]/3.0;

      su2double etaz = UnitNormal[0]*UnitNormal[1]/3.0;
      su2double etax = UnitNormal[1]*UnitNormal[2]/3.0;
      su2double etay = UnitNormal[0]*UnitNormal[2]/3.0;

      su2double pix = GridVel[0]*thetax + GridVel[1]*etaz   + GridVel[2]*etay;
      su2double piy = GridVel[0]*etaz   + GridVel[1]*thetay + GridVel[2]*etax;
      su2double piz = GridVel[0]*etay   + GridVel[1]*etax   + GridVel[2]*thetaz;

      Jacobian_i[nDim+1][0] += factor*(-pix*GridVel[0]+piy*GridVel[1]+piz*GridVel[2]);
      Jacobian_i[nDim+1][1] += factor*pix;
      Jacobian_i[nDim+1][2] += factor*piy;
      Jacobian_i[nDim+1][3] += factor*piz;
    }
  }
}

void CNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver**, CNumerics*,
                                 CNumerics*, CConfig *config, unsigned short val_marker) {

  BC_HeatFlux_Wall_Generic(geometry, config, val_marker, HEAT_FLUX);
}

void CNSSolver::BC_HeatTransfer_Wall(const CGeometry *geometry, const CConfig *config, const unsigned short val_marker) {

  BC_HeatFlux_Wall_Generic(geometry, config, val_marker, HEAT_TRANSFER);
}

void CNSSolver::BC_HeatFlux_Wall_Generic(const CGeometry* geometry, const CConfig* config, unsigned short val_marker,
                                         unsigned short kind_boundary) {
  /*--- Identify the boundary by string name and get the specified wall
   heat flux from config as well as the wall function treatment. ---*/

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux, temperature or heat transfer coefficient from config ---*/

  su2double Wall_HeatFlux = 0.0, Tinfinity = 0.0, Transfer_Coefficient = 0.0;

  if (kind_boundary == HEAT_FLUX) {
    Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag) / config->GetHeat_Flux_Ref();
    if (config->GetIntegrated_HeatFlux()) {
      Wall_HeatFlux /= geometry->GetSurfaceArea(config, val_marker);
    }
  } else if (kind_boundary == HEAT_TRANSFER) {
    /*--- The required heatflux will be computed for each iPoint individually based on local Temperature. ---*/
    Transfer_Coefficient = config->GetWall_HeatTransfer_Coefficient(Marker_Tag) * config->GetTemperature_Ref() /
                           config->GetHeat_Flux_Ref();
    Tinfinity = config->GetWall_HeatTransfer_Temperature(Marker_Tag) / config->GetTemperature_Ref();
  }

//  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
//  if (Wall_Function != WALL_FUNCTION::NONE) {
//    SU2_MPI::Error("Wall function treatment not implemented yet", CURRENT_FUNCTION);
//  }

  /*--- Jacobian, initialized to zero if needed. ---*/
  su2double **Jacobian_i = nullptr;
  if ((dynamic_grid || (kind_boundary == HEAT_TRANSFER)) && implicit) {
    Jacobian_i = new su2double* [nVar];
    for (auto iVar = 0u; iVar < nVar; iVar++)
      Jacobian_i[iVar] = new su2double [nVar] ();
  }

  /*--- Loop over all of the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- If it is a customizable patch, retrieve the specified wall heat flux. ---*/

    if (config->GetMarker_All_PyCustom(val_marker))
      Wall_HeatFlux = geometry->GetCustomBoundaryHeatFlux(val_marker, iVertex);
    else if (kind_boundary == HEAT_TRANSFER) {
      const su2double Twall = nodes->GetTemperature(iPoint);
      Wall_HeatFlux = Transfer_Coefficient * (Tinfinity - Twall);
    }

    /*--- Compute dual-grid area and boundary normal ---*/

    const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

    su2double Area = GeometryToolbox::Norm(nDim, Normal);

    su2double UnitNormal[MAXNDIM] = {0.0};
    for (auto iDim = 0u; iDim < nDim; iDim++)
      UnitNormal[iDim] = -Normal[iDim]/Area;

    /*--- Apply a weak boundary condition for the energy equation.
     Compute the residual due to the prescribed heat flux.
     The convective part will be zero if the grid is not moving. ---*/

    su2double Res_Conv = 0.0;
    su2double Res_Visc = Wall_HeatFlux * Area;

    /*--- Impose the value of the velocity as a strong boundary
     condition (Dirichlet). Fix the velocity and remove any
     contribution to the residual at this node. ---*/

    if (dynamic_grid) {
      nodes->SetVelocity_Old(iPoint, geometry->nodes->GetGridVel(iPoint));
    }
    else {
      su2double zero[MAXNDIM] = {0.0};
      nodes->SetVelocity_Old(iPoint, zero);
    }

    for (auto iDim = 0u; iDim < nDim; iDim++)
      LinSysRes(iPoint, iDim+1) = 0.0;
    nodes->SetVel_ResTruncError_Zero(iPoint);

    /*--- If the wall is moving, there are additional residual contributions
     due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/

    if (dynamic_grid) {
      if (implicit) {
        for (auto iVar = 0u; iVar < nVar; ++iVar)
          Jacobian_i[nDim+1][iVar] = 0.0;
      }

      const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      AddDynamicGridResidualContribution(iPoint, Point_Normal, geometry, UnitNormal,
                                         Area, geometry->nodes->GetGridVel(iPoint),
                                         Jacobian_i, Res_Conv, Res_Visc);
    }

    /*--- Convective and viscous contributions to the residual at the wall ---*/

    LinSysRes(iPoint, nDim+1) += Res_Conv - Res_Visc;

    /*--- Enforce the no-slip boundary condition in a strong way by
     modifying the velocity-rows of the Jacobian (1 on the diagonal).
     And add the contributions to the Jacobian due to energy. ---*/

    if (implicit) {
      if (kind_boundary == HEAT_TRANSFER){

        /*--- It is necessary to zero the jacobian entries of the energy equation. ---*/
        if (!dynamic_grid)
          for (auto iVar = 0u; iVar < nVar; ++iVar)
            Jacobian_i[nDim+1][iVar] = 0.0;

        const su2double oneOnRho = 1.0 / nodes->GetDensity(iPoint);
        const su2double oneOnCv = (Gamma - 1.0) / config->GetGas_ConstantND();
        const su2double Vel2 = nodes->GetVelocity2(iPoint);
        const su2double dTdrho = oneOnRho * ( -Tinfinity + oneOnCv * 0.5 * Vel2);
        const su2double dTdrhoe = oneOnCv * oneOnRho;

        /*--- Total specific energy: e=c_v*T+1/2*v^2 => T=1/c_v(rho*e/rho - 1/2||rho v||^2/rho^2).
        Together with cv=R/(gamma-1) the following Jacobian contributions for the energy equation can be derived. ---*/
        Jacobian_i[nDim+1][0] += Transfer_Coefficient * dTdrho * Area;

        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          Jacobian_i[nDim+1][iDim+1] -= Transfer_Coefficient * dTdrhoe * nodes->GetVelocity(iPoint, iDim) * Area;

        Jacobian_i[nDim+1][nDim+1] += Transfer_Coefficient * dTdrhoe * Area;

      }
      if (dynamic_grid || (kind_boundary == HEAT_TRANSFER)) {
        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
      }

      for (auto iVar = 1u; iVar <= nDim; iVar++) {
        auto total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    }
  }
  END_SU2_OMP_FOR

  if (Jacobian_i)
    for (auto iVar = 0u; iVar < nVar; iVar++)
      delete [] Jacobian_i[iVar];
  delete [] Jacobian_i;

}

void CNSSolver::BC_WallModel(CGeometry      *geometry,
                                CSolver        **solver_container,
                                CNumerics      *conv_numerics,
                                CNumerics      *visc_numerics,
                                CConfig        *config,
                                unsigned short val_marker) {

  unsigned short iDim, iVar;
  unsigned long iVertex, iPoint, total_index;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool HeatFlux_Prescribed = false;

  /*--- Allocation of variables necessary for convective fluxes. ---*/
  su2double Area, ProjVelocity_i, Wall_HeatFlux;
  su2double *V_reflected, *V_domain;
  su2double *Normal     = new su2double[nDim];
  su2double *UnitNormal = new su2double[nDim];

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  if(config->GetMarker_All_KindBC(val_marker) == HEAT_FLUX) {
    HeatFlux_Prescribed = true;
  }

  /*--- Loop over all the vertices on this boundary marker. ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*-------------------------------------------------------------------------------*/
      /*--- Step 1: For the convective fluxes, create a reflected state of the      ---*/
      /*---         Primitive variables by copying all interior values to the       ---*/
      /*---         reflected. Only the velocity is mirrored for the wall model     ---*/
      /*---         and negative for wall functions (weakly impose v = 0)           ---*/
      /*---         axis. Based on the Upwind_Residual routine.                     ---*/
      /*-------------------------------------------------------------------------------*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++) {
        UnitNormal[iDim] = -Normal[iDim]/Area;
      }

      /*--- Allocate the reflected state at the symmetry boundary. ---*/
      V_reflected = GetCharacPrimVar(val_marker, iVertex);

      /*--- Grid movement ---*/
      if (config->GetGrid_Movement())
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      /*--- Normal vector for this vertex (negate for outward convention). ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Get current solution at this boundary node ---*/
      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Set the reflected state based on the boundary node. ---*/
      for(iVar = 0; iVar < nPrimVar; iVar++)
        V_reflected[iVar] = nodes->GetPrimitive(iPoint,iVar);

      /*--- Compute velocity in normal direction (ProjVelcity_i=(v*n)) and substract from
       velocity in normal direction: v_r = v - (v*n)n ---*/
      ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity_i += nodes->GetVelocity(iPoint,iDim)*UnitNormal[iDim];

      if (nodes->GetTauWall_Flag(iPoint)){
        /*--- Scalars are copied and the velocity is mirrored along the wall boundary,
         i.e. the velocity in normal direction is substracted twice. ---*/

        /*--- Force the velocity to be tangential ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          V_reflected[iDim+1] = nodes->GetVelocity(iPoint,iDim) - 2.0 * ProjVelocity_i*UnitNormal[iDim];

        /*--- Set Primitive and Secondary for numerics class. ---*/

        conv_numerics->SetPrimitive(V_domain, V_reflected);
        conv_numerics->SetSecondary(nodes->GetSecondary(iPoint), nodes->GetSecondary(iPoint));

        /*--- Compute the residual using an upwind scheme. ---*/

        auto residual = conv_numerics->ComputeResidual(config);

        LinSysRes.AddBlock(iPoint, residual);

        /*--- Jacobian contribution for implicit integration. ---*/
        if (implicit)
          Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      }
      else{

        for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar] = 0.0;

        /*--- Store the corrected velocity at the wall which will
         be zero (v = 0), unless there are moving walls (v = u_wall)---*/
        /*--- TODO: Before adding the moving walls capability make sure that it is running correctly. ---*/

        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;

        /*--- Impose the value of the velocity as a strong boundary
         condition (Dirichlet). Fix the velocity and remove any
         contribution to the residual at this node. ---*/


        nodes->SetVelocity_Old(iPoint,Vector);

        for (iDim = 0; iDim < nDim; iDim++)
          LinSysRes.SetBlock_Zero(iPoint, iDim+1);
        nodes->SetVel_ResTruncError_Zero(iPoint);

        /*--- Update residual value ---*/
        LinSysRes.AddBlock(iPoint, Res_Conv);

        /*--- Enforce the no-slip boundary condition in a strong way by
         modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/
        if (implicit){
          for (iVar = 1; iVar <= nDim; iVar++) {
            total_index = iPoint*nVar+iVar;
            Jacobian.DeleteValsRowi(total_index);
          }
        }


      }


      /*-------------------------------------------------------*/
      /*-------------------------------------------------------*/
      /*--- Viscous residual contribution of the wall model ---*/
      /*--- TODO: Build the jacobian contribution of the WM ---*/
      /*-------------------------------------------------------*/
      /*-------------------------------------------------------*/

      if (nodes->GetTauWall_Flag(iPoint)){

        /*--- Weakly enforce the WM heat flux for the energy equation---*/
        su2double velWall_tan = 0.;
        su2double DirTanWM[3] = {0.,0.,0.};

        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          DirTanWM[iDim] = GetFlowDirTan_WMLES(val_marker,iVertex,iDim);

        const su2double TauWall = GetTauWall_WMLES(val_marker,iVertex);
        const su2double Wall_HeatFlux = GetHeatFlux_WMLES(val_marker, iVertex);

        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          velWall_tan +=  nodes->GetVelocity(iPoint,iDim) * DirTanWM[iDim];

        Res_Visc[0] = 0.0;
        Res_Visc[nDim+1] = 0.0;
        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          Res_Visc[iDim+1] = 0.0;

        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          Res_Visc[iDim+1] = - TauWall * DirTanWM[iDim] * Area;

        Res_Visc[nDim+1] = (Wall_HeatFlux - TauWall * velWall_tan) * Area;
      }
      else{

        /*--- If it is isothermal wall, calculate the heat flux---*/
        if (!HeatFlux_Prescribed){

          su2double Prandtl_Lam  = config->GetPrandtl_Lam();
          su2double Prandtl_Turb = config->GetPrandtl_Turb();
          su2double Gas_Constant = config->GetGas_ConstantND();
          su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

          /*--- Retrieve the specified wall temperature from config
                as well as the wall function treatment.---*/
          su2double Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

          /*--- Compute closest normal neighbor ---*/
          unsigned long Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

          /*--- Get coordinates of i & nearest normal and compute distance ---*/
          su2double *Coord_i = geometry->nodes->GetCoord(iPoint);
          su2double *Coord_j = geometry->nodes->GetCoord(Point_Normal);
          su2double dist_ij = 0;
          for (iDim = 0; iDim < nDim; iDim++)
            dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
          dist_ij = sqrt(dist_ij);


          /*--- Compute the normal gradient in temperature using Twall ---*/
          su2double dTdn = -(nodes->GetTemperature(Point_Normal) - Twall)/dist_ij;

          /*--- Get transport coefficients ---*/
          su2double laminar_viscosity    = nodes->GetLaminarViscosity(iPoint);
          su2double eddy_viscosity       = nodes->GetEddyViscosity(iPoint);
          su2double thermal_conductivity = Cp * ( laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb);

          Wall_HeatFlux = thermal_conductivity * dTdn;
        }else{

          Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag) /config->GetHeat_Flux_Ref();
        }
        for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar] = 0.0;

        /*--- Weakly impose the WM heat flux for the energy equation---*/
        Res_Visc[nDim+1] = Wall_HeatFlux * Area;

      }

      LinSysRes.SubtractBlock(iPoint, Res_Visc);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  delete [] UnitNormal;

}


su2double CNSSolver::GetCHTWallTemperature(const CConfig* config, unsigned short val_marker,
                                           unsigned long iVertex, su2double thermal_conductivity,
                                           su2double dist_ij, su2double There,
                                           su2double Temperature_Ref) const {

  /*--- Compute the normal gradient in temperature using Twall ---*/

  const su2double Tconjugate = GetConjugateHeatVariable(val_marker, iVertex, 0) / Temperature_Ref;

  su2double Twall = 0.0;

  if ((config->GetKind_CHT_Coupling() == CHT_COUPLING::AVERAGED_TEMPERATURE_NEUMANN_HEATFLUX) ||
      (config->GetKind_CHT_Coupling() == CHT_COUPLING::AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)) {

    /*--- Compute wall temperature from both temperatures ---*/

    su2double HF_FactorHere = thermal_conductivity*config->GetViscosity_Ref()/dist_ij;
    su2double HF_FactorConjugate = GetConjugateHeatVariable(val_marker, iVertex, 2);

    Twall = (There*HF_FactorHere + Tconjugate*HF_FactorConjugate)/(HF_FactorHere + HF_FactorConjugate);
  }
  else if ((config->GetKind_CHT_Coupling() == CHT_COUPLING::DIRECT_TEMPERATURE_NEUMANN_HEATFLUX) ||
           (config->GetKind_CHT_Coupling() == CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX)) {

    /*--- (Directly) Set wall temperature to conjugate temperature. ---*/

    Twall = Tconjugate;
  }
  else {
    SU2_MPI::Error("Unknown CHT coupling method.", CURRENT_FUNCTION);
  }

  return Twall;
}

void CNSSolver::BC_Isothermal_Wall_Generic(CGeometry *geometry, CSolver **solver_container,
                                           CNumerics *conv_numerics, CNumerics *visc_numerics,
                                           CConfig *config, unsigned short val_marker, bool cht_mode) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const su2double Temperature_Ref = config->GetTemperature_Ref();
  const su2double Prandtl_Lam = config->GetPrandtl_Lam();
  const su2double Prandtl_Turb = config->GetPrandtl_Turb();
  const su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  /*--- Identify the boundary and retrieve the specified wall temperature from
   the config (for non-CHT problems) as well as the wall function treatment. ---*/

  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  su2double Twall = 0.0;
  if (!cht_mode) {
    Twall = config->GetIsothermal_Temperature(Marker_Tag) / Temperature_Ref;
  }

//  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
//  if (Wall_Function != WALL_FUNCTION::NONE) {
//    SU2_MPI::Error("Wall function treatment not implemented yet", CURRENT_FUNCTION);
//  }

  su2double **Jacobian_i = nullptr;
  if (implicit) {
    Jacobian_i = new su2double* [nVar];
    for (auto iVar = 0u; iVar < nVar; iVar++)
      Jacobian_i[iVar] = new su2double [nVar] ();
  }

  /*--- Loop over boundary points ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Compute dual-grid area and boundary normal ---*/

    const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

    su2double Area = GeometryToolbox::Norm(nDim, Normal);

    su2double UnitNormal[MAXNDIM] = {0.0};
    for (auto iDim = 0u; iDim < nDim; iDim++)
      UnitNormal[iDim] = -Normal[iDim]/Area;

    /*--- Compute closest normal neighbor ---*/

    const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

    /*--- Get coordinates of i & nearest normal and compute distance ---*/

    const auto Coord_i = geometry->nodes->GetCoord(iPoint);
    const auto Coord_j = geometry->nodes->GetCoord(Point_Normal);

    su2double dist_ij = GeometryToolbox::Distance(nDim, Coord_i, Coord_j);

    /*--- Store the corrected velocity at the wall which will
     be zero (v = 0), unless there is grid motion (v = u_wall)---*/

    if (dynamic_grid) {
      nodes->SetVelocity_Old(iPoint, geometry->nodes->GetGridVel(iPoint));
    }
    else {
      su2double zero[MAXNDIM] = {0.0};
      nodes->SetVelocity_Old(iPoint, zero);
    }

    for (auto iDim = 0u; iDim < nDim; iDim++)
      LinSysRes(iPoint, iDim+1) = 0.0;
    nodes->SetVel_ResTruncError_Zero(iPoint);

    /*--- Get transport coefficients ---*/

    su2double laminar_viscosity    = nodes->GetLaminarViscosity(iPoint);
    su2double eddy_viscosity       = nodes->GetEddyViscosity(iPoint);
    su2double thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb);

    // work in progress on real-gases...
    //thermal_conductivity = nodes->GetThermalConductivity(iPoint);
    //Cp = nodes->GetSpecificHeatCp(iPoint);
    //thermal_conductivity += Cp*eddy_viscosity/Prandtl_Turb;

    /*--- If it is a customizable or CHT patch, retrieve the specified wall temperature. ---*/

    const su2double There = nodes->GetTemperature(Point_Normal);

    if (cht_mode) {
      Twall = GetCHTWallTemperature(config, val_marker, iVertex, dist_ij,
                                    thermal_conductivity, There, Temperature_Ref);
    }
    else if (config->GetMarker_All_PyCustom(val_marker)) {
      Twall = geometry->GetCustomBoundaryTemperature(val_marker, iVertex);
    }

    /*--- Compute the normal gradient in temperature using Twall ---*/

    su2double dTdn = -(There - Twall)/dist_ij;

    /*--- Apply a weak boundary condition for the energy equation.
     Compute the residual due to the prescribed heat flux. ---*/

    su2double Res_Conv = 0.0;
    su2double Res_Visc = thermal_conductivity * dTdn * Area;

    /*--- Calculate Jacobian for implicit time stepping ---*/

    if (implicit) {

      /*--- Add contributions to the Jacobian from the weak enforcement of the energy equations. ---*/

      su2double Density = nodes->GetDensity(iPoint);
      su2double Vel2 = GeometryToolbox::SquaredNorm(nDim, &nodes->GetPrimitive(iPoint)[prim_idx.Velocity()]);
      su2double dTdrho = 1.0/Density * ( -Twall + (Gamma-1.0)/Gas_Constant*(Vel2/2.0) );

      Jacobian_i[nDim+1][0] = thermal_conductivity/dist_ij * dTdrho * Area;

      for (auto jDim = 0u; jDim < nDim; jDim++)
        Jacobian_i[nDim+1][jDim+1] = 0.0;

      Jacobian_i[nDim+1][nDim+1] = thermal_conductivity/dist_ij * (Gamma-1.0)/(Gas_Constant*Density) * Area;
    }

    /*--- If the wall is moving, there are additional residual contributions
     due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/

    if (dynamic_grid) {
      AddDynamicGridResidualContribution(iPoint, Point_Normal, geometry, UnitNormal,
                                         Area, geometry->nodes->GetGridVel(iPoint),
                                         Jacobian_i, Res_Conv, Res_Visc);
    }

    /*--- Convective and viscous contributions to the residual at the wall ---*/

    LinSysRes(iPoint, nDim+1) += Res_Conv - Res_Visc;

    /*--- Enforce the no-slip boundary condition in a strong way by
     modifying the velocity-rows of the Jacobian (1 on the diagonal).
     And add the contributions to the Jacobian due to energy. ---*/

    if (implicit) {
      Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

      for (auto iVar = 1u; iVar <= nDim; iVar++) {
        auto total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    }
  }
  END_SU2_OMP_FOR

  if (Jacobian_i)
    for (auto iVar = 0u; iVar < nVar; iVar++)
      delete [] Jacobian_i[iVar];
  delete [] Jacobian_i;

}

void CNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                   CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  BC_Isothermal_Wall_Generic(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CNSSolver::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                           CConfig *config, unsigned short val_marker) {
  BC_Isothermal_Wall_Generic(geometry, solver_container, conv_numerics, nullptr, config, val_marker, true);
}

void CNSSolver::SetTau_Wall_WF(CGeometry *geometry, CSolver **solver_container, const CConfig *config) {
  /*---
   The wall function implemented herein is based on Nichols and Nelson, AIAA J. v32 n6 2004.
   ---*/

  unsigned long notConvergedCounter = 0;  /*--- counts the number of wall cells that are not converged ---*/
  unsigned long smallYPlusCounter = 0;    /*--- counts the number of wall cells where y+ < 5 ---*/

  const su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  const unsigned short max_iter = config->GetwallModel_MaxIter();
  const su2double relax = config->GetwallModel_RelFac();

  /*--- Compute the recovery factor
   * use Molecular (Laminar) Prandtl number (see Nichols & Nelson, nomenclature ) ---*/

  const su2double Recovery = pow(config->GetPrandtl_Lam(), (1.0/3.0));

  /*--- Typical constants from boundary layer theory ---*/

  const su2double kappa = config->GetwallModel_Kappa();
  const su2double B = config->GetwallModel_B();

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {

    if (!config->GetViscous_Wall(iMarker)) continue;

    /*--- Identify the boundary by string name ---*/

    const auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);

    /*--- Jump to another BC if it is not wall function ---*/

    if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::STANDARD_WALL_FUNCTION)
      continue;

    /*--- Loop over all of the vertices on this boundary marker ---*/

    SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->nVertex[iMarker]; iVertex++) {

      const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      const auto Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

      /*--- Check if the node belongs to the domain (i.e, not a halo node)
       *    and the neighbor is not part of the physical boundary ---*/

      if (!geometry->nodes->GetDomain(iPoint)) continue;

      /*--- Get coordinates of the current vertex and nearest normal point ---*/

      const auto Coord = geometry->nodes->GetCoord(iPoint);
      const auto Coord_Normal = geometry->nodes->GetCoord(Point_Normal);

      /*--- Compute dual-grid area and boundary normal ---*/

      const auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

      const su2double Area = GeometryToolbox::Norm(nDim, Normal);

      su2double UnitNormal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;

      /*--- Get the velocity, pressure, and temperature at the nearest
       (normal) interior point. ---*/

      su2double Vel[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Vel[iDim] = nodes->GetVelocity(Point_Normal,iDim);

      /*--- Compute the wall-parallel velocity at first point off the wall ---*/

      const su2double VelNormal = GeometryToolbox::DotProduct(int(MAXNDIM), Vel, UnitNormal);

      su2double VelTang[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

      const su2double VelTangMod = GeometryToolbox::Norm(int(MAXNDIM), VelTang);

      /*--- Compute normal distance of the interior point from the wall ---*/

      su2double WallDist[MAXNDIM] = {0.0};
      GeometryToolbox::Distance(nDim, Coord, Coord_Normal, WallDist);

      const su2double WallDistMod = GeometryToolbox::Norm(int(MAXNDIM), WallDist);

      su2double T_Wall = nodes->GetTemperature(iPoint);
      const su2double Conductivity_Wall = nodes->GetThermalConductivity(iPoint);

      /*--- If a wall temperature was given, we compute the local heat flux using k*dT/dn ---*/

      su2double q_w = 0.0;

      if (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) {
        q_w = config->GetWall_HeatFlux(Marker_Tag) / config->GetHeat_Flux_Ref();
      }

      /*--- Extrapolate the pressure from the interior & compute the
       wall density using the equation of state ---*/

      const su2double P_Normal = nodes->GetPressure(Point_Normal);
      const su2double T_Normal = nodes->GetTemperature(Point_Normal);
      const su2double P_Wall = P_Normal;

      /*--- Compressible formulation ---*/

      su2double Density_Wall = P_Wall / (Gas_Constant * T_Wall);
      const su2double Lam_Visc_Normal = nodes->GetLaminarViscosity(Point_Normal);

      /*--- Compute the shear stress at the wall in the regular fashion
       *    by using the stress tensor on the surface ---*/

      su2double tau[MAXNDIM][MAXNDIM] = {{0.0}};
      const su2double Lam_Visc_Wall = nodes->GetLaminarViscosity(iPoint);
      su2double Eddy_Visc_Wall = nodes->GetEddyViscosity(iPoint);

      CNumerics::ComputeStressTensor(nDim, tau, nodes->GetVelocityGradient(iPoint), Lam_Visc_Wall);

      su2double TauTangent[MAXNDIM] = {0.0};
      GeometryToolbox::TangentProjection(nDim, tau, UnitNormal, TauTangent);

      const su2double WallShearStress = GeometryToolbox::Norm(int(MAXNDIM), TauTangent);

      /*--- Calculate the quantities from boundary layer theory and
       *    iteratively solve for a new wall shear stress. Use the current wall
       *    shear stress as a starting guess for the wall function. ---*/

      unsigned long counter = 0;
      su2double diff = 1.0;
      su2double U_Tau = max(1.0e-6,sqrt(WallShearStress/Density_Wall));
      /*--- Use minimum y+ as defined in the config, in case the routine below for computing y+ does not converge ---*/
      su2double Y_Plus = 0.99*config->GetwallModel_MinYPlus(); // use clipping value as minimum

      const su2double Y_Plus_Start = Density_Wall * U_Tau * WallDistMod / Lam_Visc_Wall;

      /*--- Automatic switch off when y+ < "limit" according to Nichols & Nelson (2004) ---*/

      if (Y_Plus_Start < config->GetwallModel_MinYPlus()) {
        smallYPlusCounter++;
        continue;
      }

      /*--- Convergence criterium for the Newton solver, note that 1e-10 is too large ---*/
      const su2double tol = 1e-12;
      while (fabs(diff) > tol) {

        /*--- Friction velocity and u+ ---*/

        const su2double U_Plus = VelTangMod/U_Tau;

        /*--- Gamma, Beta, Q, and Phi, defined by Nichols & Nelson (2004) page 1108 ---*/

        const su2double Gam  = Recovery*U_Tau*U_Tau/(2.0*Cp*T_Wall);
        const su2double Beta = q_w*Lam_Visc_Wall/(Density_Wall*T_Wall*Conductivity_Wall*U_Tau);
        const su2double Q    = sqrt(Beta*Beta + 4.0*Gam);
        const su2double Phi  = asin(-1.0*Beta/Q);

        /*--- Crocco-Busemann equation for wall temperature (eq. 11 of Nichols and Nelson) ---*/
        /*--- update T_Wall due to aerodynamic heating, unless the wall is isothermal      ---*/

        if (config->GetMarker_All_KindBC(iMarker) != ISOTHERMAL) {
          const su2double denum = (1.0 + Beta*U_Plus - Gam*U_Plus*U_Plus);
          if (denum > EPS){
            T_Wall = T_Normal / denum;
            nodes->SetTemperature(iPoint,T_Wall);
          }
          else {
            SU2_OMP_CRITICAL
            {
              cout << "Warning: T_Wall < 0 " << endl;
            }
            END_SU2_OMP_CRITICAL
          }
        }

        /*--- update of wall density using the wall temperature ---*/
        Density_Wall = P_Wall/(Gas_Constant*T_Wall);

        /*--- Y+ defined by White & Christoph (compressibility and heat transfer) negative value for (2.0*Gam*U_Plus - Beta)/Q ---*/

        const su2double Y_Plus_White = exp((kappa/sqrt(Gam))*(asin((2.0*Gam*U_Plus - Beta)/Q) - Phi))*exp(-1.0*kappa*B);

        /*--- Spalding's universal form for the BL velocity with the
         *    outer velocity form of White & Christoph above. ---*/
        const su2double kUp = kappa*U_Plus;
        Y_Plus = U_Plus + Y_Plus_White - (exp(-1.0*kappa*B)* (1.0 + kUp + 0.5*kUp*kUp + kUp*kUp*kUp/6.0));

        const su2double dypw_dyp = 2.0*Y_Plus_White*(kappa*sqrt(Gam)/Q)*sqrt(1.0 - pow(2.0*Gam*U_Plus - Beta,2.0)/(Q*Q));

        Eddy_Visc_Wall = Lam_Visc_Wall*(1.0 + dypw_dyp - kappa*exp(-1.0*kappa*B)*
                                         (1.0 + kappa*U_Plus + kappa*kappa*U_Plus*U_Plus/2.0)
                                         - Lam_Visc_Normal/Lam_Visc_Wall);
        Eddy_Visc_Wall = max(1.0e-6, Eddy_Visc_Wall);

        /* --- Define function for Newton method to zero --- */

        diff = (Density_Wall * U_Tau * WallDistMod / Lam_Visc_Wall) - Y_Plus;

        /* --- Gradient of function defined above --- */

        const su2double grad_diff = Density_Wall * WallDistMod / Lam_Visc_Wall + VelTangMod / (U_Tau * U_Tau) +
                  kappa /(U_Tau * sqrt(Gam)) * asin(U_Plus * sqrt(Gam)) * Y_Plus_White -
                  exp(-1.0 * B * kappa) * (0.5 * pow(VelTangMod * kappa / U_Tau, 3) +
                  pow(VelTangMod * kappa / U_Tau, 2) + VelTangMod * kappa / U_Tau) / U_Tau;

        /* --- Newton Step --- */

        U_Tau = U_Tau - relax*(diff / grad_diff);

        counter++;
        if (counter > max_iter) {
          notConvergedCounter++;
          // use some safe values for convergence
          Y_Plus = 30.0;
          Eddy_Visc_Wall = 1.0;
          U_Tau = 1.0;
          break;
        }
      }

      /*--- Calculate an updated value for the wall shear stress
       *    using the y+ value, the definition of y+, and the definition of
       *    the friction velocity. ---*/

      YPlus[iMarker][iVertex] = Y_Plus;
      EddyViscWall[iMarker][iVertex] = Eddy_Visc_Wall;
      UTau[iMarker][iVertex] = U_Tau;

      const su2double Tau_Wall = (1.0/Density_Wall)*pow(Y_Plus*Lam_Visc_Wall/WallDistMod,2.0);

      /*--- Store this value for the wall shear stress at the node.  ---*/

      nodes->SetTau_Wall(iPoint, Tau_Wall);

    }
    END_SU2_OMP_FOR
  }

  if (config->GetComm_Level() == COMM_FULL) {
    static unsigned long globalCounter1, globalCounter2;

    ompMasterAssignBarrier(globalCounter1,0, globalCounter2,0);

    SU2_OMP_ATOMIC
    globalCounter1 += notConvergedCounter;

    SU2_OMP_ATOMIC
    globalCounter2 += smallYPlusCounter;

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      SU2_MPI::Allreduce(&globalCounter1, &notConvergedCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&globalCounter2, &smallYPlusCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

      if (rank == MASTER_NODE) {
        if (notConvergedCounter)
          cout << "Warning: Computation of wall coefficients (y+) did not converge in "
               << notConvergedCounter << " points." << endl;

        if (smallYPlusCounter)
          cout << "Warning: y+ < " << config->GetwallModel_MinYPlus() << " in " << smallYPlusCounter
               << " points, for which the wall model is not active." << endl;
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

}

void CNSSolver::SetEddyViscFirstPoint(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iDim;
  unsigned long iPoint, iVertex, Point_Normal;
  su2double Lam_Visc_Normal, Vel[3] = {0.,0.,0.}, P_Normal, T_Normal;
  su2double VelTang[3] = {0.,0.,0.}, VelTangMod, VelNormal;
  su2double WallDist[3] = {0.,0.,0.}, WallDistMod;
  su2double *Coord, *Coord_Normal, Tau_Wall;
  su2double *Normal, Area, UnitNormal[3] = {0.,0.,0.};
  su2double T_Wall, P_Wall, Density_Wall, Lam_Visc_Wall;
  su2double dypw_dyp, Y_Plus_White, Eddy_Visc;
  su2double U_Plus, U_Tau, Gam, Beta, Q, Phi;
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  /*--- Typical constants from boundary layer theory ---*/

  su2double kappa = 0.41;
  su2double B = 5.0;

  /*--- Compute the recovery factor ---*/
  // su2double-check: laminar or turbulent Pr for this?
  su2double Recovery = pow(config->GetPrandtl_Lam(),(1.0/3.0));

  /* Loop over the markers and select the ones for which a wall model
    treatment is carried out. */

  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
      case HEAT_FLUX:{
        const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);

        /* Set the Eddy Viscosity at the first point off the wall */

        if ((config->GetWallFunction_Treatment(Marker_Tag) == WALL_FUNCTIONS::STANDARD_WALL_FUNCTION)){

          for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
            iPoint       = geometry->vertex[iMarker][iVertex]->GetNode();
            Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

            Coord        = geometry->nodes->GetCoord(iPoint);
            Coord_Normal = geometry->nodes->GetCoord(Point_Normal);

            /*--- Compute dual-grid area and boundary normal ---*/

            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt (Area);

            for (iDim = 0; iDim < nDim; iDim++)
              UnitNormal[iDim] = -Normal[iDim]/Area;

            /*--- Check if the node belongs to the domain (i.e, not a halo node)
             and the neighbor is not part of the physical boundary ---*/

            if (geometry->nodes->GetDomain(iPoint)) {

              Tau_Wall = nodes->GetTauWall(iPoint);

              /*--- Verify the wall function flag on the node. If false,
               jump to the next iPoint.---*/

              if (!nodes->GetTauWall_Flag(iPoint)) continue;

              /*--- Get the velocity, pressure, and temperature at the nearest
               (normal) interior point. ---*/

              for (iDim = 0; iDim < nDim; iDim++)
                Vel[iDim]    = nodes->GetVelocity(Point_Normal,iDim);
              P_Normal       = nodes->GetPressure(Point_Normal);
              T_Normal       = nodes->GetTemperature(Point_Normal);

              /*--- Compute the wall-parallel velocity at first point off the wall ---*/

              VelNormal = 0.0;
              for (iDim = 0; iDim < nDim; iDim++)
                VelNormal += Vel[iDim] * UnitNormal[iDim];
              for (iDim = 0; iDim < nDim; iDim++)
                VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

              VelTangMod = 0.0;
              for (iDim = 0; iDim < nDim; iDim++)
                VelTangMod += VelTang[iDim]*VelTang[iDim];
              VelTangMod = sqrt(VelTangMod);

              /*--- Compute normal distance of the interior point from the wall ---*/

              for (iDim = 0; iDim < nDim; iDim++)
                WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);

              WallDistMod = 0.0;
              for (iDim = 0; iDim < nDim; iDim++)
                WallDistMod += WallDist[iDim]*WallDist[iDim];
              WallDistMod = sqrt(WallDistMod);

              /*--- Compute the wall temperature using the Crocco-Buseman equation ---*/

              //T_Wall = T_Normal * (1.0 + 0.5*Gamma_Minus_One*Recovery*M_Normal*M_Normal);
              T_Wall = T_Normal + Recovery*pow(VelTangMod,2.0)/(2.0*Cp);

              /*--- Extrapolate the pressure from the interior & compute the
               wall density using the equation of state ---*/

              P_Wall       = P_Normal;
              Density_Wall = P_Wall/(Gas_Constant*T_Wall);
              Lam_Visc_Wall = nodes->GetLaminarViscosity(iPoint);

              /*--- Friction velocity and u+ ---*/

              U_Tau  = sqrt(Tau_Wall / Density_Wall);
              U_Plus = VelTangMod/U_Tau;

              Gam  = Recovery*U_Tau*U_Tau/(2.0*Cp*T_Wall);
              Beta = 0.0; // For adiabatic flows only
              Q    = sqrt(Beta*Beta + 4.0*Gam);
              Phi  = asin(-1.0*Beta/Q);

              /*--- Y+ defined by White & Christoph (compressibility and heat transfer) ---*/

              Y_Plus_White = exp((kappa/sqrt(Gam))*(asin((2.0*Gam*U_Plus - Beta)/Q) - Phi))*exp(-1.0*kappa*B);

              /*--- If the Y+ defined by White & Christoph is too high so the eddy viscosity.
                Disable the wall function calculation on this point and do not set the eddy viscosity
                at the first point off the wall. ---*/

              if (Y_Plus_White > 1e4){
                cout << "WARNING: Y+ is too high (>1e4). The muT computation at the 1st point off the wall is disable." << endl;
                cout << rank << " " << iPoint;
                for (iDim = 0; iDim < nDim; iDim++)
                  cout << " " << Coord[iDim];
                cout << endl;
                nodes->SetTauWall_Flag(iPoint,false);
                continue;
              }

              /*--- Now compute the Eddy viscosity at the first point off of the wall ---*/

              Lam_Visc_Normal = nodes->GetLaminarViscosity(Point_Normal);

              dypw_dyp = 2.0*Y_Plus_White*(kappa*sqrt(Gam)/Q)*pow(1.0 - pow(2.0*Gam*U_Plus - Beta,2.0)/(Q*Q), -0.5);

              Eddy_Visc = Lam_Visc_Wall*(1.0 + dypw_dyp - kappa*exp(-1.0*kappa*B)*
                                         (1.0 + kappa*U_Plus + kappa*kappa*U_Plus*U_Plus/2.0)
                                         - Lam_Visc_Normal/Lam_Visc_Wall);

              nodes->SetEddyViscosity(Point_Normal,Eddy_Visc);

            }
          }
        }
      }
        break;
      default:
        break;
    }
  }
}

void CNSSolver::Setmut_LES(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned long iPoint;
  su2double Grad_Vel[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  su2double lenScale, muTurb, rho;

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /* Get Density */
    rho = nodes->GetSolution(iPoint, 0);

    /* Velocity Gradients */
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        Grad_Vel[iDim][jDim] = nodes->GetGradient_Primitive(iPoint, iDim+1, jDim);

    /* Distance to the wall. */
    su2double dist = geometry->nodes->GetWall_Distance(iPoint); // Is the distance to the wall used in any SGS calculation?

    /* Length Scale for the SGS model: Cubic root of the volume. */
    su2double Vol = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
    lenScale = pow(Vol,1./3.);

    /* Compute the eddy viscosity. */
    if (nDim == 2){
      muTurb = SGSModel->ComputeEddyViscosity_2D(rho, Grad_Vel[0][0], Grad_Vel[1][0],
                                                 Grad_Vel[0][1], Grad_Vel[1][1],
                                                 lenScale, dist);
    }
    else{
      muTurb = SGSModel->ComputeEddyViscosity_3D(rho, Grad_Vel[0][0], Grad_Vel[1][0], Grad_Vel[2][0],
                                               Grad_Vel[0][1], Grad_Vel[1][1], Grad_Vel[2][1],
                                               Grad_Vel[0][2], Grad_Vel[1][2], Grad_Vel[2][2],
                                               lenScale, dist);
    }
    /* Set eddy viscosity. */
    nodes->SetEddyViscosity(iPoint, muTurb);
  }

  /*--- MPI parallelization ---*/

//  InitiateComms(geometry, config, SGS_MODEL);
//  CompleteComms(geometry, config, SGS_MODEL);

}

void CNSSolver::SetTauWallHeatFlux_WMLES1stPoint(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep) {

  /*---
  List TODO here:
   - For each vertex (point):
   - Load the interpolation coefficients.
   - Extract the LES quantities at the exchange points.
   - Call the Wall Model: Calculate Tau_Wall and Heat_Flux.
   - Set Tau_Wall and Heat_Flux in the node structure for future use.
  ---*/

  unsigned short iDim, iMarker;
  unsigned long iVertex, iPoint, Point_Normal;
  bool CalculateWallModel = false;
  bool WMLESFirstPoint = config->GetWMLES_First_Point();

  su2double Vel[3], VelNormal, VelTang[3], VelTangMod, WallDist[3], WallDistMod;
  su2double GradP[3], GradP_TangMod;
  su2double T_Normal, P_Normal, mu_Normal;
  su2double *Coord, *Coord_Normal, UnitNormal[3], *Normal, Area;
  su2double TimeFilter = config->GetDelta_UnstTimeND()/ (config->GetTimeFilter_WMLES() / config->GetTime_Ref());

  su2double Yplus_Max_Local = 0.0;
  su2double Yplus_Min_Local = 1e9;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {

     /*--- Identify the boundary by string name ---*/
     string Marker_Tag = config->GetMarker_All_TagBound(iMarker);


     /*--- Identify if this marker is a wall model one---*/
     /*
     switch (config->GetWallFunction_Treatment(Marker_Tag)) {
       case WALL_FUNCTIONS::EQUILIBRIUM_MODEL:
       case WALL_FUNCTIONS::LOGARITHMIC_MODEL:
       case WALL_FUNCTIONS::NONEQUILIBRIUM_MODEL:
       case WALL_FUNCTIONS::ADAPTIVE_FUNCTION:
         break;

       case WALL_FUNCTIONS::NONE:
       case WALL_FUNCTIONS::STANDARD_FUNCTION:
         CalculateWallModel = false;
       default:
         break;
     }
    */
     /*--- Identify if this marker is a wall model one---*/
     switch (config->GetWallFunction_Treatment(Marker_Tag)) {
       case WALL_FUNCTIONS::EQUILIBRIUM_WALL_MODEL:
       case WALL_FUNCTIONS::LOGARITHMIC_WALL_MODEL:
       case WALL_FUNCTIONS::ALGEBRAIC_WALL_MODEL:
       case WALL_FUNCTIONS::APGLL_WALL_MODEL:
       case WALL_FUNCTIONS::TEMPLATE_WALL_MODEL:
         CalculateWallModel = true;
         break;

       case WALL_FUNCTIONS::NO_WALL_FUNCTION:
       case WALL_FUNCTIONS::STANDARD_WALL_FUNCTION:
         CalculateWallModel = false;
       default:
         break;
     }

     /*--- If not just continue to the next---*/
     if (!CalculateWallModel) continue;

     /*--- Determine the prescribed heat flux or prescribed temperature. ---*/
     bool HeatFlux_Prescribed = false, Temperature_Prescribed = false;
     su2double Wall_HeatFlux = 0.0, Wall_Temperature = 0.0;

     if(config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) {
       HeatFlux_Prescribed = true;
       Wall_HeatFlux       = config->GetWall_HeatFlux(Marker_Tag);
     }
     else {
       Temperature_Prescribed = true;
       Wall_Temperature       = config->GetIsothermal_Temperature(Marker_Tag);
     }

     /*--- Loop over all of the vertices on this boundary marker ---*/
     for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

       iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
       Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

       /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
       if (!geometry->nodes->GetDomain(iPoint)) continue;

       /*--- Check if the node has all boundary neighbors ---*/
       const auto nNeigh = geometry->nodes->GetnPoint(iPoint);
       bool found_fluid = false;
       for (unsigned short iNeigh = 0; iNeigh <= nNeigh; iNeigh++) {
       
         auto jPoint = iPoint;
         if (iNeigh < nNeigh) jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
         bool boundary_j = geometry->nodes->GetPhysicalBoundary(jPoint);
         if (!boundary_j){
           found_fluid = true;
           break;
         }
       }

       if(!found_fluid){
        nodes->SetTauWall_Flag(iPoint,false);
        continue;
       }

       /*--- Get coordinates of the current vertex and nearest normal point ---*/

       Coord = geometry->nodes->GetCoord(iPoint);
       Coord_Normal = geometry->nodes->GetCoord(Point_Normal);

       /*--- Compute dual-grid area and boundary normal ---*/

       Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

       Area = 0.0;
       for (iDim = 0; iDim < nDim; iDim++)
         Area += Normal[iDim]*Normal[iDim];
       Area = sqrt (Area);

       for (iDim = 0; iDim < nDim; iDim++)
         UnitNormal[iDim] = -Normal[iDim]/Area;

       /*--- If an exchange location was found (donor element) use this information as the input
       for the wall model. Otherwise, the information of the 1st point off the wall is used. ---*/

       if (geometry->vertex[iMarker][iVertex]->GetDonorFound() && !WMLESFirstPoint){

         const su2double *doubleInfo = config->GetWallFunction_DoubleInfo(Marker_Tag);
         WallDistMod = doubleInfo[0];

         /*--- Load the coefficients and interpolate---*/
         unsigned short nDonors = geometry->vertex[iMarker][iVertex]->GetnDonorPoints();
         su2double rho_Normal = 0.0, e_Normal   = 0.0;

         for (iDim = 0; iDim < nDim; iDim++ ){
           Vel[iDim]   = 0.0;
         }

         for (unsigned short iNode = 0; iNode < nDonors; iNode++) {
           unsigned long donorPoint = geometry->vertex[iMarker][iVertex]->GetInterpDonorPoint(iNode);
           su2double donnorCoeff    = geometry->vertex[iMarker][iVertex]->GetDonorCoeff(iNode);

           rho_Normal += donnorCoeff*nodes->GetSolution(donorPoint,0);
           e_Normal   += donnorCoeff*nodes->GetSolution(donorPoint,nVar-1)/nodes->GetSolution(donorPoint,0);

           for (iDim = 0; iDim < nDim; iDim++ ){
             Vel[iDim] += donnorCoeff*nodes->GetSolution(donorPoint,iDim+1)/nodes->GetSolution(donorPoint,0);
           }
         }

         /*--- Load the fluid model to have pressure and temperature at exchange location---*/
         GetFluidModel()->SetTDState_rhoe(rho_Normal, e_Normal);
         P_Normal = GetFluidModel()->GetPressure();
         T_Normal = GetFluidModel()->GetTemperature();
         mu_Normal= GetFluidModel()->GetLaminarViscosity();

       }
       else{
 
        /*--- Compute normal distance of the interior point from the wall ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);

        WallDistMod = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          WallDistMod += WallDist[iDim]*WallDist[iDim];
        WallDistMod = sqrt(WallDistMod);

        /*--- Get the velocity, pressure, and temperature at the nearest
        (normal) interior point. ---*/

        for (iDim = 0; iDim < nDim; iDim++){
          Vel[iDim]   = nodes->GetVelocity(Point_Normal,iDim);
        }

        P_Normal  = nodes->GetPressure(Point_Normal);
        T_Normal  = nodes->GetTemperature(Point_Normal);
        mu_Normal = nodes->GetLaminarViscosity(Point_Normal);
       }

       for (iDim = 0; iDim < nDim; iDim++){
         GradP[iDim] = nodes->GetGradient_Primitive(iPoint, nDim+1, iDim);
       }

       /*--- Filter the input LES velocity ---*/

       long curAbsTimeIter = (config->GetTimeIter() - config->GetRestart_Iter());
       if (curAbsTimeIter > 0){

         /*--- Old input LES velocity and GradP---*/
         su2double Vel_old[3]   = {0.,0.,0.};
         for (iDim = 0; iDim < nDim; iDim++){
           Vel_old[iDim]   = VelTimeFilter_WMLES[iMarker][iVertex][iDim];
         }
         /*--- Now filter the LES velocity ---*/
         for (iDim = 0; iDim < nDim; iDim++){
           Vel[iDim] = (1.0 - TimeFilter) * Vel_old[iDim] + TimeFilter * Vel[iDim];
         }
       }

       /*--- Update input LES velocity if it is the 1st inner iteration---*/
       if (config->GetInnerIter() == 0){
         for (iDim = 0; iDim < nDim; iDim++){
           VelTimeFilter_WMLES[iMarker][iVertex][iDim] = Vel[iDim];
         }
       }

       /*--- Compute dimensional variables before calling the Wall Model ---*/
       for (iDim = 0; iDim < nDim; iDim++ ){
         Vel[iDim] *= config->GetVelocity_Ref();
         GradP[iDim] *= config->GetPressure_Ref();
       }
       P_Normal *= config->GetPressure_Ref();
       T_Normal *= config->GetTemperature_Ref();
       mu_Normal *= (config->GetPressure_Ref()/config->GetVelocity_Ref());

       /*--- Compute the wall-parallel velocity ---*/

       VelNormal = 0.0;
       for (iDim = 0; iDim < nDim; iDim++)
         VelNormal += Vel[iDim] * UnitNormal[iDim];
       for (iDim = 0; iDim < nDim; iDim++)
         VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

       VelTangMod = 0.0;
       for (iDim = 0; iDim < nDim; iDim++)
         VelTangMod += VelTang[iDim]*VelTang[iDim];
       VelTangMod = sqrt(VelTangMod);
       VelTangMod = max(VelTangMod,1.e-25);

       su2double dirTan[3] = {0.0, 0.0, 0.0};
       for(iDim = 0; iDim<nDim; iDim++) dirTan[iDim] = VelTang[iDim]/VelTangMod;

       /*--- If it is pressure gradient driven flow
        subtract the body force in all directions. ---*/
       if (config->GetBody_Force()){
         for (iDim = 0; iDim < nDim; iDim++)
           GradP[iDim] -= config->GetBody_Force_Vector()[iDim];
       }

       /*--- Pressure gradient in the tangent direction: ---*/
       GradP_TangMod = 0.0;
       for (iDim = 0; iDim < nDim; iDim++)
         GradP_TangMod += GradP[iDim]*dirTan[iDim];

       /* Compute the wall shear stress and heat flux vector using
        the wall model. */
       su2double tauWall, qWall, ViscosityWall, kOverCvWall;
       bool converged;
       WallModel->UpdateExchangeLocation(WallDistMod);
       WallModel->WallShearStressAndHeatFlux(T_Normal, VelTangMod, mu_Normal, P_Normal, GradP_TangMod,
                                             Wall_HeatFlux, HeatFlux_Prescribed,
                                             Wall_Temperature, Temperature_Prescribed,
                                             GetFluidModel(), tauWall, qWall, ViscosityWall,
                                             kOverCvWall, converged);
       if (!converged || std::isnan(tauWall)){
          nodes->SetTauWall_Flag(iPoint,false);
         continue;
       }
       
       su2double rho    = nodes->GetDensity(iPoint) * config->GetDensity_Ref();
       su2double u_tau  = sqrt(tauWall/rho);
       su2double y_plus = rho * u_tau * WallDistMod / ViscosityWall; 

 
       if (config->GetWMLES_Monitoring()){
        if(y_plus < 0.1){
         nodes->SetTauWall_Flag(iPoint,false);
         continue;          
        }
       }

       Yplus_Max_Local = max(Yplus_Max_Local, y_plus);
       Yplus_Min_Local = min(Yplus_Min_Local, y_plus);

       /*--- Compute the non-dimensional values if necessary. ---*/
       tauWall /= config->GetPressure_Ref();
       qWall   /= (config->GetPressure_Ref() * config->GetVelocity_Ref());
       ViscosityWall /= (config->GetPressure_Ref()/config->GetVelocity_Ref());
       nodes->SetLaminarViscosity(iPoint, ViscosityWall);


       /*--- Set tau wall value and flag for flux computation---*/
       nodes->SetTauWall_Flag(iPoint,true);
       nodes->SetTauWall(iPoint, tauWall);

       /*--- Set tau wall projected to the flow direction for pos-processing only---*/
       for(iDim = 0; iDim<nDim; iDim++)
         nodes->SetTauWallDir(iPoint, iDim, tauWall*dirTan[iDim]);


       /*--- Set tau wall value and heat flux for boundary conditions---*/
       TauWall_WMLES[iMarker][iVertex] = tauWall;
       HeatFlux_WMLES[iMarker][iVertex] = qWall;
       for (iDim = 0; iDim < nDim; iDim++)
         FlowDirTan_WMLES[iMarker][iVertex][iDim] = dirTan[iDim];

     }
   }
 }
 su2double Yplus_Max_Global = Yplus_Max_Local;
 su2double Yplus_Min_Global = Yplus_Min_Local;

 SU2_MPI::Allreduce(&Yplus_Max_Local, &Yplus_Max_Global, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
 SU2_MPI::Allreduce(&Yplus_Min_Local, &Yplus_Min_Global, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());

 if ((rank == MASTER_NODE) && (config->GetInnerIter()==0)){
  cout << endl   << "------------------------ WMLES -----------------------" << endl;
  cout << "Y+ (Max): " << setprecision(6) << Yplus_Max_Global << endl;
  cout << "Y+ (Min): " << setprecision(6) << Yplus_Min_Global << endl;
 }

}
