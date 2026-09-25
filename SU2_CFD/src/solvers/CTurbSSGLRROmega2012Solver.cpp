/*!
 * \file CTurbSSGLRROmega2012Solver.cpp
 * \brief Main subroutines of CTurbSSGLRROmega2012Solver class
 * \author F. Palacios, A. Bueno
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/solvers/CTurbSSGLRROmega2012Solver.hpp"
#include "../../include/variables/CTurbSSGLRROmega2012Variable.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"


CTurbSSGLRROmega2012Solver::CTurbSSGLRROmega2012Solver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver(geometry, config, true) {
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;

  bool multizone = config->GetMultizone_Problem();

  /*--- Dimension of the problem --> dependent on the turbulence model. ---*/

  nVar = 7;
  nPrimVar = 7;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0 || config->GetMGCycle() == FULLMG_CYCLE) {

    /*--- Define some auxiliary vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SSG/LRR-RSM-omega2012)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    System.SetxIsZero(true);

    if (ReducerStrategy) {
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
      EdgeFluxesDiff.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
    }

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS.resize(nVar,0.0);
      Residual_Max_BGS.resize(nVar,0.0);
      Point_Max_BGS.resize(nVar,0);
      Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
    }

  }

  /*--- Far-field flow state quantities and initialization. ---*/
  su2double rhoInf, *VelInf, muLamInf, Intensity, viscRatio, muT_Inf;

  rhoInf    = config->GetDensity_FreeStreamND();
  VelInf    = config->GetVelocity_FreeStreamND();
  muLamInf  = config->GetViscosity_FreeStreamND();
  Intensity = config->GetTurbulenceIntensity_FreeStream();
  viscRatio = config->GetTurb2LamViscRatio_FreeStream();

  su2double VelMag2 = GeometryToolbox::SquaredNorm(nDim, VelInf);

  su2double kine_Inf  = 3.0/2.0*(VelMag2*Intensity*Intensity);
  su2double omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);

  for (unsigned short v=0; v<7; ++v) Solution_Inf[v]=0.0;
  Solution_Inf[0]=Solution_Inf[1]=Solution_Inf[2]=(2.0/3.0)*kine_Inf;
  Solution_Inf[6]=omega_Inf;

  /*--- Constants to use for lower limit of turbulence variable. ---*/
  su2double Ck = config->GetKFactor_LowerLimit();
  su2double Cw = config->GetOmegaFactor_LowerLimit();

  /*--- Initialize lower and upper limits. ---*/
  for (unsigned short v=0; v<7; ++v) {
    lowerlimit[v]= v<3 ? 0.0 : -1.0e15;
    upperlimit[v]= 1.0e15;
  }
  lowerlimit[6]=1.0e-4;

  /*--- Eddy viscosity, initialized without stress limiter at the infinity ---*/
  muT_Inf = rhoInf*kine_Inf/omega_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CTurbSSGLRROmega2012Variable(kine_Inf, omega_Inf, muT_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION_EDDY);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION_EDDY);

  /*--- Initialize quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
    due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_TurbVars[iMarker].resize(nVertex[iMarker],nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      for (unsigned short v=0; v<7; ++v)
        Inlet_TurbVars[iMarker](iVertex,v)=Solution_Inf[v];
    }
  }

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Turb();
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name. ---*/
  SolverName = "SSGLRR_OMEGA2012";

}

void CTurbSSGLRROmega2012Solver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
         unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)

  /*--- Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);
}

void CTurbSSGLRROmega2012Solver::Postprocessing(CGeometry *geometry, CSolver **solver_container,
                                                CConfig *config, unsigned short iMesh) {
  flowState=solver_container[FLOW_SOL]->GetNodes();
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
    SetSolution_Gradient_GG(geometry, config, -1);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
    SetSolution_Gradient_LS(geometry, config, -1);
  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint=0; iPoint<nPoint; ++iPoint) {
    const su2double rho=flowNodes->GetDensity(iPoint);
    const su2double mu=flowNodes->GetLaminarViscosity(iPoint);
    const su2double k=0.5*(nodes->GetSolution(iPoint,0)+nodes->GetSolution(iPoint,1)
                           +nodes->GetSolution(iPoint,2));
    const su2double omega=nodes->GetSolution(iPoint,6);
    const su2double d=geometry->nodes->GetWall_Distance(iPoint);
    if (d>0.0 && k>0.0 && omega>0.0)
      nodes->SetBlendingFunc(iPoint,mu,d,rho,TURB_TRANS_MODEL::NONE);
    // Only the turbulent heat-conductivity surrogate uses muT. The momentum
    // flux must use the six transported stresses, never this Boussinesq value.
    nodes->SetmuT(iPoint,rho*k/omega);
  }
  END_SU2_OMP_FOR
}

void CTurbSSGLRROmega2012Solver::Viscous_Residual(const unsigned long iEdge, const CGeometry* geometry, CSolver** solver_container,
                                     CNumerics* numerics, const CConfig* config) {

  /*--- Define an object to set solver specific numerics contribution. ---*/
  auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) {
    /*--- Menter's first blending function (only SST)---*/
    numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(jPoint));
  };

  /*--- Now instantiate the generic non-conservative implementation with the functor above. ---*/
  Viscous_Residual_NonCons(iEdge, geometry, solver_container, numerics, config, SolverSpecificNumerics);

}

void CTurbSSGLRROmega2012Solver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  bool axisymmetric = config->GetAxisymmetric();

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  /*--- Pick one numerics object per thread. ---*/
  auto* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Loop over all points. ---*/

  AD::StartNoSharedReading();

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/

    numerics->SetScalarVar(nodes->GetSolution(iPoint), nullptr);
    numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Set distance to the surface ---*/

    numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

    /*--- Menter's first blending function ---*/

    numerics->SetF1blending(nodes->GetF1blending(iPoint),0.0);

    /*--- Menter's second blending function ---*/

    numerics->SetF2blending(nodes->GetF2blending(iPoint));

    /*--- Set vorticity and strain rate magnitude ---*/

    numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);

    numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);

    /*--- Cross diffusion ---*/

    numerics->SetCrossDiff(nodes->GetCrossDiff(iPoint));

    /*--- Effective Intermittency ---*/
    if (config->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      numerics->SetIntermittencyEff(solver_container[TRANS_SOL]->GetNodes()->GetIntermittencyEff(iPoint));
    }

    if (axisymmetric){
      /*--- Set y coordinate ---*/
      numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(iPoint));
    }

    /*--- Compute the source term ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Store the intermittency ---*/

    if (config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE) {
      nodes->SetIntermittency(iPoint, numerics->GetIntermittencyEff());
    }

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, residual);
    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

  /*--- Custom user defined source term (from the python wrapper) ---*/
  if (config->GetPyCustomSource()) {
    CustomSourceResidual(geometry, solver_container, numerics_container, config, iMesh);
  }

}

void CTurbSSGLRROmega2012Solver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {
}

void CTurbSSGLRROmega2012Solver::BC_HeatFlux_Wall(
    CGeometry *geometry, CSolver **solver_container, CNumerics*, CNumerics*,
    CConfig *config, unsigned short val_marker) {
  if (config->GetWall_Functions())
    SU2_MPI::Error("SSG/LRR wall functions are not implemented; use resolved walls",CURRENT_FUNCTION);
  const bool implicit=config->GetKind_TimeIntScheme()==EULER_IMPLICIT;
  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex=0u; iVertex<geometry->nVertex[val_marker]; ++iVertex) {
    const auto iPoint=geometry->vertex[val_marker][iVertex]->GetNode();
    if (!geometry->nodes->GetDomain(iPoint)) continue;
    const su2double d1=geometry->vertex[val_marker][iVertex]->GetNearestNeighborDistance();
    const su2double rho=solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    const su2double mu=solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    su2double wall[7]={};
    wall[6]=60.0*mu/(rho*0.075*d1*d1); // NASA: 10*6*nu/(beta_omega*d1^2)
    nodes->SetSolution_Old(iPoint,wall);
    nodes->SetSolution(iPoint,wall);
    LinSysRes.SetBlock_Zero(iPoint);
    if (implicit)
      for (unsigned short v=0; v<7; ++v) Jacobian.DeleteValsRowi(iPoint*nVar+v);
  }
  END_SU2_OMP_FOR
}

void CTurbSSGLRROmega2012Solver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTurbSSGLRROmega2012Solver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                              CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/

      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      su2double Inlet_Vars[MAXNVAR];
      if (config->GetInlet_Profile_From_File()) {
        for (unsigned short v=0; v<6; ++v)
          Inlet_Vars[v]=Inlet_TurbVars[val_marker][iVertex][v]/pow(config->GetVelocity_Ref(),2);
        Inlet_Vars[6]=Inlet_TurbVars[val_marker][iVertex][6]*config->GetViscosity_Ref()
                       /(config->GetDensity_Ref()*pow(config->GetVelocity_Ref(),2));
      } else {
        /*--- Obtain fluid model for computing the  kine and omega to impose at the inlet boundary. ---*/
        CFluidModel* FluidModel = solver_container[FLOW_SOL]->GetFluidModel();

        /*--- Obtain flow velocity vector at inlet boundary node ---*/

        const su2double* Velocity_Inlet = &V_inlet[prim_idx.Velocity()];
        su2double Density_Inlet;
        if (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {
          Density_Inlet = V_inlet[prim_idx.Density()];
          FluidModel->SetTDState_Prho(V_inlet[prim_idx.Pressure()], Density_Inlet);
        } else {
          const su2double* Scalar_Inlet = nullptr;
          if (config->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
            Scalar_Inlet = config->GetInlet_SpeciesVal(config->GetMarker_All_TagBound(val_marker));
          }
          FluidModel->SetTDState_T(V_inlet[prim_idx.Temperature()], Scalar_Inlet);
          Density_Inlet = FluidModel->GetDensity();
        }
        const su2double Laminar_Viscosity_Inlet = FluidModel->GetLaminarViscosity();
        const su2double* Turb_Properties = config->GetInlet_TurbVal(config->GetMarker_All_TagBound(val_marker));
        const su2double Intensity = Turb_Properties[0];
        const su2double viscRatio = Turb_Properties[1];
        const su2double VelMag2 = GeometryToolbox::SquaredNorm(nDim, Velocity_Inlet);

        const su2double k=1.5*VelMag2*Intensity*Intensity;
        for (unsigned short v=0; v<7; ++v) Inlet_Vars[v]=0.0;
        Inlet_Vars[0]=Inlet_Vars[1]=Inlet_Vars[2]=(2.0/3.0)*k;
        Inlet_Vars[6]=Density_Inlet*k/(Laminar_Viscosity_Inlet*viscRatio);
      }

      /*--- Positive covariance and omega at inflow. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Inlet_Vars);

      /*--- Set various other quantities in the solver class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      if (conv_numerics->GetBoundedScalar()) {
        const su2double* velocity = &V_inlet[prim_idx.Velocity()];
        const su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        conv_numerics->SetMassFlux(BoundedScalarBCFlux(iPoint, implicit, density, velocity, Normal));
      }

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
      //
      //      su2double Coord_Reflected[MAXNDIM];
      //      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
      //                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      //      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      //      visc_numerics->SetNormal(Normal);
      //
      //      /*--- Conservative variables w/o reconstruction ---*/
      //
      //      visc_numerics->SetPrimitive(V_domain, V_inlet);
      //
      //      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      //
      //     visc_numerics->SetScalarVar(Solution_i, Solution_j);
      //     visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      //
      //      /*--- Menter's first blending function ---*/
      //
      //      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
      //
      //      /*--- Compute residual, and Jacobians ---*/
      //
      //      auto residual = visc_numerics->ComputeResidual(config);
      //
      //      /*--- Subtract residual, and update Jacobians ---*/
      //
      //      LinSysRes.SubtractBlock(iPoint, residual);
      //      Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

    }

  }
  END_SU2_OMP_FOR
}

void CTurbSSGLRROmega2012Solver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                               CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Allocate the value at the outlet ---*/

      auto V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint),
                                nodes->GetSolution(iPoint));

      /*--- Set Normal (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
      conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                geometry->nodes->GetGridVel(iPoint));

      if (conv_numerics->GetBoundedScalar()) {
        const su2double* velocity = &V_outlet[prim_idx.Velocity()];
        const su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        conv_numerics->SetMassFlux(BoundedScalarBCFlux(iPoint, implicit, density, velocity, Normal));
      }

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      su2double Coord_Reflected[MAXNDIM];
//      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
//                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetScalarVar(Solution_i, Solution_j);
//      visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- Menter's first blending function ---*/
//
//      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      auto residual = visc_numerics->ComputeResidual(config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, residual);
//      Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR
}


void CTurbSSGLRROmega2012Solver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                          CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (auto iSpan = 0u; iSpan < nSpanWiseSections ; iSpan++){

    su2double extAverageKine = solver_container[FLOW_SOL]->GetExtAverageKine(val_marker, iSpan);
    su2double extAverageOmega = solver_container[FLOW_SOL]->GetExtAverageOmega(val_marker, iSpan);
    su2double solution_j[] = {extAverageKine, extAverageOmega};

    /*--- Loop over all the vertices on this boundary marker ---*/

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      const auto iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const auto oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const auto Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][oldVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/
      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/
      auto conv_residual = conv_numerics->ComputeResidual(config);

      /*--- Jacobian contribution for implicit integration ---*/
      LinSysRes.AddBlock(iPoint, conv_residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, conv_residual.jacobian_i);

      /*--- Viscous contribution ---*/
      su2double Coord_Reflected[MAXNDIM];
      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);
      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

      /*--- Compute residual, and Jacobians ---*/
      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, visc_residual);
      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

}

void CTurbSSGLRROmega2012Solver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Quantities for computing the  kine and omega to impose at the inlet boundary. ---*/
  CFluidModel *FluidModel = solver_container[FLOW_SOL]->GetFluidModel();

  su2double Intensity = config->GetTurbulenceIntensity_FreeStream();
  su2double viscRatio = config->GetTurb2LamViscRatio_FreeStream();

  for (auto iSpan = 0u; iSpan < nSpanWiseSections ; iSpan++){

    /*--- Compute the inflow kine and omega using the span wise averge quntities---*/

    su2double rho       = solver_container[FLOW_SOL]->GetAverageDensity(val_marker, iSpan);
    su2double pressure  = solver_container[FLOW_SOL]->GetAveragePressure(val_marker, iSpan);
    su2double kine      = solver_container[FLOW_SOL]->GetAverageKine(val_marker, iSpan);

    FluidModel->SetTDState_Prho(pressure, rho);
    su2double muLam = FluidModel->GetLaminarViscosity();

    su2double VelMag2 = GeometryToolbox::SquaredNorm(nDim,
                          solver_container[FLOW_SOL]->GetAverageTurboVelocity(val_marker, iSpan));

    su2double kine_b  = 3.0/2.0*(VelMag2*Intensity*Intensity);
    su2double omega_b = rho*kine_b/(muLam*viscRatio);

    su2double solution_j[7] = {(2.0/3.0)*kine_b,(2.0/3.0)*kine_b,
                               (2.0/3.0)*kine_b,0.0,0.0,0.0,omega_b};

    /*--- Loop over all the vertices on this boundary marker ---*/

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      const auto iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const auto oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const auto Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][oldVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/
      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states. Use average span-wise values
             values for the turbulent state at the inflow. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/
      auto conv_residual = conv_numerics->ComputeResidual(config);

      /*--- Jacobian contribution for implicit integration ---*/
      LinSysRes.AddBlock(iPoint, conv_residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, conv_residual.jacobian_i);

      /*--- Viscous contribution ---*/
      su2double Coord_Reflected[MAXNDIM];
      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);

      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

      /*--- Compute residual, and Jacobians ---*/
      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, visc_residual);
      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

}

void CTurbSSGLRROmega2012Solver::SetInletAtVertex(const su2double *val_inlet,
                                     unsigned short iMarker,
                                     unsigned long iVertex) {

  for (unsigned short v=0; v<7; ++v)
    Inlet_TurbVars[iMarker][iVertex][v]=val_inlet[nDim+2+nDim+v];

}

su2double CTurbSSGLRROmega2012Solver::GetInletAtVertex(unsigned short iMarker, unsigned long iVertex,
                                           const CGeometry* geometry, su2double* val_inlet) const {
  const auto rsm_position = nDim + 2 + nDim;
  for (unsigned short v=0; v<7; ++v)
    val_inlet[rsm_position+v]=Inlet_TurbVars[iMarker][iVertex][v];

  /*--- Compute boundary face area for this vertex. ---*/

  su2double Normal[MAXNDIM] = {0.0};
  geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
  return GeometryToolbox::Norm(nDim, Normal);}

void CTurbSSGLRROmega2012Solver::SetUniformInlet(const CConfig* config, unsigned short iMarker) {
  if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      for (unsigned short v=0; v<7; ++v)
        Inlet_TurbVars[iMarker][iVertex][v]=Solution_Inf[v];
    }
  }

}

void CTurbSSGLRROmega2012Solver::ComputeUnderRelaxationFactor(const CConfig *config) {
  const su2double allowableRatio = config->GetMaxUpdateFractionSST();
  if (!flowState) SU2_MPI::Error("Missing flow state in SSG/LRR update",CURRENT_FUNCTION);
  unsigned long limited=0, cancelled=0;
  su2double minFactor=1.0;
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint=0; iPoint<nPointDomain; ++iPoint) {
    const su2double rho=flowState->GetDensity(iPoint);
    const su2double k=0.5*(nodes->GetSolution(iPoint,0)+nodes->GetSolution(iPoint,1)
                           +nodes->GetSolution(iPoint,2));
    const su2double omega=nodes->GetSolution(iPoint,6);
    const su2double stressScale=max(k,su2double(1e-20));
    const su2double omegaScale=max(omega,su2double(1e-20));
    su2double factor=1.0;
    for (unsigned short v=0; v<7; ++v) {
      const su2double increment=LinSysSol(iPoint,v)/rho;
      const su2double scale=v==6 ? omegaScale : stressScale;
      const su2double ratio=fabs(increment)/scale;
      if (ratio>allowableRatio) factor=min(factor,allowableRatio/ratio);
    }
    // A single factor preserves the tensor update direction. This checks all
    // principal minors of the symmetric 3x3 covariance before accepting it.
    bool admissible=false;
    for (unsigned short attempt=0; attempt<60; ++attempt) {
      su2double r[6];
      for (unsigned short v=0; v<6; ++v)
        r[v]=nodes->GetSolution(iPoint,v)+factor*LinSysSol(iPoint,v)/rho;
      const su2double w=omega+factor*LinSysSol(iPoint,6)/rho;
      const su2double scale=max(k,su2double(1e-20));
      const su2double tol2=1e-12*scale*scale;
      const su2double tol3=1e-12*scale*scale*scale;
      const su2double det=r[0]*r[1]*r[2]+2*r[3]*r[4]*r[5]
                          -r[0]*r[5]*r[5]-r[1]*r[4]*r[4]-r[2]*r[3]*r[3];
      if (w>0.0 && r[0]>=0.0 && r[1]>=0.0 && r[2]>=0.0 &&
          r[0]*r[1]-r[3]*r[3]>=-tol2 &&
          r[0]*r[2]-r[4]*r[4]>=-tol2 &&
          r[1]*r[2]-r[5]*r[5]>=-tol2 && det>=-tol3) {
        admissible=true;
        break;
      }
      factor*=0.5;
    }
    if (!admissible || factor<1e-10) factor=0.0;
    nodes->SetUnderRelaxation(iPoint,factor);
  }
  END_SU2_OMP_FOR
  for (unsigned long iPoint=0; iPoint<nPointDomain; ++iPoint) {
    const su2double factor=nodes->GetUnderRelaxation(iPoint);
    if (factor<1.0) ++limited;
    if (factor==0.0) ++cancelled;
    minFactor=min(minFactor,factor);
  }
  if (config->GetInnerIter()<10 || config->GetInnerIter()%100==0) {
    unsigned long globalLimited=0,globalCancelled=0;
    double localMin=SU2_TYPE::GetValue(minFactor),globalMin=1.0;
    SU2_MPI::Allreduce(&limited,&globalLimited,1,MPI_UNSIGNED_LONG,MPI_SUM,SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&cancelled,&globalCancelled,1,MPI_UNSIGNED_LONG,MPI_SUM,SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&localMin,&globalMin,1,MPI_DOUBLE,MPI_MIN,SU2_MPI::GetComm());
    if (rank==MASTER_NODE)
      cout << "SSGLRR shared update: limited=" << globalLimited << " cancelled=" << globalCancelled
           << " minFactor=" << globalMin << endl;
  }
}
