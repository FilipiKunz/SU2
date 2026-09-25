/*!
 * \file turb_convection.hpp
 * \brief Declarations of numerics classes for discretization of
 *        convective fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
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

#pragma once

#include "../scalar/scalar_convection.hpp"

/*!
 * \class CUpwSca_TurbSA
 * \brief Class for doing a scalar upwind solver for the Spalar-Allmaras turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 */
template <class FlowIndices>
class CUpwSca_TurbSA final : public CUpwScalar<FlowIndices> {
private:
  using Base = CUpwScalar<FlowIndices>;
  using Base::a0;
  using Base::a1;
  using Base::Flux;
  using Base::Jacobian_i;
  using Base::Jacobian_j;
  using Base::ScalarVar_i;
  using Base::ScalarVar_j;
  using Base::bounded_scalar;

  /*!
   * \brief Adds any extra variables to AD.
   */
  void ExtraADPreaccIn() override {}

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override {
    Flux[0] = a0*ScalarVar_i[0] + a1*ScalarVar_j[0];
    Jacobian_i[0][0] = a0;
    Jacobian_j[0][0] = a1;
  }

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
    : CUpwScalar<FlowIndices>(val_nDim, val_nVar, config) { bounded_scalar = config->GetBounded_Turb(); }
};

/*!
 * \class CUpwSca_TurbSST
 * \brief Class for doing a scalar upwind solver for the Menter SST turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Campos.
 */
template <class FlowIndices>
class CUpwSca_TurbSST final : public CUpwScalar<FlowIndices> {
private:
  using Base = CUpwScalar<FlowIndices>;
  using Base::nDim;
  using Base::V_i;
  using Base::V_j;
  using Base::a0;
  using Base::a1;
  using Base::Flux;
  using Base::Jacobian_i;
  using Base::Jacobian_j;
  using Base::ScalarVar_i;
  using Base::ScalarVar_j;
  using Base::idx;
  using Base::bounded_scalar;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override {}

  /*!
   * \brief SST specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override {
    Flux[0] = a0*V_i[idx.Density()]*ScalarVar_i[0] + a1*V_j[idx.Density()]*ScalarVar_j[0];
    Flux[1] = a0*V_i[idx.Density()]*ScalarVar_i[1] + a1*V_j[idx.Density()]*ScalarVar_j[1];

    Jacobian_i[0][0] = a0;    Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;   Jacobian_i[1][1] = a0;

    Jacobian_j[0][0] = a1;    Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;   Jacobian_j[1][1] = a1;
  }

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
    : CUpwScalar<FlowIndices>(val_nDim, val_nVar, config) { bounded_scalar = config->GetBounded_Turb(); }
};

/* Six Favre stresses (11,22,33,12,13,23) and omega, transported conservatively. */
template <class FlowIndices>
class CUpwSca_TurbSSGLRR final : public CUpwScalar<FlowIndices> {
  using Base = CUpwScalar<FlowIndices>;
  using Base::a0;
  using Base::a1;
  using Base::V_i;
  using Base::V_j;
  using Base::idx;
  using Base::ScalarVar_i;
  using Base::ScalarVar_j;
  using Base::Flux;
  using Base::Jacobian_i;
  using Base::Jacobian_j;
  using Base::bounded_scalar;

  void ExtraADPreaccIn() override {}
  void FinishResidualCalc(const CConfig*) override {
    for (unsigned short v=0; v<7; ++v) {
      Flux[v] = a0*V_i[idx.Density()]*ScalarVar_i[v]
              + a1*V_j[idx.Density()]*ScalarVar_j[v];
      for (unsigned short w=0; w<7; ++w) {
        Jacobian_i[v][w] = v==w ? a0 : 0.0;
        Jacobian_j[v][w] = v==w ? a1 : 0.0;
      }
    }
  }

public:
  CUpwSca_TurbSSGLRR(unsigned short ndim, unsigned short nvar, const CConfig* config)
    : CUpwScalar<FlowIndices>(ndim,nvar,config) {
    if (nvar != 7) SU2_MPI::Error("SSG/LRR needs seven equations", CURRENT_FUNCTION);
    bounded_scalar = config->GetBounded_Turb();
  }
};
