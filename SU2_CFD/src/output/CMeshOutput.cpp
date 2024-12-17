/*!
 * \file CMeshOutput.cpp
 * \brief Main subroutines for the heat solver output
 * \author R. Sanchez
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


#include "../../include/output/CMeshOutput.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"

CMeshOutput::CMeshOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {

  /*--- Set the default history fields if nothing is set in the config file ---*/

  requestedVolumeFields.emplace_back("COORDINATES");
  nRequestedVolumeFields = requestedVolumeFields.size();

  /*--- Set the volume filename --- */

  volumeFilename = config->GetMesh_Out_FileName();

  /*--- Set the surface filename ---*/

  surfaceFilename = "surface_mesh";

}

CMeshOutput::~CMeshOutput(void) {}

void CMeshOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

  // Mesh quality metrics, computed in CPhysicalGeometry::ComputeMeshQualityStatistics.
  AddVolumeOutput("ORTHOGONALITY", "Orthogonality", "MESH_QUALITY", "Orthogonality Angle (deg.)");
  AddVolumeOutput("ASPECT_RATIO",  "Aspect_Ratio",  "MESH_QUALITY", "CV Face Area Aspect Ratio");
  AddVolumeOutput("VOLUME_RATIO",  "Volume_Ratio",  "MESH_QUALITY", "CV Sub-Volume Ratio");

  if (config->GetQuickSurfOut()){
    AddVolumeOutput("GLOBAL_ID", "Global_id", "PRIMITIVE", "Global_id");
  }

}

void CMeshOutput::SetProbeOutputFields(CConfig *config, unsigned int nProbe){

  // Grid coordinates
  AddProbeOutput(nProbe, "COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddProbeOutput(nProbe, "COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddProbeOutput(nProbe, "COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

  // Mesh quality metrics, computed in CPhysicalGeometry::ComputeMeshQualityStatistics.
  AddProbeOutput(nProbe, "ORTHOGONALITY", "Orthogonality", "MESH_QUALITY", "Orthogonality Angle (deg.)");
  AddProbeOutput(nProbe, "ASPECT_RATIO",  "Aspect_Ratio",  "MESH_QUALITY", "CV Face Area Aspect Ratio");
  AddProbeOutput(nProbe, "VOLUME_RATIO",  "Volume_Ratio",  "MESH_QUALITY", "CV Sub-Volume Ratio");

}

void CMeshOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CPoint*    Node_Geo  = geometry->nodes;

  if (config->GetQuickSurfOut()){
    SetVolumeOutputValue("GLOBAL_ID", iPoint, Node_Geo->GetGlobalIndex(iPoint));
  }

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(iPoint, 2));

  // Mesh quality metrics
  if (config->GetWrt_MeshQuality()) {
    SetVolumeOutputValue("ORTHOGONALITY", iPoint, geometry->Orthogonality[iPoint]);
    SetVolumeOutputValue("ASPECT_RATIO",  iPoint, geometry->Aspect_Ratio[iPoint]);
    SetVolumeOutputValue("VOLUME_RATIO",  iPoint, geometry->Volume_Ratio[iPoint]);
  }

}

void CMeshOutput::LoadProbeData(CConfig *config, CGeometry *geometry, CSolver **solver){

  auto probe_list = geometry->GetProbe_list();
  if(probe_list.size()>0 && probe_list[0].rankID==rank){
    CPoint*    Node_Geo  = geometry->nodes;
    unsigned long iPoint{0};
    vector<unsigned long> Probe_pointID = geometry->GetProbe_pointID();

    for(unsigned int index=0; index<probe_list.size(); index++){
      iPoint = probe_list[index].pointID;

      SetProbeOutputValue("COORD-X", index,  Node_Geo->GetCoord(iPoint, 0));
      SetProbeOutputValue("COORD-Y", index,  Node_Geo->GetCoord(iPoint, 1));
      if (nDim == 3)
        SetProbeOutputValue("COORD-Z", index, Node_Geo->GetCoord(iPoint, 2));

      // Mesh quality metrics
      if (config->GetWrt_MeshQuality()) {
        SetProbeOutputValue("ORTHOGONALITY", index, geometry->Orthogonality[iPoint]);
        SetProbeOutputValue("ASPECT_RATIO",  index, geometry->Aspect_Ratio[iPoint]);
        SetProbeOutputValue("VOLUME_RATIO",  index, geometry->Volume_Ratio[iPoint]);
      }
    }
  }
}
