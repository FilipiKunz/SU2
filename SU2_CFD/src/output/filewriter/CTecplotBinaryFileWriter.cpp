/*!
 * \file CTecplotBinaryFileWriter.cpp
 * \brief Filewriter class for Tecplot binary format.
 * \author T. Albring
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

#include "../../../include/output/filewriter/CTecplotBinaryFileWriter.hpp"
#ifdef HAVE_TECIO
  #include "TECIO.h"
#endif
#ifdef HAVE_SERIAL_TECIO
  #include "../../../externals/tecio/serial_writer.hpp"
#endif
#include <algorithm>
#include <limits>

const string CTecplotBinaryFileWriter::fileExt = ".szplt";

CTecplotBinaryFileWriter::CTecplotBinaryFileWriter(CParallelDataSorter *valDataSorter,
                                                   unsigned long valTimeIter, su2double valTimeStep) :
  CFileWriter(valDataSorter, fileExt), timeIter(valTimeIter), timeStep(valTimeStep){}

CTecplotBinaryFileWriter::~CTecplotBinaryFileWriter()= default;

void CTecplotBinaryFileWriter::WriteData(string val_filename){

  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);

  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  /*--- Set a timer for the binary file writing. ---*/

  startTime = SU2_MPI::Wtime();

#ifdef HAVE_TECIO

  const vector<string> fieldNames = dataSorter->GetFieldNames();

  /*--- Reduce the total number of each element. ---*/

  unsigned long nParallel_Line = dataSorter->GetnElem(LINE),
                nParallel_Tria = dataSorter->GetnElem(TRIANGLE),
                nParallel_Quad = dataSorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = dataSorter->GetnElem(TETRAHEDRON),
                nParallel_Hexa = dataSorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = dataSorter->GetnElem(PRISM),
                nParallel_Pyra = dataSorter->GetnElem(PYRAMID);

  unsigned long nTot_Line = dataSorter->GetnElemGlobal(LINE),
                nTot_Tria = dataSorter->GetnElemGlobal(TRIANGLE),
                nTot_Quad = dataSorter->GetnElemGlobal(QUADRILATERAL),
                nTot_Tetr = dataSorter->GetnElemGlobal(TETRAHEDRON),
                nTot_Hexa = dataSorter->GetnElemGlobal(HEXAHEDRON),
                nTot_Pris = dataSorter->GetnElemGlobal(PRISM),
                nTot_Pyra = dataSorter->GetnElemGlobal(PYRAMID);

  string data_set_title = "Visualization of the solution";

  ostringstream tecplot_variable_names;
  for (size_t iVar = 0; iVar < fieldNames.size()-1; ++iVar) {
    tecplot_variable_names << fieldNames[iVar] << ",";
  }
  tecplot_variable_names << fieldNames[fieldNames.size()-1];

  /*--- Determine the Tecplot zone type from the mesh elements. ---*/

  int64_t num_nodes;
  int64_t num_cells;
  int32_t zone_type;


  num_nodes = static_cast<int64_t>(dataSorter->GetnPointsGlobal());
  num_cells = static_cast<int64_t>(dataSorter->GetnElemGlobal());
  if (dataSorter->GetnDim() == 3){
    if ((nTot_Quad > 0 || nTot_Tria > 0) && (nTot_Hexa + nTot_Pris + nTot_Pyra + nTot_Tetr == 0)){
      zone_type = ZONETYPE_FEQUADRILATERAL;
    }
    else {
      zone_type = ZONETYPE_FEBRICK;
    }
  }
  else {
    if (nTot_Line > 0 && (nTot_Tria + nTot_Quad == 0)){
      zone_type = ZONETYPE_FELINESEG;

    }
    else{
      zone_type = ZONETYPE_FEQUADRILATERAL;
    }
  }

  void* file_handle = nullptr;
  int32_t err = 0;
  /*--- TecIO partitions 3D brick zones. For other zones, only the master
        opens the shared output file. ---*/
#ifdef HAVE_MPI
  if (zone_type == ZONETYPE_FEBRICK || rank == MASTER_NODE)
#endif
  {
#ifdef HAVE_SERIAL_TECIO
    if (zone_type != ZONETYPE_FEBRICK)
      err = SU2TecOpen(val_filename.c_str(), data_set_title.c_str(), tecplot_variable_names.str().c_str(), &file_handle);
    else
#endif
      err = tecFileWriterOpen(val_filename.c_str(), data_set_title.c_str(), tecplot_variable_names.str().c_str(),
        FILEFORMAT_SZL, FILETYPE_FULL, (int32_t)FieldDataType_Float, nullptr, &file_handle);
    if (err) SU2_MPI::Error("Error opening Tecplot file " + val_filename, CURRENT_FUNCTION);
  }
#ifdef HAVE_MPI
  if (zone_type == ZONETYPE_FEBRICK) {
    err = tecMPIInitialize(file_handle, SU2_MPI::GetComm(), MASTER_NODE);
    if (err) SU2_MPI::Error("Error initializing Tecplot parallel output.", CURRENT_FUNCTION);
  }
#endif

  bool is_unsteady = false;
  passivedouble solution_time = 0.0;

  if (timeStep > 0.0){
    is_unsteady = true;
    solution_time = SU2_TYPE::GetValue(timeStep)*timeIter;
  }

  int32_t zone = 0;
  vector<int32_t> value_locations(fieldNames.size(), 1); /* Nodal variables. */
  if (file_handle) {
#ifdef HAVE_SERIAL_TECIO
    if (zone_type != ZONETYPE_FEBRICK)
      err = SU2TecZone(file_handle, "Zone", zone_type, num_nodes, num_cells, value_locations.data(), &zone);
    else
#endif
      err = tecZoneCreateFE(file_handle, "Zone", zone_type, num_nodes, num_cells, nullptr, nullptr,
                            value_locations.data(), nullptr, 0, 0, 0, &zone);
    if (err) SU2_MPI::Error("Error creating Tecplot zone.", CURRENT_FUNCTION);
    if (is_unsteady) {
#ifdef HAVE_SERIAL_TECIO
      if (zone_type != ZONETYPE_FEBRICK)
        err = SU2TecTime(file_handle, zone, solution_time, timeIter + 1);
      else
#endif
        err = tecZoneSetUnsteadyOptions(file_handle, zone, solution_time, timeIter + 1);
      if (err) SU2_MPI::Error("Error setting Tecplot zone time.", CURRENT_FUNCTION);
    }
  }

#ifdef HAVE_MPI

  unsigned short iVar;
  unsigned long localPoints = dataSorter->GetnPoints();
  vector<unsigned long> pointCounts(size), nodeOffsets(size + 1, 0);
  SU2_MPI::Allgather(&localPoints, 1, MPI_UNSIGNED_LONG, pointCounts.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
  for (int iRank = 0; iRank < size; ++iRank)
    nodeOffsets[iRank + 1] = nodeOffsets[iRank] + pointCounts[iRank];
  const unsigned long localBegin = nodeOffsets[rank];
  const unsigned long localEnd = nodeOffsets[rank + 1];
  if (nodeOffsets.back() != num_nodes)
    SU2_MPI::Error("Tecplot point counts do not match the zone.", CURRENT_FUNCTION);

  vector<unsigned long> halo_nodes;
  vector<unsigned long> sorted_halo_nodes;
  vector<float> halo_var_data;
  vector<int> num_nodes_to_receive(size, 0);
  vector<int> values_to_receive_displacements(size);
  size_t num_halo_nodes = 0;

  if (zone_type == ZONETYPE_FEBRICK) {
    /* Each rank writes its own part of the single Tecplot file. */
    vector<int32_t> partition_owners;
    partition_owners.reserve(size);
    for (int32_t iRank = 0; iRank < size; ++iRank)
      partition_owners.push_back(iRank);
    err = tecZoneMapPartitionsToMPIRanks(file_handle, zone, size, partition_owners.data());
    if (err) SU2_MPI::Error("Error assigning MPI ranks for Tecplot zone partitions.", CURRENT_FUNCTION);

    /* Gather a list of nodes we refer to but are not outputting. */

    for (unsigned long i = 0; i < nParallel_Line * N_POINTS_LINE; ++i)
      if ((unsigned long)dataSorter->GetElemConnectivity(LINE, 0, i) <= localBegin ||
          localEnd < (unsigned long)dataSorter->GetElemConnectivity(LINE, 0, i))
        halo_nodes.push_back(dataSorter->GetElemConnectivity(LINE, 0, i));

    for (unsigned long i = 0; i < nParallel_Tria * N_POINTS_TRIANGLE; ++i)
      if ((unsigned long)dataSorter->GetElemConnectivity(TRIANGLE, 0, i) <= localBegin ||
          localEnd < (unsigned long)dataSorter->GetElemConnectivity(TRIANGLE, 0, i))
        halo_nodes.push_back(dataSorter->GetElemConnectivity(TRIANGLE, 0, i));

    for (unsigned long i = 0; i < nParallel_Quad * N_POINTS_QUADRILATERAL; ++i)
      if ((unsigned long)dataSorter->GetElemConnectivity(QUADRILATERAL, 0, i) <= localBegin ||
          localEnd < (unsigned long)dataSorter->GetElemConnectivity(QUADRILATERAL, 0, i))
        halo_nodes.push_back(dataSorter->GetElemConnectivity(QUADRILATERAL, 0, i));

    for (unsigned long i = 0; i < nParallel_Tetr * N_POINTS_TETRAHEDRON; ++i)
      if ((unsigned long)dataSorter->GetElemConnectivity(TETRAHEDRON, 0, i) <= localBegin ||
          localEnd < (unsigned long)dataSorter->GetElemConnectivity(TETRAHEDRON, 0, i))
        halo_nodes.push_back(dataSorter->GetElemConnectivity(TETRAHEDRON, 0, i));

    for (unsigned long i = 0; i < nParallel_Hexa * N_POINTS_HEXAHEDRON; ++i)
      if ((unsigned long)dataSorter->GetElemConnectivity(HEXAHEDRON, 0, i) <= localBegin ||
          localEnd < (unsigned long)dataSorter->GetElemConnectivity(HEXAHEDRON, 0, i))
        halo_nodes.push_back(dataSorter->GetElemConnectivity(HEXAHEDRON, 0, i));

    for (unsigned long i = 0; i < nParallel_Pris * N_POINTS_PRISM; ++i)
      if ((unsigned long)dataSorter->GetElemConnectivity(PRISM, 0, i) <= localBegin ||
          localEnd < (unsigned long)dataSorter->GetElemConnectivity(PRISM, 0, i))
        halo_nodes.push_back(dataSorter->GetElemConnectivity(PRISM, 0, i));

    for (unsigned long i = 0; i < nParallel_Pyra * N_POINTS_PYRAMID; ++i)
      if ((unsigned long)dataSorter->GetElemConnectivity(PYRAMID, 0, i) <= localBegin ||
          localEnd < (unsigned long)dataSorter->GetElemConnectivity(PYRAMID, 0, i))
        halo_nodes.push_back(dataSorter->GetElemConnectivity(PYRAMID, 0, i));

    /* Sorted list of halo nodes for this MPI rank. */
    std::sort(halo_nodes.begin(), halo_nodes.end());
    halo_nodes.erase(std::unique(halo_nodes.begin(), halo_nodes.end()), halo_nodes.end());
    sorted_halo_nodes.swap(halo_nodes);

    /* Have to include all nodes our cells refer to or TecIO will barf, so add the halo node count to the number of local nodes. */
    int64_t partition_num_nodes = localEnd - localBegin + static_cast<int64_t>(sorted_halo_nodes.size());
    int64_t partition_num_cells = dataSorter->GetnElem();

    /*--- We effectively tack the halo nodes onto the end of the node list for this partition.
      TecIO will later replace them with references to nodes in neighboring partitions. */
    num_halo_nodes = sorted_halo_nodes.size();
    vector<int64_t> halo_node_local_numbers(max((size_t)1, num_halo_nodes)); /* Min size 1 to avoid crashes when we access these vectors below. */
    vector<int32_t> neighbor_partitions(max((size_t)1, num_halo_nodes));
    vector<int64_t> neighbor_nodes(max((size_t)1, num_halo_nodes));
    for(int64_t i = 0; i < static_cast<int64_t>(num_halo_nodes); ++i) {
      halo_node_local_numbers[i] = localEnd - localBegin + i + 1;
      const unsigned long node = sorted_halo_nodes[i];
      const int owning_rank = std::upper_bound(nodeOffsets.begin(), nodeOffsets.end(), node - 1) - nodeOffsets.begin() - 1;
      const unsigned long node_number = node - nodeOffsets[owning_rank];
      neighbor_partitions[i] = owning_rank + 1; /* Partition numbers are 1-based. */
      neighbor_nodes[i] = static_cast<int64_t>(node_number);
    }
    err = tecFEPartitionCreate64(file_handle, zone, rank + 1, partition_num_nodes, partition_num_cells,
      static_cast<int64_t>(num_halo_nodes), halo_node_local_numbers.data(), neighbor_partitions.data(), neighbor_nodes.data(), 0, nullptr);
    if (err) SU2_MPI::Error("Error creating Tecplot zone partition.", CURRENT_FUNCTION);

    /* Gather halo node data. First, tell each rank how many nodes' worth of data we need from them. */
    for (size_t i = 0; i < num_halo_nodes; ++i)
      ++num_nodes_to_receive[neighbor_partitions[i] - 1];
    vector<int> num_nodes_to_send(size);
    SU2_MPI::Alltoall(num_nodes_to_receive.data(), 1, MPI_INT, num_nodes_to_send.data(), 1, MPI_INT, SU2_MPI::GetComm());

    /* Now send the global node numbers whose data we need,
       and receive the same from all other ranks.
       Each rank has globally consecutive node numbers,
       so we can just parcel out sorted_halo_nodes for send. */
    vector<int> nodes_to_send_displacements(size);
    vector<int> nodes_to_receive_displacements(size);
    nodes_to_send_displacements[0] = 0;
    nodes_to_receive_displacements[0] = 0;
    for(int iRank = 1; iRank < size; ++iRank) {
      nodes_to_send_displacements[iRank] = nodes_to_send_displacements[iRank - 1] + num_nodes_to_send[iRank - 1];
      nodes_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank - 1] + num_nodes_to_receive[iRank - 1];
    }
    int total_num_nodes_to_send = nodes_to_send_displacements[size - 1] + num_nodes_to_send[size - 1];
    vector<unsigned long> nodes_to_send(max(1, total_num_nodes_to_send));

    /* The terminology gets a bit confusing here. We're sending the node numbers
       (sorted_halo_nodes) whose data we need to receive, and receiving
       lists of nodes whose data we need to send. */
    if (sorted_halo_nodes.empty()) sorted_halo_nodes.resize(1); /* Avoid crash. */
    SU2_MPI::Alltoallv(sorted_halo_nodes.data(), num_nodes_to_receive.data(), nodes_to_receive_displacements.data(), MPI_UNSIGNED_LONG,
                       nodes_to_send.data(),     num_nodes_to_send.data(),    nodes_to_send_displacements.data(),    MPI_UNSIGNED_LONG,
                       SU2_MPI::GetComm());

    /* Now actually send and receive the data */
    vector<float> data_to_send(max(1, total_num_nodes_to_send * (int)fieldNames.size()));
    halo_var_data.resize(max((size_t)1, fieldNames.size() * num_halo_nodes));
    vector<int> num_values_to_send(size);
    vector<int> values_to_send_displacements(size);
    vector<int> num_values_to_receive(size);
    size_t index = 0;
    for(int iRank = 0; iRank < size; ++iRank) {
      /* We send and receive GlobalField_Counter values per node. */
      num_values_to_send[iRank]              = num_nodes_to_send[iRank] * fieldNames.size();
      values_to_send_displacements[iRank]    = nodes_to_send_displacements[iRank] * fieldNames.size();
      num_values_to_receive[iRank]           = num_nodes_to_receive[iRank] * fieldNames.size();
      values_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank] * fieldNames.size();
      for(iVar = 0; iVar < fieldNames.size(); ++iVar)
        for(int iNode = 0; iNode < num_nodes_to_send[iRank]; ++iNode) {
          unsigned long node_offset = nodes_to_send[nodes_to_send_displacements[iRank] + iNode] - localBegin - 1;
          data_to_send[index++] =dataSorter->GetData(iVar,node_offset);
        }
    }
    CBaseMPIWrapper::Alltoallv(data_to_send.data(),  num_values_to_send.data(),    values_to_send_displacements.data(),    MPI_FLOAT,
                       halo_var_data.data(), num_values_to_receive.data(), values_to_receive_displacements.data(), MPI_FLOAT,
                       SU2_MPI::GetComm());
  }
  /*--- Write surface and volumetric solution data. ---*/

  {
    const size_t localNodes = dataSorter->GetnPoints();
    const bool brick = zone_type == ZONETYPE_FEBRICK;
    const bool gather = !brick && num_nodes <= std::numeric_limits<int>::max();
    vector<int> counts, offsets;
    if (gather) {
      counts.resize(size);
      offsets.resize(size);
      for (int iRank = 0; iRank < size; ++iRank) {
        counts[iRank] = pointCounts[iRank];
        offsets[iRank] = nodeOffsets[iRank];
      }
    }
    std::vector<float> values_to_write(brick || gather ? localNodes : std::min<size_t>(localNodes, 1 << 20));
    std::vector<float> recv_values(brick || rank != MASTER_NODE ? 1 : std::min<size_t>(num_nodes, 1 << 20));
    std::vector<float> global_values(gather && rank == MASTER_NODE ? num_nodes : 1);
    for (iVar = 0; err == 0 && iVar < fieldNames.size(); iVar++) {
      if (brick) {
        for (size_t i = 0; i < localNodes; ++i)
          values_to_write[i] = dataSorter->GetData(iVar, i);
        err = tecZoneVarWriteFloatValues(file_handle, zone, iVar + 1, rank + 1, localNodes, values_to_write.data());
        if (err) SU2_MPI::Error("Error outputting Tecplot variable values.", CURRENT_FUNCTION);
        for (int iRank = 0; iRank < size; ++iRank) {
          const int count = num_nodes_to_receive[iRank];
          if (count > 0) {
            const int offset = values_to_receive_displacements[iRank] + count * iVar;
            err = tecZoneVarWriteFloatValues(file_handle, zone, iVar + 1, rank + 1, count,
                                               halo_var_data.data() + offset);
            if (err) SU2_MPI::Error("Error outputting Tecplot halo values.", CURRENT_FUNCTION);
          }
        }
      } else {
        if (gather) {
          for (size_t i = 0; i < localNodes; ++i)
            values_to_write[i] = dataSorter->GetData(iVar, i);
          MPI_Gatherv(values_to_write.data(), localNodes, MPI_FLOAT, global_values.data(), counts.data(),
                      offsets.data(), MPI_FLOAT, MASTER_NODE, SU2_MPI::GetComm());
        } else {
          for (int iRank = 0; iRank < size; ++iRank) {
            for (size_t start = 0; start < pointCounts[iRank]; start += 1 << 20) {
              const int count = std::min<size_t>(1 << 20, pointCounts[iRank] - start);
              if (rank == iRank) {
                for (int i = 0; i < count; ++i)
                  values_to_write[i] = dataSorter->GetData(iVar, start + i);
                if (rank != MASTER_NODE)
                  CBaseMPIWrapper::Send(values_to_write.data(), count, MPI_FLOAT, MASTER_NODE, iRank, SU2_MPI::GetComm());
              }
              if (rank == MASTER_NODE) {
                float* values = values_to_write.data();
                if (iRank != MASTER_NODE) {
                  CBaseMPIWrapper::Recv(recv_values.data(), count, MPI_FLOAT, iRank, iRank, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
                  values = recv_values.data();
                }
#ifdef HAVE_SERIAL_TECIO
                err = SU2TecVar(file_handle, zone, iVar + 1, count, values);
#else
                err = tecZoneVarWriteFloatValues(file_handle, zone, iVar + 1, 0, count, values);
#endif
                if (err) SU2_MPI::Error("Error outputting Tecplot variable values.", CURRENT_FUNCTION);
              }
            }
          }
        }
        if (gather && rank == MASTER_NODE) {
#ifdef HAVE_SERIAL_TECIO
          err = SU2TecVar(file_handle, zone, iVar + 1, num_nodes, global_values.data());
#else
          err = tecZoneVarWriteFloatValues(file_handle, zone, iVar + 1, 0, num_nodes, global_values.data());
#endif
          if (err) SU2_MPI::Error("Error outputting Tecplot variable values.", CURRENT_FUNCTION);
        }
      }
    }
  }

#else

  unsigned short iVar;

  vector<float> var_data;
  size_t var_data_size = fieldNames.size() * dataSorter->GetnPoints();
  var_data.reserve(var_data_size);


  for (iVar = 0; err == 0 && iVar <  fieldNames.size(); iVar++) {
    for(unsigned long i = 0; i < dataSorter->GetnPoints(); ++i)
      var_data.push_back(dataSorter->GetData(iVar,i));
    err = tecZoneVarWriteFloatValues(file_handle, zone, iVar + 1, 0, dataSorter->GetnPoints(), &var_data[iVar * dataSorter->GetnPoints()]);
    if (err) SU2_MPI::Error("Error outputting Tecplot variable value.", CURRENT_FUNCTION);
  }


#endif /* HAVE_MPI */

  /*--- Write connectivity data. ---*/

  unsigned long iElem;

  /*--- TecIO accepts a node count, so write bounded blocks of cells instead
   * of calling it once for every element. ---*/
  vector<int64_t> nodeBuffer;
  nodeBuffer.reserve(65536);
  int32_t bufferPartition = 0;
  int32_t bufferWidth = 0;
  auto flushNodes = [&]() -> int32_t {
    if (nodeBuffer.empty()) return 0;
    const int32_t status = tecZoneNodeMapWrite64(file_handle, zone, bufferPartition, 1,
                                                  nodeBuffer.size(), nodeBuffer.data());
    nodeBuffer.clear();
    return status;
  };
  auto writeNodes = [&](const int64_t* nodes, int32_t width, int32_t partition) -> int32_t {
    if (!nodeBuffer.empty() && (partition != bufferPartition || width != bufferWidth ||
                                nodeBuffer.size() + width > 65536)) {
      const int32_t status = flushNodes();
      if (status) return status;
    }
    bufferPartition = partition;
    bufferWidth = width;
    nodeBuffer.insert(nodeBuffer.end(), nodes, nodes + width);
    return 0;
  };

#ifdef HAVE_MPI
  if (zone_type == ZONETYPE_FEBRICK) {

    int64_t nodes[8];

    /**
     *  Each rank writes node numbers relative to the partition it is outputting (starting with node number 1).
     *  Ghost (halo) nodes identified above are numbered sequentially just beyond the end of the actual, local nodes.
     *  Note that beg_node and end_node refer to zero-based node numbering, but Conn_* contain one-based node numbers.
     */
#define MAKE_LOCAL(n) localBegin < (unsigned long)n && (unsigned long)n <= localEnd \
  ? (int64_t)((unsigned long)n - localBegin) \
  : GetHaloNodeNumber(n, localEnd - localBegin, sorted_halo_nodes)

    for (iElem = 0; err == 0 && iElem < nParallel_Tetr; iElem++) {
      nodes[0] = MAKE_LOCAL(dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 0));
      nodes[1] = MAKE_LOCAL(dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 1));
      nodes[2] = MAKE_LOCAL(dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 2));
      nodes[3] = nodes[2];
      nodes[4] = MAKE_LOCAL(dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 3));
      nodes[5] = nodes[4];
      nodes[6] = nodes[4];
      nodes[7] = nodes[4];
      err = writeNodes(nodes, 8, rank + 1);
      if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Hexa; iElem++) {
      nodes[0] = MAKE_LOCAL(dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 0));
      nodes[1] = MAKE_LOCAL(dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 1));
      nodes[2] = MAKE_LOCAL(dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 2));
      nodes[3] = MAKE_LOCAL(dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 3));
      nodes[4] = MAKE_LOCAL(dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 4));
      nodes[5] = MAKE_LOCAL(dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 5));
      nodes[6] = MAKE_LOCAL(dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 6));
      nodes[7] = MAKE_LOCAL(dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 7));
      err = writeNodes(nodes, 8, rank + 1);
      if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Pris; iElem++) {
      nodes[0] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PRISM, iElem, 0));
      nodes[1] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PRISM, iElem, 1));
      nodes[2] = nodes[1];
      nodes[3] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PRISM, iElem, 2));
      nodes[4] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PRISM, iElem, 3));
      nodes[5] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PRISM, iElem, 4));
      nodes[6] = nodes[5];
      nodes[7] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PRISM, iElem, 5));
      err = writeNodes(nodes, 8, rank + 1);
      if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Pyra; iElem++) {
      nodes[0] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PYRAMID, iElem, 0));
      nodes[1] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PYRAMID, iElem, 1));
      nodes[2] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PYRAMID, iElem, 2));
      nodes[3] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PYRAMID, iElem, 3));
      nodes[4] = MAKE_LOCAL(dataSorter->GetElemConnectivity(PYRAMID, iElem, 4));
      nodes[5] = nodes[4];
      nodes[6] = nodes[4];
      nodes[7] = nodes[4];
      err = writeNodes(nodes, 8, rank + 1);
      if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
    }
  } else {
    /* TecIO cannot partition line or quadrilateral zones. Stream their data
       through the master into one zone and one file. */
    unsigned long localCounts[3] = {nParallel_Line, nParallel_Tria, nParallel_Quad};
    vector<unsigned long> allCounts(3 * size);
    SU2_MPI::Allgather(localCounts, 3, MPI_UNSIGNED_LONG, allCounts.data(), 3, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
    vector<int64_t> cells(65536);
    const GEO_TYPE types[3] = {LINE, TRIANGLE, QUADRILATERAL};
    for (int iRank = 0; iRank < size; ++iRank) {
      for (int iType = 0; iType < 3; ++iType) {
        const int width = iType == 0 ? 2 : 4;
        const unsigned long count = allCounts[3 * iRank + iType];
        const unsigned long blockCells = cells.size() / width;
        for (unsigned long start = 0; start < count; start += blockCells) {
          const unsigned long chunk = std::min(blockCells, count - start);
          if (rank == iRank) {
            for (unsigned long i = 0; i < chunk; ++i) {
              for (int j = 0; j < width; ++j) {
                const int source = iType == 1 && j == 3 ? 2 : j;
                cells[i * width + j] = dataSorter->GetElemConnectivity(types[iType], start + i, source);
              }
            }
            if (rank != MASTER_NODE)
              CBaseMPIWrapper::Send(cells.data(), chunk * width, MPI_LONG_LONG, MASTER_NODE, iRank, SU2_MPI::GetComm());
          }
          if (rank == MASTER_NODE) {
            if (iRank != MASTER_NODE)
              CBaseMPIWrapper::Recv(cells.data(), chunk * width, MPI_LONG_LONG, iRank, iRank,
                                    SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
#ifdef HAVE_SERIAL_TECIO
            err = SU2TecNodes(file_handle, zone, chunk * width, cells.data());
#else
            err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, chunk * width, cells.data());
#endif
            if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
          }
        }
      }
    }
  }

#else

  int64_t nodes[8];

  for (iElem = 0; err == 0 && iElem < nParallel_Line; iElem++) {
    nodes[0] = dataSorter->GetElemConnectivity(LINE, iElem, 0);
    nodes[1] = dataSorter->GetElemConnectivity(LINE, iElem, 1);
    err = writeNodes(nodes, 2, rank);
  }

  for (iElem = 0; err == 0 && iElem < nParallel_Tria; iElem++) {
    nodes[0] = dataSorter->GetElemConnectivity(TRIANGLE, iElem, 0);
    nodes[1] = dataSorter->GetElemConnectivity(TRIANGLE, iElem, 1);
    nodes[2] = dataSorter->GetElemConnectivity(TRIANGLE, iElem, 2);
    nodes[3] = dataSorter->GetElemConnectivity(TRIANGLE, iElem, 2);
    err = writeNodes(nodes, 4, rank);
    if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
  }

  for (iElem = 0; err == 0 && iElem < nParallel_Quad; iElem++) {
    nodes[0] = dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 0);
    nodes[1] = dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 1);
    nodes[2] = dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 2);
    nodes[3] = dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 3);
    err = writeNodes(nodes, 4, rank);
    if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
  }

  for (iElem = 0; err == 0 && iElem < nParallel_Tetr; iElem++) {
    nodes[0] = dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 0);
    nodes[1] = dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 1);
    nodes[2] = dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 2);
    nodes[3] = dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 2);
    nodes[4] = dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 3);
    nodes[5] = dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 3);
    nodes[6] = dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 3);
    nodes[7] = dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 3);
    err = writeNodes(nodes, 8, rank);
    if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
  }

  for (iElem = 0; err == 0 && iElem < nParallel_Hexa; iElem++) {
    nodes[0] = dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 0);
    nodes[1] = dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 1);
    nodes[2] = dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 2);
    nodes[3] = dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 3);
    nodes[4] = dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 4);
    nodes[5] = dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 5);
    nodes[6] = dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 6);
    nodes[7] = dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 7);
    err = writeNodes(nodes, 8, rank);
    if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
  }

  for (iElem = 0; err == 0 && iElem < nParallel_Pris; iElem++) {
    nodes[0] = dataSorter->GetElemConnectivity(PRISM, iElem, 0);
    nodes[1] = dataSorter->GetElemConnectivity(PRISM, iElem, 1);
    nodes[2] = dataSorter->GetElemConnectivity(PRISM, iElem, 1);
    nodes[3] = dataSorter->GetElemConnectivity(PRISM, iElem, 2);
    nodes[4] = dataSorter->GetElemConnectivity(PRISM, iElem, 3);
    nodes[5] = dataSorter->GetElemConnectivity(PRISM, iElem, 4);
    nodes[6] = dataSorter->GetElemConnectivity(PRISM, iElem, 4);
    nodes[7] = dataSorter->GetElemConnectivity(PRISM, iElem, 5);
    err = writeNodes(nodes, 8, rank);
    if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
  }

  for (iElem = 0; err == 0 && iElem < nParallel_Pyra; iElem++) {
    nodes[0] = dataSorter->GetElemConnectivity(PYRAMID, iElem, 0);
    nodes[1] = dataSorter->GetElemConnectivity(PYRAMID, iElem, 1);
    nodes[2] = dataSorter->GetElemConnectivity(PYRAMID, iElem, 2);
    nodes[3] = dataSorter->GetElemConnectivity(PYRAMID, iElem, 3);
    nodes[4] = dataSorter->GetElemConnectivity(PYRAMID, iElem, 4);
    nodes[5] = dataSorter->GetElemConnectivity(PYRAMID, iElem, 4);
    nodes[6] = dataSorter->GetElemConnectivity(PYRAMID, iElem, 4);
    nodes[7] = dataSorter->GetElemConnectivity(PYRAMID, iElem, 4);
    err = writeNodes(nodes, 8, rank);
    if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);
  }



#endif

  if (err == 0) err = flushNodes();
  if (err) SU2_MPI::Error("Error outputting Tecplot node values.", CURRENT_FUNCTION);

  if (file_handle) {
#ifdef HAVE_SERIAL_TECIO
    if (zone_type != ZONETYPE_FEBRICK)
      err = SU2TecClose(&file_handle);
    else
#endif
    err = tecFileWriterClose(&file_handle);
    if (err) SU2_MPI::Error("Error finishing Tecplot file output.", CURRENT_FUNCTION);
  }
#ifdef HAVE_MPI
  SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif

#endif /* HAVE_TECIO */

  /*--- Compute and store the write time. ---*/

  stopTime = SU2_MPI::Wtime();

  usedTime = stopTime-startTime;

  fileSize = DetermineFilesize(val_filename);

  /*--- Compute and store the bandwidth ---*/

  bandwidth = fileSize/(1.0e6)/usedTime;
}


int64_t CTecplotBinaryFileWriter::GetHaloNodeNumber(unsigned long global_node_number, unsigned long last_local_node, vector<unsigned long> const &halo_node_list)
{
  auto it = lower_bound(halo_node_list.begin(), halo_node_list.end(), global_node_number);
  assert(it != halo_node_list.end());
  assert(*it == global_node_number);
  /* When C++11 is universally available, replace the following mouthful with "auto" */
  iterator_traits<vector<unsigned long>::const_iterator>::difference_type offset = distance(halo_node_list.begin(), it);
  assert(offset >= 0);
  return (int64_t)(last_local_node + offset + 1);
}
