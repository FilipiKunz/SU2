/*!
 * \file CParallelFileWriter.cpp
 * \brief Filewriter base class.
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

#include <utility>
#include <algorithm>
#include <climits>

#include "../../../include/output/filewriter/CFileWriter.hpp"

CFileWriter::CFileWriter(CParallelDataSorter *valDataSorter, string valFileExt):
  fileExt(std::move(valFileExt)),
  dataSorter(valDataSorter){

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  fileSize = 0.0;
  bandwidth = 0.0;

}

CFileWriter::CFileWriter(string valFileExt):
  fileExt(std::move(valFileExt)){

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  fileSize = 0.0;
  bandwidth = 0.0;

}

CFileWriter::~CFileWriter()= default;

bool CFileWriter::WriteMPIBinaryDataAll(const void *data, unsigned long sizeInBytes,
                                        unsigned long totalSizeInBytes, unsigned long offsetInBytes){

#ifdef HAVE_MPI

  startTime = SU2_MPI::Wtime();

  /*--- Each rank writes its own disjoint range at an explicit byte offset.
   * This avoids a collective operation for every output field. ---*/
  int ierr = MPI_SUCCESS;
  const auto* bytes = static_cast<const char*>(data);
  for (unsigned long pos = 0; pos < sizeInBytes; pos += INT_MAX) {
    const int count = static_cast<int>(std::min<unsigned long>(INT_MAX, sizeInBytes - pos));
    MPI_Status status;
    const int result = MPI_File_write_at(fhw, disp + offsetInBytes + pos,
                                         bytes + pos, count, MPI_BYTE, &status);
    int written = 0;
    if (result != MPI_SUCCESS || MPI_Get_count(&status, MPI_BYTE, &written) != MPI_SUCCESS || written != count)
      ierr = MPI_ERR_IO;
  }
  if (ierr != MPI_SUCCESS) SU2_MPI::Error("Writing binary output failed.", CURRENT_FUNCTION);

  disp      += totalSizeInBytes;
  fileSize  += sizeInBytes;

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (ierr == MPI_SUCCESS);
#else

  startTime = SU2_MPI::Wtime();

  unsigned long bytesWritten;

  /*--- Write binary data ---*/

  bytesWritten = fwrite(data, sizeof(char), sizeInBytes, fhw);
  fileSize += bytesWritten;
  if (bytesWritten != sizeInBytes) SU2_MPI::Error("Writing binary output failed.", CURRENT_FUNCTION);

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (bytesWritten == sizeInBytes);
#endif

}

bool CFileWriter::WriteMPIBinaryData(const void *data, unsigned long sizeInBytes, unsigned short processor){

#ifdef HAVE_MPI

  startTime = SU2_MPI::Wtime();

  int ierr = MPI_SUCCESS;

  if (rank == processor) {
    if (sizeInBytes > INT_MAX) SU2_MPI::Error("Binary output header exceeds MPI count limit.", CURRENT_FUNCTION);
    MPI_Status status;
    ierr = MPI_File_write_at(fhw, disp, data, int(sizeInBytes), MPI_BYTE, &status);
    int written = 0;
    if (ierr != MPI_SUCCESS || MPI_Get_count(&status, MPI_BYTE, &written) != MPI_SUCCESS || written != sizeInBytes)
      SU2_MPI::Error("Writing binary output header failed.", CURRENT_FUNCTION);
  }

  disp     += sizeInBytes;
  if (rank == processor) fileSize += sizeInBytes;

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (ierr == MPI_SUCCESS);
#else

  startTime = SU2_MPI::Wtime();

  unsigned long bytesWritten = sizeInBytes;

  /*--- Write the total size in bytes at the beginning of the binary data blob ---*/

  bytesWritten = fwrite(data, sizeof(char), sizeInBytes, fhw);
  fileSize += bytesWritten;
  if (bytesWritten != sizeInBytes) SU2_MPI::Error("Writing binary output header failed.", CURRENT_FUNCTION);

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (bytesWritten == sizeInBytes);

#endif

}

bool CFileWriter::WriteMPIString(const string &str, unsigned short processor){

#ifdef HAVE_MPI

  startTime = SU2_MPI::Wtime();

  int ierr = MPI_SUCCESS;

  if (rank == processor) {
    if (str.size() > INT_MAX) SU2_MPI::Error("Output header exceeds MPI count limit.", CURRENT_FUNCTION);
    MPI_Status status;
    ierr = MPI_File_write_at(fhw, disp, str.c_str(), int(str.size()), MPI_CHAR, &status);
    int written = 0;
    if (ierr != MPI_SUCCESS || MPI_Get_count(&status, MPI_CHAR, &written) != MPI_SUCCESS || written != str.size())
      SU2_MPI::Error("Writing output header failed.", CURRENT_FUNCTION);
  }

  disp += str.size()*sizeof(char);
  if (rank == processor) fileSize += sizeof(char)*str.size();

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (ierr == MPI_SUCCESS);

#else

  startTime = SU2_MPI::Wtime();

  unsigned long bytesWritten;
  bytesWritten = fwrite(str.c_str(), sizeof(char), str.size(), fhw);

  fileSize += bytesWritten;
  if (bytesWritten != str.size()) SU2_MPI::Error("Writing output header failed.", CURRENT_FUNCTION);

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (bytesWritten == str.size()*sizeof(char));

#endif

}

bool CFileWriter::OpenMPIFile(string val_filename){

  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);

#ifdef HAVE_MPI
  int ierr;
  disp     = 0.0;

  /*--- Truncate an existing output instead of closing an invalid handle and
   * deleting the file before reopening it. ---*/
  ierr = MPI_File_open(SU2_MPI::GetComm(), val_filename.c_str(),
                       MPI_MODE_CREATE|MPI_MODE_WRONLY,
                       MPI_INFO_NULL, &fhw);

  /*--- Error check opening the file. ---*/

  if (ierr) {
    SU2_MPI::Error(string("Unable to open file ") +
                   val_filename, CURRENT_FUNCTION);
  }
  if (MPI_File_set_size(fhw, 0) != MPI_SUCCESS) {
    SU2_MPI::Error(string("Unable to truncate file ") + val_filename, CURRENT_FUNCTION);
  }
#else
  fhw = fopen(val_filename.c_str(), "wb");
  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open file ") +
                   val_filename, CURRENT_FUNCTION);
  }
#endif

  fileSize = 0.0;
  usedTime = 0;

  return true;
}

bool CFileWriter::CloseMPIFile(){

#ifdef HAVE_MPI
  /*--- All ranks close the file after writing. ---*/

  if (MPI_File_close(&fhw) != MPI_SUCCESS)
    SU2_MPI::Error("Closing output file failed.", CURRENT_FUNCTION);
#else
  if (fclose(fhw) != 0)
    SU2_MPI::Error("Closing output file failed.", CURRENT_FUNCTION);
#endif

  /*--- Communicate the total file size for the restart ---*/

  su2double my_fileSize = fileSize;
  SU2_MPI::Allreduce(&my_fileSize, &fileSize, 1,
                     MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  /*--- Compute and store the bandwidth ---*/

  bandwidth = usedTime > 0 ? fileSize/(1.0e6)/usedTime : 0;

  return true;
}
