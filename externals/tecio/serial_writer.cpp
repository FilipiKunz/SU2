#include "serial_writer.hpp"
#include "teciosrc/TECIO.h"

extern "C" int32_t SU2TecOpen(const char* file, const char* title, const char* fields, void** handle) {
  return tecFileWriterOpen(file, title, fields, FILEFORMAT_SZL, FILETYPE_FULL, FieldDataType_Float, nullptr, handle);
}

extern "C" int32_t SU2TecZone(void* handle, const char* title, int32_t type, int64_t nodes,
                                int64_t cells, const int32_t* locations, int32_t* zone) {
  return tecZoneCreateFE(handle, title, type, nodes, cells, nullptr, nullptr, locations, nullptr, 0, 0, 0, zone);
}

extern "C" int32_t SU2TecTime(void* handle, int32_t zone, double time, int32_t strand) {
  return tecZoneSetUnsteadyOptions(handle, zone, time, strand);
}

extern "C" int32_t SU2TecVar(void* handle, int32_t zone, int32_t var, int64_t count, const float* values) {
  return tecZoneVarWriteFloatValues(handle, zone, var, 0, count, values);
}

extern "C" int32_t SU2TecNodes(void* handle, int32_t zone, int64_t count, const int64_t* nodes) {
  return tecZoneNodeMapWrite64(handle, zone, 0, 1, count, nodes);
}

extern "C" int32_t SU2TecClose(void** handle) {
  return tecFileWriterClose(handle);
}
