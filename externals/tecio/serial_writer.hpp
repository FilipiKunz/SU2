#pragma once

#include <cstdint>

extern "C" {
int32_t SU2TecOpen(const char* file, const char* title, const char* fields, void** handle);
int32_t SU2TecZone(void* handle, const char* title, int32_t type, int64_t nodes,
                   int64_t cells, const int32_t* locations, int32_t* zone);
int32_t SU2TecTime(void* handle, int32_t zone, double time, int32_t strand);
int32_t SU2TecVar(void* handle, int32_t zone, int32_t var, int64_t count, const float* values);
int32_t SU2TecNodes(void* handle, int32_t zone, int64_t count, const int64_t* nodes);
int32_t SU2TecClose(void** handle);
}
