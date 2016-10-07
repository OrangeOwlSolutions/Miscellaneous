#ifndef __PACKED_DATA_TYPES_h__
#define __PACKED_DATA_TYPES_h__

#include <stdint.h>

namespace nv_types{

	#include "packed.h"
    //struct __attribute__((__packed__)) nv_int2 {
    struct nv_int2 {
        int32_t x, y;
    //};
    } PACKED;
	#include "endpacked.h"

	#include "packed.h"
    //struct __attribute__((__packed__)) nv_int4 {
    struct nv_int4 {
        int32_t x, y, z, w;
    //};
    } PACKED;
	#include "endpacked.h"

	#include "packed.h"
    //struct __attribute__((__packed__)) nv_float3{
    struct nv_float3 {
        float x, y, z;
    //};
    } PACKED;
	#include "endpacked.h"

	#include "packed.h"
    //struct __attribute__((__packed__)) nv_float4{
    struct nv_float4 {
        float x, y, z, w;
    //};
    } PACKED;
	#include "endpacked.h"
    
}

#endif
