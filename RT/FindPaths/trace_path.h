#ifndef __TRACE_PATH_H__
#define __TRACE_PATH_H__

#include <vector>

#include "ray_path.h"
#include "load_nbin.h"
#include "RXTX.h"

void trace_path(
    std::vector<path> &global_path_list, 
    const  nvmesh  &mesh,
    const std::vector<TX> &txs,
    const std::vector<RX> &rxs
);

#endif
