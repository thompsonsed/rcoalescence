// This file is part of rcoalescence project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file RGenericSpatialTreeSimulation.h
 * @brief Wraps the SpatialTree class for running spatially-explicit neutral models from R.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#ifndef RCOALESCENCE_R_SPATIAL_TREE_H
#define RCOALESCENCE_R_SPATIAL_TREE_H

#ifndef WIN_INSTALL
#include <unistd.h>
#endif // WIN_INSTALL

#include <Rcpp.h>
#include "necsim/SpatialTree.h"
#include "RGenericSpatialTreeSimulation.h"
using namespace necsim;

namespace rcoalescence
{
/**
 * @brief Wrapper for spatial simulations.
 */
  using RSpatialTreeSimulation = RGenericSpatialTreeSimulation<necsim::SpatialTree>;
}

#endif // RCOALESCENCE_R_SPATIAL_TREE_H