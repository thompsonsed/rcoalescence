// This file is part of rcoalescence project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file RProtractedSpatialTreeSimulation.h
 * @brief Wraps the ProtractedSpatialTree class for running protracted spatially-explicit neutral models from R.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef RCOALESCENCE_RPROTRACTEDSPATIALTREE_H
#define RCOALESCENCE_RPROTRACTEDSPATIALTREE_H


#ifndef WIN_INSTALL
#include <unistd.h>
#endif // WIN_INSTALL

#include <Rcpp.h>
#include "necsim/ProtractedSpatialTree.h"
#include "RGenericSpatialTreeSimulation.h"
using namespace necsim;

namespace rcoalescence
{
    /**
     * @brief Wrapper for spatial simulations.
     */
    using RProtractedSpatialTreeSimulation = RGenericSpatialTreeSimulation<necsim::ProtractedSpatialTree>;
}
#endif //RCOALESCENCE_RPROTRACTEDSPATIALTREE_H
