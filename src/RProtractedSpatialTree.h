//
// Created by Sam Thompson on 18/05/2018.
//

#ifndef RCOALESCENCE_RPROTRACTEDSPATIALTREE_H
#define RCOALESCENCE_RPROTRACTEDSPATIALTREE_H

#include <necsim/ProtractedTree.h>
#include <necsim/ProtractedSpatialTree.h>
#include "RSpatialTree.h"

class RProtractedSpatialTree : virtual public RSpatialTree, virtual public ProtractedSpatialTree
{

};

#endif //RCOALESCENCE_RPROTRACTEDSPATIALTREE_H
