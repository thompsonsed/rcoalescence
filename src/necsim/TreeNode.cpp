//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Samuel Thompson
 * @file TreeNode.cpp
 * @brief Contains the TreeNode class for storing the coalescence tree.
 *
 * TreeNode objects are used both during simulation runs and afterwards, when different calculations need to be
 * performed on the coalescence tree.
 */
#include "TreeNode.h"
#include "Logging.h"

void TreeNode::setup(bool z, unsigned long xp, unsigned long yp, long xi, long yi)
{
    tip = z;
    parent = 0;
    speciated = false;

    species_id = 0;

    xpos = xp;
    ypos = yp;
    xwrap = xi;
    ywrap = yi;
    speciation_probability = 0;
    generations_existed = 0;
    generation_added = 0;
}

void TreeNode::setup(bool z)
{
    setup(z, 0, 0, 0, 0);
}

void TreeNode::setup(const bool &is_tip, const long &xp, const long &yp, const long &xi, const long &yi,
                     const long double &generation)
{
    tip = is_tip;
    parent = 0;
    speciated = false;
    species_id = 0;
    xpos = xp;
    ypos = yp;
    xwrap = xi;
    ywrap = yi;
    speciation_probability = 0;
    generations_existed = 0;
    generation_added = generation;
}

void TreeNode::setExistence(bool b)
{
    does_exist = b;
}

void TreeNode::setParent(unsigned long x)
{
    parent = x;
}

void TreeNode::qReset()
{
    species_id = 0;
    does_exist = false;
    speciated = false;
}

void TreeNode::setPosition(long x, long y, long xw, long yw)
{
    xpos = x;
    ypos = y;
    xwrap = xw;
    ywrap = yw;
}

void TreeNode::setSpec(long double d)
{
    speciation_probability = d;
}

void TreeNode::setGenerationRate(unsigned long g)
{
    generations_existed = g;
}

void TreeNode::setGeneration(long double d)
{
    generation_added = d;
}

void TreeNode::setSpeciation(bool s)
{
    speciated = s;
}

void TreeNode::burnSpecies(unsigned long idin)
{
    if(species_id == 0)
    {
        species_id = idin;
    }
}

void TreeNode::setTip(bool b)
{
    tip = b;
}

void TreeNode::resetSpecies()
{
    species_id = 0;
}

void TreeNode::increaseGen()
{
    generations_existed++;
}

bool TreeNode::exists() const
{
    return does_exist;
}

bool TreeNode::isTip() const
{
    return tip;
}

unsigned long TreeNode::getParent() const
{
    return parent;
}

unsigned long TreeNode::getXpos() const
{
    return xpos;
}

unsigned long TreeNode::getYpos() const
{
    return ypos;
}

long TreeNode::getXwrap() const
{
    return xwrap;
}

long TreeNode::getYwrap() const
{
    return ywrap;
}

bool TreeNode::hasSpeciated() const
{
    return speciated;
}

unsigned long TreeNode::getSpeciesID() const
{
    return species_id;
}

long double TreeNode::getSpecRate() const
{
    return speciation_probability;
}

unsigned long TreeNode::getGenRate() const
{
    return generations_existed;
}

long double TreeNode::getGeneration() const
{
    return generation_added;
}

void TreeNode::speciate()
{
    speciated = true;
}

ostream &operator<<(ostream &os, const TreeNode &t)
{
    os << setprecision(64);
    os << t.tip << "," << t.parent << "," << t.speciated << "," << t.does_exist << "," << t.species_id << "," << t.xpos
       << "," << t.ypos << "," << t.xwrap << ",";
    os << t.ywrap << "," << t.speciation_probability << "," << t.generations_existed << "," << t.generation_added
       << "\n";
    return os;
}

istream &operator>>(istream &is, TreeNode &t)
{
    //is << m.num_rows<<" , "<<m.num_cols<<" , "<<endl;
    char delim;
    is >> t.tip >> delim >> t.parent >> delim >> t.speciated >> delim >> t.does_exist >> delim >> t.species_id >> delim
       >> t.xpos >> delim;
    is >> t.ypos >> delim >> t.xwrap >> delim >> t.ywrap >> delim >> t.speciation_probability >> delim
       >> t.generations_existed >> delim >> t.generation_added;
    return is;
}

TreeNode &TreeNode::operator=(const TreeNode &t)
{
    tip = t.tip;
    parent = t.parent;
    speciated = t.speciated;
    does_exist = t.does_exist;
    species_id = t.species_id;
    xpos = t.xpos;
    ypos = t.ypos;
    xwrap = t.xwrap;
    ywrap = t.ywrap;
    speciation_probability = t.speciation_probability;
    generation_added = t.generation_added;
    generations_existed = t.generations_existed;
    return *this;
}

#ifdef DEBUG

void TreeNode::logLineageInformation(const int &level)
{
    writeLog(level, "Logging lineage information");
    writeLog(level, "parent: " + to_string(parent));
    writeLog(level, "tip: " + to_string(tip));
    writeLog(level, "speciated: " + to_string(speciated));
    writeLog(level, "existance: " + to_string(does_exist));
    writeLog(level, "x, y, (x wrap, y wrap): " + to_string(xpos) + ", " + to_string(ypos) + ", " +
                    to_string(xwrap) + ", " + to_string(ywrap));
    writeLog(level, "speciation rate: " + to_string(speciation_probability));
    writeLog(level, "generations (added, existed): " + to_string(generation_added) + ", " +
                    to_string(generations_existed));
}

#endif // DEBUG