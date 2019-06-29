//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Samuel Thompson
 * @file TreeNode.h
 * @brief Contains the TreeNode class for storing the coalescence tree.
 * 
 * TreeNode objects are used both during simulation runs and afterwards, when different calculations need to be performed on the coalescence tree.
 */
#ifndef TREENODE
#define TREENODE

#include <cstdio>
#include <iostream>
#include <iomanip>
#include "Logging.h"

using namespace std;

/**
 * @class TreeNode
 * @brief The TreeNode class that acts as a data storage object for the phylogenetic tree.
 * 
 * Also contains all the necessary routines for changes to a lineage's attributes, called by TreeList objects when
 * generating new coalescence trees.
 */
class TreeNode
{

private:
    // 0 means that this node is just here to mark a coalescense
    // and therefore this node of no real other relevance
    // 1 means that this node is a leaf node and counts towards diversity
    bool tip;
    // this stores the parent of the individual
    // 0 means there is no parent - we are at the end of the tree
    // (as far as has been calculated)
    unsigned long parent;
    // true if this lineage has speciated in which case it should not have a parent
    // because under the present implementation lineages are not traced beyond speciation
    // boolean for checking whether the lineage actually exists at the end. If all children of the lineages have
    // speciated, then the lineage no longer exists.
    bool speciated;
    // the species identity of the node
    bool does_exist;
    // the following 4 variables describe the position of the lineage in the present day.
    unsigned long species_id;
    // x position
    unsigned long xpos;
    // y position
    unsigned long ypos;
    // number of wraps of x around the torus
    long xwrap;
    // number of wraps of y around the torus
    long ywrap;
    // the speciation probability. This needs to be multiplied by the number of generations in order to generate the
    // actual probability.
    long double speciation_probability;
    // Number of generations this lineage has existed since coalescence (or tracking began).
    unsigned long generations_existed;
    // Simulation generation timer that the lineage was created at
    long double generation_added;
public:
    /**
     * @brief The default constructor.
     */
    TreeNode() : tip(false), parent(0), speciated(false), does_exist(false), species_id(0), xpos(0), ypos(0), xwrap(0),
                 ywrap(0), speciation_probability(0.0), generations_existed(0), generation_added(0.0)
    {

    }

    /**
     * @brief The default destructor.
     */
    ~TreeNode()
    = default;

    /**
     * @brief Sets up variables with initial conditions.
     * @param z whether this lineage is a tip or not (represents the end of a tree).
     * @param xp the x position on the grid.
     * @param yp the y position on the grid.
     * @param xi the number of wraps in the x dimension.
     * @param yi the number of wraps in the y dimension.
     */
    void setup(bool z, unsigned long xp, unsigned long yp, long xi, long yi);

    /**
     * @brief Sets up variables with initial conditions.
     * @param z whether this lineage is a tip or not (represents the end of a tree)
     */
    void setup(bool z);

    /**
     * @brief Overloaded setup() function, additionally taking a generation time point.
     *
     * Used when creating lineages after the start of the simulation, when the generation is not 0.
     * @param is_tip whether this lineage is a tip or not (represents the end of a tree).
     * @param xp the x position on the grid.
     * @param yp the y position on the grid.
     * @param xi the number of wraps in the x dimension.
     * @param yi the number of wraps in the y dimension.
     * @param generation the current generation at creation time.
     */
    void setup(const bool &is_tip, const long &xp, const long &yp, const long &xi, const long &yi,
               const long double &generation);

    /**
     * @brief Setter for the existence of the lineage.
     * @param b existence boolean.
     */
    void setExistence(bool b);

    /**
     * @brief Setter for the parent reference.
     * @param x a reference for the parent location.
     */
    void setParent(unsigned long x);

    /**
     * @brief Resets the lineage.
     * Remove any species ID, existence and speciation record.
     */
    void qReset();

    /**
     * @brief Set a new position for the lineage.
     * @param x the x position on the grid.
     * @param y the y position on the grid.
     * @param xw the number of wraps in the x dimension.
     * @param yw the number of wraps in the y dimension.
     */
    void setPosition(long x, long y, long xw, long yw);

    /**
     * @brief Setter for the randomly generated number (from NRrand.d0()) for the speciation probability for this lineage.
     * @param d the speciation probability.
     */
    void setSpec(long double d);

    /**
     * @brief Setter for the number of generations this lineage has existed for.
     * @param g the number of generations that the lineage has existed.
     */
    void setGenerationRate(unsigned long g);

    /**
     * @brief Setter for the birth generation timer for the lineage.
     *
     * Note that moves that don't involve coalescence do not create a new Treenode object,
     * and therefore the generation_added does not get updated.
     * However, coalescence events will cause a new Treenode object creation.
     * The lineage birth generation is generally only important for calculating the age of tips.
     * @param d
     */
    void setGeneration(long double d);

    /**
     * @brief Setter for the speciation boolean.
     * @param s the speciation boolean.
     */
    void setSpeciation(bool s);

    /**
     * @brief Setter for the species ID.
     * Once set to something other than 0, this cannot be changed with a call to qReset() or resetSpecies().
     * @param idin the species ID.
     */
    void burnSpecies(unsigned long idin);

    /**
     * @brief Setter for the tip boolean.
     * @param b the tip boolean.
     */
    void setTip(bool b);

    /**
     * @brief Reset the species ID to 0.
     */
    void resetSpecies();

    /**
     * @brief Increases the generation counter by one.
     */
    void increaseGen();
    // we don't allow the other variables to be changed
    // because they only need to be set once at the start of the coalescence
    // it's actually safer to leave out setters.
    // similarly we don't allow speciation to be changed once it has been set.

    // standard getters

    /**
     * @brief Getter for the existence boolean.
     * @return the existence boolean.
     */
    bool exists() const;

    /**
     * @brief Getter for the tip boolean.
     * @return the tip boolean.
     */
    bool isTip() const;

    /**
     * @brief Getter for the parent location.
     * @return the parent reference.
     */
    unsigned long getParent() const;

    /**
     * @brief Getter for the x position.
     * @return the x position.
     */
    unsigned long getXpos() const;

    /**
     * @brief Getter for the y position.
     * @return the y position.
     */
    unsigned long getYpos() const;

    /**
     * @brief Getter for the number of times the lineage is wrapped in the x dimension.
     * @return the number of times the lineage is wrapped in the x dimension.
     */
    long getXwrap() const;

    /**
     * @brief Getter for the number of times the lineage is wrapped in the y dimension.
     * @return the number of times the lineage is wrapped in the y dimension.
     */
    long getYwrap() const;

    /**
     * @brief Getter for the speciation boolean.
     * @return the speciation boolean.
     */
    bool hasSpeciated() const;

    /**
     * @brief Getter for the species ID.
     * @return the species ID.
     */
    unsigned long getSpeciesID() const;

    /**
     * @brief Getter for the randomly generated speciation probability.
     * @return the speciation probability.
     */
    long double getSpecRate() const;

    /**
     * @brief Getter for the number of generations the lineage has existed.
     * @return the number of generations of existence.
     */
    unsigned long getGenRate() const;

    /**
     * @brief Getter for the generation the lineage was created.
     * @return the generation counter the lineage was created.
     */
    long double getGeneration() const;

    /**
     * @brief Sets the speciation boolean to true.
     */
    void speciate();

    /**
     * @brief Overloading the << operator for outputting a Treenode object to an output stream.
     * @param os the output stream.
     * @param t a Treenode object to output.
     * @return the output stream.
     */
    friend ostream &operator<<(ostream &os, const TreeNode &t);

    /**
     * @brief Overloading the >> operator for inputting the Treenode object from an input stream.
     * @param is the input stream.
     * @param t a Treenode object to input to.
     * @return the input stream.
     */
    friend istream &operator>>(istream &is, TreeNode &t);

    /**
     * Overloading the equality operator for TreeNodes
     * @param t the input TreeNode
     */
    TreeNode &operator=(const TreeNode &t);

#ifdef DEBUG

    /**
     * @brief Logs lineage information at the level provided.
     *
     * Logs all lineage variables on separate logging lines.
     * @param level the level to be logged at
     */
    void logLineageInformation(const int &level);

#endif // DEBUG
};

#endif