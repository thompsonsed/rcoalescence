// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file Metacommunity.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Contains the Metacommunity class for generating a neutral metacommunity.
 *
 * For use with completed simulations from necsim, using the SpeciationCounter program. Individuals will be drawn from
 * the metacommunity for each speciation event, instead of creating a new species each time.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#include <utility>
#include <unordered_map>
#include "Metacommunity.h"
#include "neutral_analytical.h"
#include "LogFile.h"
#include "SpeciesAbundancesHandler.h"
#include "SimulatedSpeciesAbundancesHandler.h"
#include "AnalyticalSpeciesAbundancesHandler.h"

namespace na = neutral_analytical;

Metacommunity::Metacommunity() : seed(0), job_type(0), parameters_checked(false),
                                 species_abundances_handler(make_unique<SimulatedSpeciesAbundancesHandler>()),
                                 random(make_shared<RNGController>()),
                                 metacommunity_tree(make_unique<Tree>())
{

}

void Metacommunity::setCommunityParameters(shared_ptr<MetacommunityParameters> metacommunity_parameters)
{
    current_metacommunity_parameters = std::move(metacommunity_parameters);
    // Check if no metacommunity option has been supplied - in which case use sensible defaults
    if(current_metacommunity_parameters->option == "none" || current_metacommunity_parameters->option == "null")
    {
        // Simulate only for smaller community sizes
        if(current_metacommunity_parameters->metacommunity_size > 1000)
        {
            current_metacommunity_parameters->option = "analytical";
        }
        else
        {
            current_metacommunity_parameters->option = "simulated";
        }
    }
}

void Metacommunity::checkSimulationParameters()
{
    if(!parameters_checked)
    {
        if(!database->isOpen())
        {
            throw FatalException(
                    "Cannot read simulation metacommunity parameters as database is null pointer.");
        }
        // Now do the same for times
        string sql_call = "SELECT seed, job_type from SIMULATION_PARAMETERS";
        auto stmt = database->prepare(sql_call);
        database->step();
        seed = static_cast<unsigned long>(sqlite3_column_int(stmt->stmt, 0));
        random->setSeed(seed);
        job_type = static_cast<unsigned long>(sqlite3_column_int(stmt->stmt, 1));
        database->finalise();
        parameters_checked = true;
    }
}

void Metacommunity::addSpecies(unsigned long &species_count, TreeNode *tree_node, set<unsigned long> &species_list)
{

    auto species_id = species_abundances_handler->getRandomSpeciesID();
    if(species_id == 0)
    {
        throw FatalException("Obtained species id was 0 in metacommunity application - please report this bug.");
    }

    if(species_list.empty() || species_list.find(species_id) == species_list.end())
    {
        species_list.insert(species_id);
        species_count++;
    }
#ifdef DEBUG
    if(tree_node->getSpeciesID() != 0)
    {
        throw FatalException("Trying to add species for lineages with non-zero species id. Please report this bug.");
    }
#endif // DEBUG
    tree_node->burnSpecies(species_id);
}

void Metacommunity::createMetacommunityNSENeutralModel()
{
#ifdef DEBUG
    writeLog(10, "Running spatially implicit model for metacommunity generation.");
#endif //DEBUG
    // First set up a non-spatial coalescence simulation to generate our metacommunity
    shared_ptr<SimParameters> temp_parameters = make_shared<SimParameters>();
    // Generate a new unique seed by adding 1073741823 if the seed is 0 - this ensures that 0 and 1 never appear as
    // random seeds, which cause autocorrelation in simulation outputs.
    if(seed == 0)
    {
        seed = 1073741823;
    }
    temp_parameters->setMetacommunityParameters(current_metacommunity_parameters->metacommunity_size,
                                                current_metacommunity_parameters->speciation_rate, seed,
                                                job_type);
    // Dispose of any previous Tree object and create a new one
    metacommunity_tree = make_unique<Tree>();
    metacommunity_tree->internalSetup(temp_parameters);
    // Run our simulation and calculate the species abundance distribution (as this is all that needs to be stored).
    if(!metacommunity_tree->runSimulation())
    {
        throw FatalException("Completion of the non-spatial coalescence simulation "
                             "to create the metacommunity did not finish in time.");
    }
    metacommunity_tree->applySpecRateInternal(current_metacommunity_parameters->speciation_rate, 0.0);
    // species_abundances now contains the number of individuals per species
    // Make it cumulative to increase the speed of indexing using binary search.
    species_abundances_handler = make_unique<SimulatedSpeciesAbundancesHandler>();
    species_abundances_handler->setup(random, current_metacommunity_parameters->metacommunity_size,
                                      current_metacommunity_parameters->speciation_rate, nodes->size());
    auto tmp_species_abundances = metacommunity_tree->getSpeciesAbundances();
    if(tmp_species_abundances->empty())
    {
        throw FatalException("Simulated species abundance list is empty. Please report this bug.");
    }
    // Remove the 0 at the start
    species_abundances_handler->setAbundanceList(tmp_species_abundances);
#ifdef DEBUG
    writeLog(10, "Spatially implicit simulation completed.");
#endif //DEBUG
}

void Metacommunity::applyNoOutput(shared_ptr<SpecSimParameters> sp, shared_ptr<vector<TreeNode>> tree_data)
{
#ifdef DEBUG
    writeLog(10, "********************");
    writeLog(10, "Metacommunity application");
#endif //DEBUG

    // Make sure that the connection is opened to file.
    if(!bSqlConnection)
    {
        openSqlConnection(sp->filename);
    }
    checkSimulationParameters();
    for(const auto &item: sp->metacommunity_parameters)
    {
        setCommunityParameters(item);
        printMetacommunityParameters();
        setupApplication(sp, tree_data);
        if(current_metacommunity_parameters->option == "simulated")
        {
            createMetacommunityNSENeutralModel();
        }
        else if(current_metacommunity_parameters->option == "analytical")
        {
            // Use approximation for the SAD from Chisholm and Pacala (2010)
            approximateSAD();
        }
        else
        {
            // Use the file path provided
            readSAD();
        }
#ifdef DEBUG
        writeLog(10, "Creating coalescence tree from metacommunity...");
#endif //DEBUG
        calculateTree();
    }
}

void Metacommunity::approximateSAD()
{
    species_abundances_handler.reset();
    species_abundances_handler = make_unique<AnalyticalSpeciesAbundancesHandler>();
    species_abundances_handler->setup(random, current_metacommunity_parameters->metacommunity_size,
                                      current_metacommunity_parameters->speciation_rate, nodes->size());
}

void Metacommunity::readSAD()
{
    Community external_metacommunity;
    external_metacommunity.openSqlConnection(current_metacommunity_parameters->option);
    shared_ptr<map<unsigned long, unsigned long>> sad = external_metacommunity.getSpeciesAbundances(
            current_metacommunity_parameters->external_reference);
    species_abundances_handler.reset();
    species_abundances_handler = make_unique<SimulatedSpeciesAbundancesHandler>();
    species_abundances_handler->setup(random, current_metacommunity_parameters->metacommunity_size,
                                      current_metacommunity_parameters->speciation_rate, nodes->size());
    species_abundances_handler->setAbundanceList(sad);
}

void Metacommunity::printMetacommunityParameters()
{
    stringstream ss;
    ss << "Metacommunity current_metacommunity_parameters:" << endl;
    ss << "Metacommunity size: " << current_metacommunity_parameters->metacommunity_size << endl;
    ss << "Speciation rate: " << current_metacommunity_parameters->speciation_rate << endl;
    ss << "Option: " << current_metacommunity_parameters->option << endl;
    ss << "External reference: " << current_metacommunity_parameters->external_reference << endl;
    writeInfo(ss.str());
}






