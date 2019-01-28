//
// Created by Sam Thompson on 28/01/2019.
//

#ifndef RCOALESCENCE_RMETACOMMUNITY_H
#define RCOALESCENCE_RMETACOMMUNITY_H
#include <string>

#include "Metacommunity.h"
#include "SpecSimParameters.h"

using namespace std;

class RMetacommunity : Metacommunity
{
public:
    SpecSimParameters spec_sim_parameters;

    RMetacommunity() : spec_sim_parameters()
    {

    }

    /**
     * @brief Sets up the metacommunity object with the base parameters.
     * @param database the database to read in from
     * @param record_spatial if true, record the full spatial data of the simulation
     * @param sample_file a sample file to provide sampling routines from
     * @param use_fragments_in either a path to a fragments file, or True/False.
     */
    void setup(const string &database, const bool &record_spatial, const string &sample_file,
               const string &use_fragments_in)
    {
        spec_sim_parameters.setup(database, record_spatial, sample_file, use_fragments_in);
    }

    void addSpeciationRate(const double &spec_rate)
    {
        spec_sim_parameters.addSpeciationRate(spec_rate);
    }

    void addTime(const double &time)
    {
        spec_sim_parameters.addTime(time);
    }

    void addMetacommunityParameters(const unsigned long &metacommunity_size_in,
                                    const double &metacommunity_speciation_rate_in,
                                    const string &metacommunity_option_in,
                                    const unsigned long &metacommunity_reference_in)

    {
        spec_sim_parameters.addMetacommunityParameters(metacommunity_size_in, metacommunity_speciation_rate_in,
                                                       metacommunity_option_in, metacommunity_reference_in);
    }

    void apply()
    {

    }

};

#endif //RCOALESCENCE_RMETACOMMUNITY_H
