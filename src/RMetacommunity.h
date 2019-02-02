////
//// Created by Sam Thompson on 28/01/2019.
////
//
//#ifndef RCOALESCENCE_RMETACOMMUNITY_H
//#define RCOALESCENCE_RMETACOMMUNITY_H
//#include <string>
//#ifndef WIN_INSTALL
//
//#include <unistd.h>
//
//#endif
//
//#ifdef CXX14_SUPPORT
//#include "memory.h"
//#else
//
//#include <memory>
//
//#endif
//#include "necsim/Metacommunity.h"
//#include "necsim/SpecSimParameters.h"
//
//using namespace std;
//
///**
// * @brief Provides the Metacommunity object for R.
// */
//class RMetacommunity : public virtual Metacommunity
//{
//public:
//    shared_ptr<SpecSimParameters> spec_sim_parameters;
//
//    /**
//     * @brief Default constructor for the RMetacommunity.
//     */
//    RMetacommunity();
//
//    /**
//     * @brief Sets up the metacommunity object with the base parameters.
//     * @param database the database to read in from
//     * @param record_spatial if true, record the full spatial data of the simulation
//     * @param sample_file a sample file to provide sampling routines from
//     * @param use_fragments_in either a path to a fragments file, or True/False
//     */
//    void setup(const string &database, const bool &record_spatial, const string &sample_file,
//               const string &use_fragments_in);
//
//    /**
//     * @brief Adds a speciation rate for the community to be calculated for.
//     * @param spec_rate the speciation rate for determining the community
//     */
//    void addSpeciationRate(const double &spec_rate);
//
//    /**
//     * @brief Adds a time to determine the community at.
//     * @param time the time in generations for community reconstruction
//     */
//    void addTime(const double &time);
//
//    /**
//     * @brief Adds metacommunity parameters to calculate the community for.
//     * @param metacommunity_size_in the number of individuals in the metacommunity
//     * @param metacommunity_speciation_rate_in the speciation rate of the metacommunity
//     * @param metacommunity_option_in the metacommunity option
//     * @param metacommunity_reference_in the external metacommunity reference
//     */
//    void addMetacommunityParameters(const unsigned long &metacommunity_size_in,
//                                    const double &metacommunity_speciation_rate_in,
//                                    const string &metacommunity_option_in,
//                                    const unsigned long &metacommunity_reference_in);
//
//    /**
//     * @brief Adds protracted speciation parameters to the metacommunity.
//     * @param proc_spec_min the protracted speciation minimum generation
//     * @param proc_spec_max the protracted speciation maximum generation
//     */
//    void addProtractedParameters(const double &proc_spec_min, const double &proc_spec_max);
//
//    /**
//     * @brief Applies the speciation parameters to the curently community and outputs the results to the database.
//     */
//    void apply();
//
//};
//
//#endif //RCOALESCENCE_RMETACOMMUNITY_H
