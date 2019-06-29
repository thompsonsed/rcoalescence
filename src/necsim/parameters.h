//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @file parameters.h
 *
 * @brief Contains various objects used for storing parameters related to tree reconstruction.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef NECSIM_PARAMETERS_H
#define NECSIM_PARAMETERS_H

#include <cstring>
#include <stdexcept>
#include <string>
#include <memory>
#include <vector>
#include "double_comparison.h"

using namespace std;

struct ProtractedSpeciationParameters
{
    double min_speciation_gen;
    double max_speciation_gen;

    ProtractedSpeciationParameters() : min_speciation_gen(0), max_speciation_gen(0){};

    bool operator==(const ProtractedSpeciationParameters &p1) const
    {
        return (doubleCompare(p1.min_speciation_gen, min_speciation_gen, 0.00000001) &&
                doubleCompare(p1.max_speciation_gen, max_speciation_gen, 0.00000001));
    }

};

/**
 * @struct CommunityParameters
 * @brief A struct for containing pairs of previous calculations to make sure that aren't repeated.
 */
struct CommunityParameters
{
    unsigned long reference;
    long double speciation_rate;
    long double time;
    bool fragment;
    unsigned long metacommunity_reference; // will be 0 if no metacommunity used.
    // protracted speciation current_metacommunity_parameters
    /**
     * @brief The protracted speciation parameters for this object
     */
    ProtractedSpeciationParameters protracted_parameters;
    bool updated; // set to true if the fragment reference needs updating in the database

    /**
     * @brief Default constructor.
     */
    CommunityParameters() : reference(0), speciation_rate(0.0), time(0), fragment(false), metacommunity_reference(0),
                            protracted_parameters(), updated(false){}

    /**
     * @brief Trivial destructor.
     */
    ~CommunityParameters() = default;

    /**
     * @brief Constructor for CommunityParameters, for storing a pairs of previous calculations, requiring a speciation rate
     * and a time.
     * Overloaded version with setup routines.
     * @param reference_in the reference to set for this CommunityParameters set
     * @param speciation_rate_in the speciation rate of the previous calculation
     * @param time_in the time of the previous calculation
     * @param fragment_in bool of whether fragments  were used in the previous calculation
     * @param metacommunity_reference_in the metacommunity reference, or 0 for no metacommunity
     * @param protracted_params protracted speciation parameters to add
     */
    CommunityParameters(unsigned long reference_in, long double speciation_rate_in, long double time_in,
                        bool fragment_in, unsigned long metacommunity_reference_in,
                        const ProtractedSpeciationParameters &protracted_params);

    /**
     * @brief Sets up the CommunityParameters object
     * @param reference_in the reference to set for this CommunityParameters set
     * @param speciation_rate_in the speciation rate of the previous calculation
     * @param time_in the time of the previous calculation
     * @param fragment_in bool of whether fragments  were used in the previous calculation
     * @param metacommunity_reference_in the metacommunity reference, or 0 for no metacommunity
     * @param protracted_params protracted speciation parameters to add
     */
    void setup(unsigned long reference_in, long double speciation_rate_in, long double time_in, bool fragment_in,
               unsigned long metacommunity_reference_in, const ProtractedSpeciationParameters &protracted_params);

    /**
     * @brief Compare these set of parameters with the input set. If they match, return true, otherwise return false
     * @param speciation_rate_in speciation rate to compare with stored community parameter
     * @param time_in time to compare with stored community parameter
     * @param fragment_in if fragments are being used on this database
     * @param metacommunity_reference_in metacommunity reference to compare with stored community parameter
     * @param protracted_params the minimum number of generations required for existance before speciation
     * @param protracted_params protracted speciation parameters to add
     * @return
     */
    bool compare(long double speciation_rate_in, long double time_in, bool fragment_in,
                 unsigned long metacommunity_reference_in,
                 const ProtractedSpeciationParameters &protracted_params);

    /**
     * @brief Compare these set of parameters with the input set. If they match, return true, otherwise return false
     * Overloaded version ignoring the fragments parameter
     * @param speciation_rate_in speciation rate to compare with stored community parameter
     * @param time_in time to compare with stored community parameter
     * @param metacommunity_reference_in metacommunity reference to compare with stored community parameter
     * @param protracted_params protracted speciation parameters to add
     * @return
     */
    bool compare(long double speciation_rate_in, long double time_in,
                 unsigned long metacommunity_reference_in,
                 const ProtractedSpeciationParameters &protracted_params);

    /**
     * @brief Checks if the supplied reference is the same in the community parameter
     * @param reference_in
     * @return
     */
    bool compare(unsigned long reference_in);
};

/**
 * @brief A structure for containing an array of previous calculation information, including which fragments have been
 * already calculated for.
 */
struct CommunitiesArray
{
    /**
     * @brief The array of CommunityParameters which have been stored.
     */
    vector<shared_ptr<CommunityParameters>> comm_parameters;

    /**
     * @brief Default constructor.
     */
    CommunitiesArray();

    /**
     * @brief Trivial destructor
     */
    ~CommunitiesArray() = default;

    /**
     * @brief Adds an extra CommunityParameters object to the calc_array vector with the supplied variables
     * @param reference the reference for this set of community parameters
     * @param speciation_rate the speciation rate of the past calculation
     * @param time the time of the past calculation
     * @param fragment bool of whether fragments were used in the past calculation
     * @param metacommunity_reference reference for the metacommunity parameters, or 0 if no metacommunity
     * @param protracted_params protracted speciation parameters to add
     */
    void pushBack(unsigned long reference, long double speciation_rate, long double time, bool fragment,
                  unsigned long metacommunity_reference,
                  const ProtractedSpeciationParameters &protracted_params);

    /**
     * @brief Adds the provided CommunityParameters object to the calc_array vector
     * @param tmp_param the set of community parameters to add
     */
    void pushBack(shared_ptr<CommunityParameters> tmp_param);

    /**
     * @brief Adds a new communities calculation paremeters reference, with a new unique reference
     * @param speciation_rate the speciation rate of the new calculation
     * @param time the time used in the new calculation
     * @param fragment true if fragments were used in the new calculation
     * @param metacommunity_reference the reference to the set of metacommunity parameters (0 for none)
     * @param protracted_params protracted speciation parameters to add
     * @return reference to the new CommunityParameters object added
     */
    shared_ptr<CommunityParameters> addNew(long double speciation_rate, long double time, bool fragment,
                                           unsigned long metacommunity_reference,
                                           const ProtractedSpeciationParameters &protracted_params);

    /**
     * @brief Checks whether the calculation with the supplied variables has already been performed.
     *
     * @note
     * @param speciation_rate the speciation rate to check for
     * @param time the time to check for
     * @param fragment bool for checking if fragments were used
     * @param metacommunity_reference the reference to the set of metacommunity parameters (0 for none)
     * @param protracted_params protracted speciation parameters to add
     * @return true if the reference exists in past community parameters
     */
    bool hasPair(long double speciation_rate, double time, bool fragment,
                 unsigned long metacommunity_reference,
                 const ProtractedSpeciationParameters &protracted_params);

};

/**
 * @brief Contains a set of metacommunity parameters that have been applied, or are to be applied, to the coalescence
 * tree.
 */
struct MetacommunityParameters
{
    unsigned long reference{};
    unsigned long metacommunity_size{};
    long double speciation_rate{};
    string option;
    unsigned long external_reference{};

    /**
     * @brief Default constructor.
     */
    MetacommunityParameters();

    /**
     * @brief Trivial destructor
     */
    ~MetacommunityParameters() = default;

    /**
     * @brief Constructor for MetacommunityParameters, storing a previously applied metacommunity
     * @param reference_in the metacommunity reference number
     * @param metacommunity_size_in the speciation rate used for metacommunity generation
     * @param speciation_rate_in size of the tested metacommunity
     * @param option the metacommunity option ("simulated", "analytical" or a path to a database)
     * @param external_reference the reference for the external database
     */
    MetacommunityParameters(const unsigned long &reference_in, const unsigned long &metacommunity_size_in,
                            const long double &speciation_rate_in, const string &option_in,
                            const unsigned long &external_reference_in);

    /**
     * @brief Compare these set of parameters with the input set. If they match, return true, otherwise return false
     * @param speciation_rate_in speciation rate to compare with stored community parameter
     * @param metacommunity_size_in size of the tested metacommunity
     * @param option_in the metacommunity option ("simulated", "analytical" or a path to a database)
     * @param external_reference_in the reference for the external database
     * @return true if the two parameter sets are identical
     */
    bool compare(unsigned long metacommunity_size_in, long double speciation_rate_in, const string &option_in,
                 const unsigned long &ext_reference_in);

    /**
     * @brief Compare these set of parameters with the input set. If they match, return true, otherwise return false
     * @param metacomm_in
     * @return true if the two parameter sets are identical
     */
    bool compare(const MetacommunityParameters &metacomm_in);

    /**
     * @brief Checks if the supplied reference is the same in the metacommunity reference
     * @param reference_in the reference to check against
     * @return
     */
    bool compare(unsigned long reference_in);

    /**
     * @brief Checks if this combination of metacommunity parameters is not a null option
     * @return true if the combination is valid (i.e. not null)
     */
    bool isMetacommunityOption() const;

    /**
     * @brief Wipes the metacommunity parameters.
     */
    void clear();

    MetacommunityParameters &operator=(const MetacommunityParameters &parameters);

    bool operator==(const MetacommunityParameters &parameters);

    bool operator!=(const MetacommunityParameters &parameters);
};

/**
 * @brief Contains an array of MetacommunityParameters that have been applied to the coalescence tree.
 */
struct MetacommunitiesArray
{
    vector<shared_ptr<MetacommunityParameters>> metacomm_parameters;

    /**
     * @brief Default constructor.
     */
    MetacommunitiesArray();

    /**
     * @brief Trivial destructor.
     */
    ~MetacommunitiesArray() = default;

    /**
     * @brief Support for range-based for loop iteration over metacomm_parameters.
     * @return the start of the parameters vector
     */
    vector<shared_ptr<MetacommunityParameters>>::iterator begin();

    /**
     * @brief Support for range-based for loop iteration over metacomm_parameters.
     * @return the end of the parameters vector
     */
    vector<shared_ptr<MetacommunityParameters>>::iterator end();

    /**
     * @brief Support for range-based for loop iteration over metacomm_parameters.
     * @return the start of the parameters vector
     */
    vector<shared_ptr<MetacommunityParameters>>::const_iterator begin() const;

    /**
     * @brief Support for range-based for loop iteration over metacomm_parameters.
     * @return the end of the parameters vector
     */
    vector<shared_ptr<MetacommunityParameters>>::const_iterator end() const;

    /**
     * @brief Adds an extra CommunityParameters object to the calc_array vector with the supplied variables
     * @param reference the reference for this set of metacommunity parameters
     * @param speciation_rate the speciation rate used in generation of the metacommunity
     * @param metacommunity_size the size of the metacommunity used
     * @param option the metacommunity option ("simulated", "analytical" or a path to a database)
     * @param external_reference the reference for the external database
     */
    void pushBack(const unsigned long &reference, const unsigned long &metacommunity_size,
                  const long double &speciation_rate, const string &option,
                  const unsigned long &external_reference);

    /**
     * @brief Adds the provided PastMetacommunityParameters object to the calc_array vector
     * @param tmp_param the set of metacommunity parameters to add
     */
    void pushBack(shared_ptr<MetacommunityParameters> tmp_param);

    /**
     * @brief Wipes the internal vector.
     */
    void clear();

    /**
     * @brief Gets the size of the array
     * @return the number of metacommunity parameters that have been added
     */
    unsigned long size();

    /**
     * @brief Check if the array is empty
     * @return true if the array has no elements in
     */
    bool empty();

    /**
     * @brief Adds a new metacommunities calculation paremeters reference, with a new unique reference
     * @param speciation_rate the speciation rate of the new calculation
     * @param metacommunity_size the size of the metacommunity in the new calculation
     * @param option the metacommunity option ("simulated", "analytical" or a path to a database)
     * @param external_reference the reference for the external database
     * @return the new reference number, which should be unique
     */
    unsigned long addNew(const unsigned long &metacommunity_size, const long double &speciation_rate,
                         const string &option, const unsigned long &external_reference);

    /**
     * @brief Adds a new set of metacommunity parameters to the array.
     * @param metacomm_in the new metacommunity parameters to add
     * @return the new reference number, which should be unique
     */
    unsigned long addNew(const MetacommunityParameters &metacomm_in);

    /**
     * @brief Checks whether the calculation with the supplied variables has already been performed.
     * @param speciation_rate the speciation rate to check for
     * @param metacommunity_size the size of metacommunity to check for
     * @param option the metacommunity option ("simulated", "analytical" or a path to a database)
     * @param external_reference the reference for the external database
     * @return true if the reference exists in past metacommunity parameters
     */
    bool hasOption(const unsigned long &metacommunity_size, const long double &speciation_rate,
                   const string &option, const unsigned long &external_reference);

    /**
     * @brief Checks whether the calculation with the supplied reference has already been performed.
     * Overloaded version for checking references.
     * @param reference the reference to check for in past metacommunity parameters
     * @return true if the reference exists in past metacommunity parameters
     */
    bool hasOption(unsigned long reference);

    /**
     * @brief Checks whether the calculation with the supplied reference has already been performed.
     * @param metacomm_in the metacommunity parameters to check for
     * @return true if the metacommunity already exists in past metacommunity parameters
     */
    bool hasOption(const MetacommunityParameters &metacomm_in);

    /**
     * @brief Gets the metacommunity reference for the provided parameters, or returns 0 if it doesn't exist
     * @param speciation_rate the metacommunity speciation rate to obtain for
     * @param metacommunity_size the metacommunity size to apply for
     * @param fragment bool for checking if fragments were used
     * @param option the metacommunity option ("simulated", "analytical" or a path to a database)
     * @param external_reference the reference for the external database
     * @return the metacommunity reference number, or 0 if it doesn't exist
     */
    unsigned long getReference(const unsigned long &metacommunity_size, const long double &speciation_rate,
                               const string &option, const unsigned long &external_reference);

    /**
     * @brief Gets the metacommunity reference for the provided parameters, or returns 0 if it doesn't exist
     * @param metacomm_parameters the parameters to check for
     * @return the metacommunity reference number, or 0 if it doesn't exist
     */
    unsigned long getReference(const MetacommunityParameters &metacomm_parameters);

    /**
     * @brief Adds an empty metacommunity parameters option.
     */
    void addNull();

    /**
     * @brief Determines if the array has any Metacommunity options (other than a null option).
     * @return
     */
    bool hasMetacommunityOption();
};

#endif //NECSIM_PARAMETERS_H
