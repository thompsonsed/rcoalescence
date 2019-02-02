//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @file parameters.cpp
 *
 * @brief Contains various objects used for storing parameters related to tree reconstruction.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include "parameters.h"
#include "custom_exceptions.h"

CommunityParameters::CommunityParameters(unsigned long reference_in, long double speciation_rate_in,
                                         long double time_in,
                                         bool fragment_in, unsigned long metacommunity_reference_in,
                                         const ProtractedSpeciationParameters &protracted_params)
{
    setup(reference_in, speciation_rate_in, time_in, fragment_in, metacommunity_reference_in, protracted_params);
}

void CommunityParameters::setup(unsigned long reference_in, long double speciation_rate_in, long double time_in,
                                bool fragment_in,
                                unsigned long metacommunity_reference_in,
                                const ProtractedSpeciationParameters &protracted_params)
{
    time = time_in;
    speciation_rate = speciation_rate_in;
    fragment = fragment_in;
    reference = reference_in;
    metacommunity_reference = metacommunity_reference_in;
    protracted_parameters = protracted_params;
    updated = false;

}

bool CommunityParameters::compare(long double speciation_rate_in, long double time_in, bool fragment_in,
                                  unsigned long metacommunity_reference_in,
                                  const ProtractedSpeciationParameters &protracted_params)
{
    if(doubleCompare(double(time_in), double(0.0), 0.00001))
    {
#ifdef DEBUG
        stringstream os;
        os << "Detected generation at t=0.0." << endl;
        writeLog(10, os);
#endif // DEBUG
        return doubleCompare(speciation_rate, speciation_rate_in, speciation_rate * 0.000001) &&
               fragment == fragment_in && metacommunity_reference == metacommunity_reference_in &&
               protracted_params == protracted_parameters;
    }
    return doubleCompare(speciation_rate, speciation_rate_in, speciation_rate * 0.000001) &&
           doubleCompare(time, time_in, time * 0.0000000001) && fragment == fragment_in &&
           metacommunity_reference == metacommunity_reference_in &&
           protracted_params == protracted_parameters;
}

bool CommunityParameters::compare(long double speciation_rate_in, long double time_in,
                                  unsigned long metacommunity_reference_in,
                                  const ProtractedSpeciationParameters &protracted_params)
{
    return doubleCompare(speciation_rate, speciation_rate_in, speciation_rate * 0.000001) &&
           doubleCompare(time, time_in, 0.0000000001) && metacommunity_reference == metacommunity_reference_in &&
           protracted_params == protracted_parameters;
}

bool CommunityParameters::compare(unsigned long reference_in)
{
    return reference == reference_in;
}

CommunitiesArray::CommunitiesArray() : comm_parameters(){};

void CommunitiesArray::pushBack(unsigned long reference, long double speciation_rate, long double time, bool fragment,
                                unsigned long metacommunity_reference,
                                const ProtractedSpeciationParameters &protracted_params)
{
    auto tmp_param = make_shared<CommunityParameters>(reference, speciation_rate, time, fragment,
                                                      metacommunity_reference, protracted_params);
    comm_parameters.emplace_back(tmp_param);
}

void CommunitiesArray::pushBack(shared_ptr<CommunityParameters> tmp_param)
{
    comm_parameters.emplace_back(tmp_param);
}

shared_ptr<CommunityParameters> CommunitiesArray::addNew(long double speciation_rate, long double time, bool fragment,
                                                         unsigned long metacommunity_reference,
                                                         const ProtractedSpeciationParameters &protracted_params)
{
    unsigned long max_reference = 1;
    for(auto &i : comm_parameters)
    {
        if(i->compare(speciation_rate, time, metacommunity_reference, protracted_params))
        {
            if(i->fragment == fragment || !fragment)
            {
                stringstream ss;
                ss << "Non-unique parameter set: " << endl;
                ss << "-speciation rate: " << speciation_rate << endl;
                ss << "-time: " << time << endl;
                ss << "-fragment: " << fragment << endl;
                ss << "-metacommunity reference: " << metacommunity_reference << endl;
                ss << "-protracted speciation min: " << protracted_params.min_speciation_gen << endl;
                ss << "-protracted speciation max: " << protracted_params.max_speciation_gen << endl;
                writeCritical(ss.str());
                throw FatalException("Tried to get reference for non-unique parameter set in communities. "
                                     "Please report this bug.");
            }
            else
            {
                i->fragment = true;
                i->updated = true;
                return i;
            }
        }
        else
        {
            if(i->reference >= max_reference)
            {
                max_reference = i->reference + 1;
            }
        }
    }
    auto tmp_param = make_shared<CommunityParameters>(max_reference, speciation_rate, time, fragment,
                                                      metacommunity_reference, protracted_params);
    comm_parameters.emplace_back(tmp_param);
    return comm_parameters.back();
}

bool CommunitiesArray::hasPair(long double speciation_rate, double time, bool fragment,
                               unsigned long metacommunity_reference,
                               const ProtractedSpeciationParameters &protracted_params)
{
    for(auto &i : comm_parameters)
    {
        if(i->compare(speciation_rate, time, fragment, metacommunity_reference, protracted_params))
        {
            return true;
        }
    }
    return false;
}

MetacommunityParameters::MetacommunityParameters()
{
    reference = 0;
    metacommunity_size = 0;
    speciation_rate = 0.0;
    option = "none";
    external_reference = 0;
}

MetacommunityParameters::MetacommunityParameters(const unsigned long &reference_in,
                                                 const unsigned long &metacommunity_size_in,
                                                 const long double &speciation_rate_in, const string &option_in,
                                                 const unsigned long &external_reference_in)
{
    reference = reference_in;
    metacommunity_size = metacommunity_size_in;
    speciation_rate = speciation_rate_in;
    option = option_in;
    external_reference = external_reference_in;
}

bool MetacommunityParameters::compare(unsigned long metacommunity_size_in, long double speciation_rate_in,
                                      const string &option_in, const unsigned long &ext_reference_in)
{
    if(option_in == "none" && metacommunity_size == 0)
    {
        return option == "none" && metacommunity_size == 0;
    }
    if(option_in == "simulated" || option_in == "analytical")
    {
        return doubleCompare(speciation_rate, speciation_rate_in, speciation_rate * 0.000001) &&
               metacommunity_size == metacommunity_size_in && option == option_in;
    }
    else
    {
        return option == option_in && external_reference == ext_reference_in;
    }
}

bool MetacommunityParameters::compare(const MetacommunityParameters &metacomm_in)
{
    return compare(metacomm_in.metacommunity_size, metacomm_in.speciation_rate, metacomm_in.option,
                   metacomm_in.external_reference);
}

bool MetacommunityParameters::compare(unsigned long reference_in)
{
    return reference == reference_in;
}

bool MetacommunityParameters::isMetacommunityOption() const
{
    if(metacommunity_size > 0 && speciation_rate > 0)
    {
        if(option != "simulated" && option != "analytical")
        {
            throw FatalException("Must use simulation or analytical solutions when supplying a metacommunity size and "
                                 "speciation rate.");
        }
        return true;
    }
    if(option != "simulated" && option != "analytical" && option != "none")
    {
        if(external_reference == 0)
        {
            throw FatalException("External reference cannot be 0 when reading metacommunity from external database.");
        }
        return true;
    }
    return false;
}

void MetacommunityParameters::clear()
{
    reference = 0;
    metacommunity_size = 0;
    speciation_rate = 0.0;
    option = "none";
    external_reference = 0;
}

MetacommunityParameters &MetacommunityParameters::operator=(const MetacommunityParameters &parameters)
{
    reference = parameters.reference;
    metacommunity_size = parameters.metacommunity_size;
    speciation_rate = parameters.speciation_rate;
    option = parameters.option;
    external_reference = parameters.external_reference;
    return *this;
}

bool MetacommunityParameters::operator==(const MetacommunityParameters &parameters)
{
    return compare(parameters);
}

bool MetacommunityParameters::operator!=(const MetacommunityParameters &parameters)
{
    return !compare(parameters);
}

MetacommunitiesArray::MetacommunitiesArray() : metacomm_parameters()
{

}

vector<shared_ptr<MetacommunityParameters>>::iterator MetacommunitiesArray::begin()
{
    return metacomm_parameters.begin();
}

vector<shared_ptr<MetacommunityParameters>>::iterator MetacommunitiesArray::end()
{
    return metacomm_parameters.end();
}

vector<shared_ptr<MetacommunityParameters>>::const_iterator MetacommunitiesArray::begin() const
{
    return metacomm_parameters.begin();
}

vector<shared_ptr<MetacommunityParameters>>::const_iterator MetacommunitiesArray::end() const
{
    return metacomm_parameters.end();
}

void MetacommunitiesArray::pushBack(const unsigned long &reference, const unsigned long &metacommunity_size,
                                    const long double &speciation_rate, const string &option,
                                    const unsigned long &external_reference)
{
    auto tmp_param = make_shared<MetacommunityParameters>(reference, speciation_rate, metacommunity_size, option,
                                                          external_reference);
    metacomm_parameters.emplace_back(tmp_param);
}

void MetacommunitiesArray::pushBack(shared_ptr<MetacommunityParameters> tmp_param)
{
    metacomm_parameters.emplace_back(tmp_param);
}

void MetacommunitiesArray::clear()
{
    metacomm_parameters.clear();
}

unsigned long MetacommunitiesArray::size()
{
    return metacomm_parameters.size();
}

bool MetacommunitiesArray::empty()
{
    return metacomm_parameters.empty();
}

unsigned long MetacommunitiesArray::addNew(const unsigned long &metacommunity_size, const long double &speciation_rate,
                                           const string &option, const unsigned long &external_reference)
{
    unsigned long max_reference = 1;
    for(auto &i : metacomm_parameters)
    {
        if(i->compare(metacommunity_size, speciation_rate, option, external_reference))
        {
            stringstream ss;
            ss << "Duplicate metacommunity parameters: " << endl;
            ss << "Size: " << metacommunity_size << endl;
            ss << "Speciation rate: " << speciation_rate << endl;
            ss << "Option: " << option << endl;
            ss << "External reference: " << external_reference << endl;
            writeCritical(ss.str());
            throw FatalException("Tried to get reference for non-unique parameter set in metacommunities. "
                                 "Please report this bug.");
        }
        else
        {
            if(i->reference >= max_reference)
            {
                max_reference = i->reference + 1;
            }
        }
    }
    auto tmp_param = make_shared<MetacommunityParameters>(max_reference, metacommunity_size, speciation_rate, option,
                                                          external_reference);
    metacomm_parameters.emplace_back(tmp_param);
    return max_reference;
}

unsigned long MetacommunitiesArray::addNew(const MetacommunityParameters &metacomm_in)
{
    return addNew(metacomm_in.metacommunity_size, metacomm_in.speciation_rate, metacomm_in.option,
                  metacomm_in.external_reference);
}

bool MetacommunitiesArray::hasOption(const unsigned long &metacommunity_size, const long double &speciation_rate,
                                     const string &option, const unsigned long &external_reference)
{
    for(auto &i : metacomm_parameters)
    {
        if(i->compare(metacommunity_size, speciation_rate, option, external_reference))
        {
            return true;
        }
    }
    return false;
}

bool MetacommunitiesArray::hasOption(unsigned long reference)
{
    for(auto &i : metacomm_parameters)
    {
        if(i->compare(reference))
        {
            return true;
        }
    }
    return false;
}

bool MetacommunitiesArray::hasOption(const MetacommunityParameters &metacomm_in)
{
    return hasOption(metacomm_in.metacommunity_size, metacomm_in.speciation_rate, metacomm_in.option,
                     metacomm_in.external_reference);
}

unsigned long MetacommunitiesArray::getReference(const unsigned long &metacommunity_size,
                                                 const long double &speciation_rate, const string &option,
                                                 const unsigned long &external_reference)
{
    if(metacommunity_size == 0)
    {
        if(option == "simulated" || option == "analytical" || option == "none")
        {
            return 0;
        }
        else
        {
            if(external_reference == 0)
            {
                return 0;
            }
        }
    }
    for(auto &i : metacomm_parameters)
    {
        if(i->compare(metacommunity_size, speciation_rate, option, external_reference))
        {
            return i->reference;
        }
    }
    return 0;
}

unsigned long MetacommunitiesArray::getReference(const MetacommunityParameters &metacomm_parameters)
{
    return getReference(metacomm_parameters.metacommunity_size, metacomm_parameters.speciation_rate,
                        metacomm_parameters.option, metacomm_parameters.external_reference);
}

void MetacommunitiesArray::addNull()
{
    MetacommunityParameters null_parameters;
    addNew(null_parameters);
}

bool MetacommunitiesArray::hasMetacommunityOption()
{
    MetacommunityParameters null_parameters;
    for(const auto & item: metacomm_parameters)
    {
        if(*item != null_parameters)
        {
            return true;
        }
    }
    return false;
}

