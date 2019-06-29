//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SpeciesList.h
 * @brief  Contains the SpeciesList class for usage in coalescence simulations.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#include <iostream>
#include "SpeciesList.h"

SpeciesList::SpeciesList() : list_size(0), maxsize(0), next_active(0), species_id_list(), nwrap(0)
{

}

void SpeciesList::initialise(unsigned long maxsizein)
{
    maxsize = maxsizein;
    nwrap = 0;
    list_size = 0;
}

void SpeciesList::setMaxsize(unsigned long maxsizein)
{
    maxsize = maxsizein;
}

void SpeciesList::setSpecies(unsigned long index, unsigned long new_val)
{
#ifdef DEBUG
    if(index >= species_id_list.size())
    {
        throw out_of_range("List index to change value is out of range of vector size. Please report this bug.");
    }
#endif //DEBUG
    if(species_id_list[index] == 0)
    {
        throw runtime_error("List position to be replaced is zero. Check species_id_list assignment.");
    }
    species_id_list[index] = new_val;
}

void SpeciesList::setSpeciesEmpty(unsigned long index, unsigned long new_val)
{
    if(index >= species_id_list.size())
    {
        species_id_list.resize(index + 1, 0);
        species_id_list[index] = 0;
    }
    if(species_id_list[index] != 0)
    {
        throw runtime_error("List position to be replaced is not zero. Check species_id_list assignment.");
    }
    species_id_list[index] = new_val;
    list_size++;
}

void SpeciesList::setNext(unsigned long n)
{
    next_active = n;
}

void SpeciesList::setNwrap(unsigned long nr)
{
    nwrap = nr;
}

unsigned long SpeciesList::addSpecies(const unsigned long &new_spec)
{
#ifdef DEBUG
    if(list_size + 1 > maxsize)
    {
        stringstream ss;
        ss << "species_id_list size: " << list_size << endl;
        ss << "max size: " << maxsize << endl;
        writeLog(10, ss.str());
        throw out_of_range("Could not add species - species_id_list size greater than max size.");
    }
#endif
    // Check if there are empty spaces
    if(species_id_list.size() > list_size)
    {
        // First loop from the species_id_list size value
        for(unsigned long i = list_size; i < species_id_list.size(); i++)
        {
            if(species_id_list[i] == 0)
            {
                list_size++;
                species_id_list[i] = new_spec;
                return i;
            }
        }
        // Now loop over the rest of the lineages
        for(unsigned long i = 0; i < list_size; i++)
        {
            if(species_id_list[i] == 0)
            {
                list_size++;
                species_id_list[i] = new_spec;
                return i;
            }
        }
    }
    else
    {
        // Just need to append to the vector
        species_id_list.push_back(new_spec);
        list_size++;
        return species_id_list.size() - 1;
    }

    throw out_of_range("Could not add species - no empty space");
}

void SpeciesList::deleteSpecies(unsigned long index)
{
    species_id_list[index] = 0;
    list_size--;
}

void SpeciesList::decreaseNwrap()
{
    if(nwrap == 0)
    {
        throw runtime_error("Nwrap should never be decreased less than 0");
    }
    else if(nwrap == 1)
    {
        if(next_active != 0)
        {
            throw runtime_error("Nwrap is being set at 0 when an wrapped lineage is still present");
        }
    }
    nwrap--;
}

void SpeciesList::increaseListSize()
{
    list_size++;
}

void SpeciesList::increaseNwrap()
{
    nwrap++;
}

void SpeciesList::changePercentCover(unsigned long newmaxsize)
{
    maxsize = newmaxsize;
}

unsigned long SpeciesList::getRandLineage(shared_ptr<RNGController> rand_no)
{
    double rand_index;
    if(maxsize <= list_size)
    {
        // Then the species_id_list size is larger than the actual size. This means we must return a lineage.
        try
        {
            do
            {
                rand_index = rand_no->d01();
                rand_index *= species_id_list.size();
                //os << "ref: " << rand_index << ", " << species_id_list[round(rand_index)] << endl;
            }
            while(species_id_list[floor(rand_index)] == 0);
            //os << "RETURNING!" << endl;
            return (species_id_list[floor(rand_index)]);
        }
        catch(out_of_range &oor)
        {
            throw runtime_error("Listpos outside maxsize.");
        }
    }
    else
    {
        rand_index = rand_no->d01();
//		os << "rand_index: " << rand_index << endl;
        rand_index *= maxsize;
        if(rand_index >= species_id_list.size())
        {
            return 0;
        }

        auto i = static_cast<unsigned long>(floor(rand_index));

#ifdef DEBUG
        if(rand_index>maxsize)
            {
                stringstream ss;
                ss << "Random index is greater than the max size. Fatal error, please report this bug." << endl;
                throw runtime_error(ss.str());
            }
#endif // DEBUG
        return species_id_list[i];
    }
}

unsigned long SpeciesList::getSpecies(unsigned long index)
{
    return species_id_list[index];
}

unsigned long SpeciesList::getNext()
{
    return next_active;
}

unsigned long SpeciesList::getNwrap()
{
    return nwrap;
}

unsigned long SpeciesList::getListSize()
{
    return list_size;
}

unsigned long SpeciesList::getMaxSize()
{
    return maxsize;
}

unsigned long SpeciesList::getListLength()
{
    return species_id_list.size();
}

void SpeciesList::wipeList()
{
    next_active = 0;
    nwrap = 0;
    list_size = 0;
}

ostream &operator<<(ostream &os, const SpeciesList &r)
{
    os << r.species_id_list.size();
    return os;
}

istream &operator>>(istream &is, SpeciesList &r)
{
    unsigned int size;
    is >> size;
    r.species_id_list.resize(size, 0);
    return is;
}
