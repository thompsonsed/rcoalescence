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

namespace necsim
{
    SpeciesList::SpeciesList() : list_size(0), max_size(0), next_active(0), lineage_indices(), nwrap(0)
    {
    }

    SpeciesList &SpeciesList::operator=(const SpeciesList &other) noexcept
    {
        list_size = other.list_size;
        max_size = other.max_size;
        next_active = other.next_active;
        lineage_indices = other.lineage_indices;
        nwrap = other.nwrap;
        return *this;
    }

    SpeciesList &SpeciesList::operator=(SpeciesList &&other) noexcept
    {
        list_size = other.list_size;
        max_size = other.max_size;
        next_active = other.next_active;
        if(other.lineage_indices.empty())
        {
            lineage_indices.clear();
        }
        else
        {
            lineage_indices = std::move(other.lineage_indices);
        }
        nwrap = other.nwrap;
        return *this;
    }

    void SpeciesList::initialise(unsigned long maxsizein)
    {
        max_size = maxsizein;
        nwrap = 0;
        list_size = 0;
    }

    void SpeciesList::setMaxsize(unsigned long maxsizein)
    {
        max_size = maxsizein;
    }

    void SpeciesList::setSpecies(unsigned long index, unsigned long new_val)
    {
#ifdef DEBUG
        if(index >= lineage_indices.size())
        {
            throw std::out_of_range ("List index to change value is out of range of vector size. Please report this bug.");
        }
#endif //DEBUG
        if(lineage_indices[index] == 0)
        {
            throw std::runtime_error("List position to be replaced is zero. Check lineage_indices assignment.");
        }
        lineage_indices[index] = new_val;
    }

    void SpeciesList::setSpeciesEmpty(unsigned long index, unsigned long new_val)
    {
        if(index >= lineage_indices.size())
        {
            lineage_indices.resize(index + 1, 0);
            lineage_indices[index] = 0;
        }
        if(lineage_indices[index] != 0)
        {
            throw std::runtime_error("List position to be replaced is not zero. Check lineage_indices assignment.");
        }
        lineage_indices[index] = new_val;
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
        if(list_size + 1 > max_size)
        {
            std::stringstream ss;
            ss << "lineage_indices size: " << list_size << std::endl;
            ss << "max size: " << max_size << std::endl;
            writeLog(10, ss.str());
            throw std::out_of_range ("Could not add species - lineage_indices size greater than max size.");
        }
#endif
        // Check if there are empty spaces
        if(lineage_indices.size() > list_size)
        {
            // First loop from the lineage_indices size value
            for(unsigned long i = list_size; i < lineage_indices.size(); i++)
            {
                if(lineage_indices[i] == 0)
                {
                    list_size++;
                    lineage_indices[i] = new_spec;
                    return i;
                }
            }
            // Now loop over the rest of the lineages
            for(unsigned long i = 0; i < list_size; i++)
            {
                if(lineage_indices[i] == 0)
                {
                    list_size++;
                    lineage_indices[i] = new_spec;
                    return i;
                }
            }
        }
        else
        {
            // Just need to append to the vector
            lineage_indices.push_back(new_spec);
            list_size++;
            return lineage_indices.size() - 1;
        }

        throw std::out_of_range("Could not add species - no empty space");
    }

    void SpeciesList::deleteSpecies(unsigned long index)
    {
        lineage_indices[index] = 0;
        list_size--;
    }

    void SpeciesList::decreaseNwrap()
    {
        if(nwrap == 0)
        {
            throw std::runtime_error("Nwrap should never be decreased less than 0");
        }
        else if(nwrap == 1)
        {
            if(next_active != 0)
            {
                throw std::runtime_error("Nwrap is being set at 0 when an wrapped lineage is still present");
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
        max_size = newmaxsize;
    }

    unsigned long SpeciesList::getRandLineage(const shared_ptr<RNGController> &rand_no)
    {
        double rand_index;
        if(max_size <= list_size)
        {
            // Then the lineage_indices size is larger than the actual size. This means we must return a lineage.
            try
            {
                do
                {
                    rand_index = rand_no->d01();
                    rand_index *= lineage_indices.size();
                    //os << "ref: " << rand_index << ", " << lineage_indices[round(rand_index)] << std::endl;
                }
                while(lineage_indices[floor(rand_index)] == 0);
                //os << "RETURNING!" << std::endl;
                return (lineage_indices[floor(rand_index)]);
            }
            catch(std::out_of_range &oor)
            {
                throw std::runtime_error("Listpos outside max_size.");
            }
        }
        else
        {
            rand_index = rand_no->d01();
            //		os << "rand_index: " << rand_index << std::endl;
            rand_index *= max_size;
            if(rand_index >= lineage_indices.size())
            {
                return 0;
            }

            auto i = static_cast<unsigned long>(floor(rand_index));

#ifdef DEBUG
            if(rand_index > max_size)
            {
                std::stringstream ss;
                ss << "Random index is greater than the max size. Fatal error, please report this bug." << std::endl;
                throw std::runtime_error(ss.str());
            }
#endif // DEBUG
            return lineage_indices[i];
        }
    }

    unsigned long SpeciesList::getLineageIndex(unsigned long index) const
    {
        return lineage_indices[index];
    }

    unsigned long SpeciesList::getNext() const
    {
        return next_active;
    }

    unsigned long SpeciesList::getNwrap() const
    {
        return nwrap;
    }

    unsigned long SpeciesList::getListSize() const
    {
        return list_size;
    }

    unsigned long SpeciesList::getMaxSize() const
    {
        return max_size;
    }

    unsigned long SpeciesList::getListLength() const
    {
        return lineage_indices.size();
    }

    void SpeciesList::wipeList()
    {
        next_active = 0;
        nwrap = 0;
        list_size = 0;
    }

    double SpeciesList::getCoalescenceProbability() const
    {
        return std::min(double(list_size) / double(max_size), 1.0);
    }

    std::ostream &operator<<(std::ostream &os, const SpeciesList &r)
    {
        os << r.lineage_indices.size();
        return os;
    }

    std::istream &operator>>(std::istream &is, SpeciesList &r)
    {
        unsigned int size;
        is >> size;
        r.lineage_indices.resize(size, 0);
        return is;
    }
}