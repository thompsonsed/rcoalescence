#include <utility>

//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 07/08/2017
 * @file DispersalCoordinator.cpp
 * @brief Contains the DispersalCoordinator, which contains all routines related to dispersal
 * including utilisation of density maps and dispersal probability maps.
 * 
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#include "DispersalCoordinator.h"

namespace necsim
{
    DispersalCoordinator::DispersalCoordinator() : dispersal_prob_map(), raw_dispersal_prob_map(), NR(nullptr),
                                                   landscape(make_shared<Landscape>()),
                                                   reproduction_map(make_shared<ActivityMap>()), generation(nullptr),
                                                   doDispersal(nullptr), checkEndPointFptr(nullptr), xdim(0), ydim(0),
                                                   full_dispersal_map(false)
    {

    }

    DispersalCoordinator::~DispersalCoordinator() = default;

    void DispersalCoordinator::setRandomNumber(shared_ptr<RNGController> NR_ptr)
    {
        NR = std::move(NR_ptr);
    }

    bool DispersalCoordinator::isFullDispersalMap() const
    {
        return full_dispersal_map;
    }

    void DispersalCoordinator::setMaps(const shared_ptr<Landscape> &landscape_ptr, shared_ptr<ActivityMap> repr_map_ptr)
    {
        landscape = landscape_ptr;
        reproduction_map = std::move(repr_map_ptr);
        if(!landscape_ptr)
        {
            throw FatalException("Attempting to set map pointer to null pointer in DispersalCoordinator.");
        }
        xdim = landscape->getSimParameters()->fine_map_x_size;
        ydim = landscape->getSimParameters()->fine_map_y_size;
    }

    void DispersalCoordinator::setMaps(const shared_ptr<Landscape> &landscape_ptr)
    {
        setMaps(landscape_ptr, make_shared<ActivityMap>());
    }

    void DispersalCoordinator::setGenerationPtr(double* generation_ptr)
    {
        generation = generation_ptr;
    }

    void DispersalCoordinator::setDispersal(const string &dispersal_method,
                                            const string &dispersal_file,
                                            const unsigned long &dispersal_x,
                                            const unsigned long &dispersal_y,
                                            const double &m_probin,
                                            const double &cutoffin,
                                            const double &sigmain,
                                            const double &tauin,
                                            const bool &restrict_self)
    {
        full_dispersal_map = false;
        // Open our file connection
        if((dispersal_file == "none") | dispersal_file.empty())
        {
            if(!NR)
            {
                throw FatalException("Random number generator pointer has not been set in DispersalCoordinator.");
            }
            writeInfo("Using dispersal kernel.\n");
            setEndPointFptr(restrict_self);
            NR->setDispersalParams(sigmain, tauin);
            NR->setDispersalMethod(dispersal_method, m_probin, cutoffin);
            doDispersal = &DispersalCoordinator::disperseDensityMap;
            reproduction_map->standardiseValues();
        }
        else if(dispersal_file == "null")
        {
            writeInfo("Using null dispersal file.\n");
            doDispersal = &DispersalCoordinator::disperseNullDispersalMap;
            reproduction_map->standardiseValues();
        }
        else
        {
            writeInfo("Using dispersal file.\n");
            doDispersal = &DispersalCoordinator::disperseDispersalMap;
            importDispersal(dispersal_x * dispersal_y, dispersal_file);
            full_dispersal_map = true;
        }
    }

    void DispersalCoordinator::setDispersal(shared_ptr<SimParameters> simParameters)
    {
        if(!simParameters)
        {
            throw FatalException(
                    "Simulation current_metacommunity_parameters pointer has not been set for DispersalCoordinator.");
        }
        setDispersal(simParameters->dispersal_method,
                     simParameters->dispersal_file,
                     simParameters->fine_map_x_size,
                     simParameters->fine_map_y_size,
                     simParameters->m_prob,
                     simParameters->cutoff,
                     simParameters->sigma,
                     simParameters->tau,
                     simParameters->restrict_self);
    }

    void DispersalCoordinator::importDispersal(const unsigned long &dispersal_dim, const string &dispersal_file)
    {
        // Check file existence
        ifstream infile(dispersal_file);
        if(!infile.good())
        {
            string msg =
                    "Could not access dispersal map file " + dispersal_file + ". Check file exists and is readable.";
            throw FatalException(msg);
        }
        infile.close();
        dispersal_prob_map.setSize(dispersal_dim, dispersal_dim);
        dispersal_prob_map.import(dispersal_file);
        if(landscape->hasHistorical())
        {
            setRawDispersalMap();
        }
        *generation = 0.0;
        addDensity();
        addReproduction();
        fixDispersal();
        dispersal_prob_map.close();
        verifyDispersalMapSetup();
    }

    void DispersalCoordinator::setRawDispersalMap()
    {
        raw_dispersal_prob_map = dispersal_prob_map;
    }

    void DispersalCoordinator::addDensity()
    {
        for(unsigned long i = 0; i < ydim; i++)
        {
            for(unsigned long j = 0; j < xdim; j++)
            {
                unsigned long index = j + i * xdim;
                for(unsigned long k = 0; k < dispersal_prob_map.getRows(); k++)
                {
                    auto density = landscape->getValFine(j, i, *generation);
                    if(dispersal_prob_map.get(k, index) > 0.0 && density == 0)
                    {
                        Step origin_step;
                        calculateCellCoordinates(origin_step, k);
                        stringstream ss;
                        Step destination_step;
                        calculateCellCoordinates(destination_step, j + (i * xdim));
                        ss << "Dispersal from " << origin_step.x << ", " << origin_step.y << " (";
                        ss << origin_step.xwrap << ", " << origin_step.ywrap << ") to ";
                        ss << destination_step.x << ", " << destination_step.y << " (" << destination_step.xwrap;
                        ss << ", " << destination_step.ywrap << ")" << endl;
                        ss << "Source row: " << k << " destination row: " << index << endl;
                        ss << "Dispersal map value: " << dispersal_prob_map.get(k, index) << endl;
                        ss << "Origin density: "
                           << landscape->getVal(origin_step.x, origin_step.y, origin_step.xwrap, origin_step.ywrap, 0.0)
                           << endl;
                        ss << "Destination density: " << landscape->getValFine(j, i, *generation) << endl;
                        writeError(ss.str());
                        throw FatalException("Dispersal map is non zero where density is 0.");
                    }
                    dispersal_prob_map.get(k, index) *= density;
                }
            }
        }
    }

    void DispersalCoordinator::addReproduction()
    {
        if(reproduction_map != nullptr)
        {
            if(!reproduction_map->isNull())
            {
                for(unsigned long i = 0; i < ydim; i++)
                {
                    for(unsigned long j = 0; j < xdim; j++)
                    {
                        unsigned long index = j + i * xdim;
                        for(unsigned long k = 0; k < dispersal_prob_map.getRows(); k++)
                        {

                            dispersal_prob_map.get(k, index) *= reproduction_map->get(i, j);
                        }
                    }
                }
            }
        }
    }

    void DispersalCoordinator::fixDispersal()
    {
        for(unsigned long row = 0; row < dispersal_prob_map.getRows(); row++)
        {
            fixDispersalRow(row);
        }
    }

    void DispersalCoordinator::fixDispersalRow(unsigned long row)
    {
        // check if the row needs to be fixed
        if(checkDispersalRow(row))
        {
            double total_value = 0.0;
            for(unsigned long i = 0; i < dispersal_prob_map.getCols(); i++)
            {
                total_value += dispersal_prob_map.get(row, i);
            }
            if(total_value == 0.0)
            {
                return;
            }
            dispersal_prob_map.get(row, 0) = dispersal_prob_map.get(row, 0) / total_value;
            for(unsigned long i = 1; i < dispersal_prob_map.getCols(); i++)
            {
                dispersal_prob_map.get(row, i) =
                        dispersal_prob_map.get(row, i - 1) + (dispersal_prob_map.get(row, i) / total_value);
            }
#ifdef DEBUG
            if(checkDispersalRow(row))
            {
                throw FatalException("Dispersal probability map not correctly fixed to sum to 1.0.");
            }
#endif // DEBUG
        }
    }

    bool DispersalCoordinator::checkDispersalRow(unsigned long row)
    {
        if(abs(dispersal_prob_map.get(row, dispersal_prob_map.getCols() - 1) - 1.0) > 0.00000001)
        {
            return true;
        }
        for(unsigned long i = 0; i < dispersal_prob_map.getCols() - 1; i++)
        {
            if(dispersal_prob_map.get(row, i) > dispersal_prob_map.get(row, i + 1))
            {
                return true;
            }
        }
        return false;
    }

    void DispersalCoordinator::verifyDispersalMapSetup()
    {
        if(dispersal_prob_map.getCols() > 0)
        {
            writeInfo("Verifying dispersal setup...\n");
            if(dispersal_prob_map.getCols() != dispersal_prob_map.getRows())
            {
                throw FatalException("Dispersal probability map dimensions do not match.");
            }
            bool has_printed = false;
            for(unsigned long y = 0; y < dispersal_prob_map.getRows(); y++)
            {
                Step origin_step;
                calculateCellCoordinates(origin_step, y);
#ifdef DEBUG
                assertReferenceMatches(y);
#endif // DEBUG
                bool origin_value =
                        landscape->getVal(origin_step.x, origin_step.y, origin_step.xwrap, origin_step.ywrap, 0.0) > 0;
                double dispersal_total = 0.0;
                for(unsigned long x = 0; x < dispersal_prob_map.getCols(); x++)
                {
                    Step destination_step;
                    calculateCellCoordinates(destination_step, x);
#ifdef DEBUG
                    assertReferenceMatches(x);
#endif // DEBUG
                    bool destination_value = landscape->getVal(destination_step.x,
                                                               destination_step.y,
                                                               destination_step.xwrap,
                                                               destination_step.ywrap,
                                                               0.0) > 0;
                    double dispersal_prob;
                    if(x == 0)
                    {
                        dispersal_prob = dispersal_prob_map.get(y, 0);
                    }
                    else
                    {
                        dispersal_prob = dispersal_prob_map.get(y, x) - dispersal_prob_map.get(y, x - 1);
                    }
                    dispersal_total += dispersal_prob;
                    if(dispersal_prob > 0.0)
                    {
                        if(!destination_value && origin_value)
                        {
                            stringstream ss;
                            ss << "Dispersal from " << origin_step.x << ", " << origin_step.y << " (";
                            ss << origin_step.xwrap << ", " << origin_step.ywrap << ") to ";
                            ss << destination_step.x << ", " << destination_step.y << " (" << destination_step.xwrap;
                            ss << ", " << destination_step.ywrap << ")" << endl;
                            ss << "Source row: " << y << " destination row: " << x << endl;
                            ss << "Dispersal map value: " << dispersal_prob << endl;
                            ss << "Origin density: " << landscape->getVal(origin_step.x,
                                                                          origin_step.y,
                                                                          origin_step.xwrap,
                                                                          origin_step.ywrap,
                                                                          0.0) << endl;
                            ss << "Destination density: " << landscape->getVal(destination_step.x,
                                                                               destination_step.y,
                                                                               destination_step.xwrap,
                                                                               destination_step.ywrap,
                                                                               0.0) << endl;
                            writeError(ss.str());
                            throw FatalException("Dispersal map is non zero where density is 0.");
                        }
                        if(!origin_value && !has_printed)
                        {
                            has_printed = true;
                            writeWarning("Dispersal values exist for non-zero density values.");
                        }
                    }
                }
                if(dispersal_total == 0.0 && origin_value)
                {
                    stringstream ss;
                    ss << "No dispersal probabilities from cell at " << origin_step.x << ", " << origin_step.y;
                    ss << " (" << origin_step.xwrap << ", " << origin_step.ywrap;
                    ss << ") to any other cell, despite non-zero density." << endl;
                    throw FatalException(ss.str());
                }
            }
        }
    }

    void DispersalCoordinator::updateDispersalMap()
    {
        if(dispersal_prob_map.getRows() > 0)
        {
            dispersal_prob_map = raw_dispersal_prob_map;
            addDensity();
            addReproduction();
            fixDispersal();
        }
    }

#ifdef DEBUG

    void DispersalCoordinator::assertReferenceMatches(unsigned long expected)
    {
        unsigned long row_ref = expected;
        Step step;
        calculateCellCoordinates(step, row_ref);
        auto actual = calculateCellReference(step);
        if(actual != expected)
        {
            stringstream ss;
            ss << "Expected reference " << expected << endl;
            ss << "Actual reference " << actual << endl;
            ss << "Coordinates: " << step.x << ", " << step.y << "(" << step.xwrap << ", ";
            ss << step.ywrap << ")" << endl;
            ss << "Converted values: " << landscape->convertSampleXToFineX(step.x, step.xwrap);
            ss << ", " << landscape->convertSampleYToFineY(step.y, step.ywrap) << endl;
            throw FatalException(ss.str());
        }
    }

        void DispersalCoordinator::validateNoSelfDispersalInDispersalMap()
    {
        Cell cell;
        try
        {
            if(dispersal_prob_map.get(0, 0) > 0.0)
            {
                cell.x = 0;
                cell.y = 0;
                throw FatalException("Self dispersal non-zero.");
            }
            for(unsigned long i = 1; i < dispersal_prob_map.getCols(); i++)
            {
                Step tmp_step;
                calculateCellCoordinates(tmp_step, i);
                if(dispersal_prob_map.get(i, i - 1) > dispersal_prob_map.get(i, i)
                   && dispersal_prob_map.get(i, i) > 0.0)
                {
                    cell.x = tmp_step.x;
                    cell.y = tmp_step.y;
                    throw FatalException("Self dispersal non-zero.");
                }
                const Step orig_step = tmp_step;
                // Generate 1000 random events to check for self-dispersal
                for(unsigned long j = 0; j < 1000; j++)
                {
                    disperseDispersalMap(tmp_step);
                    if(tmp_step == orig_step)
                    {
                        stringstream ss2;
                        ss2 << "Self dispersal possible: " << orig_step << " to " << tmp_step << endl;
                        throw FatalException(ss2.str());
                    }
                    tmp_step = orig_step;
                }
            }
        }
        catch(FatalException &fe)
        {
            stringstream ss;
            unsigned long index = calculateCellIndex(cell);
            ss << "Cell at " << cell.x << ", " << cell.y << " has incorrect self-dispersal assignment: " << fe.what()
               << endl;
            ss << "Dispersal value: " << dispersal_prob_map.get(index, index) << endl;
            if(index > 0)
            {
                ss << "Prior value: " << dispersal_prob_map.get(index, index - 1) << endl;
            }
            throw FatalException(ss.str());
        }

    }


#endif // DEBUG

    void DispersalCoordinator::disperseNullDispersalMap(Step &this_step)
    {
        // Pick a random cell - that's all we need
        long rand_x;
        long rand_y;
        // rejection sample based on reproduction values.
        do
        {
            rand_x = floor(NR->d01() * (xdim - 1));
            rand_y = floor(NR->d01() * (ydim - 1));
        }
        while(!reproduction_map->actionOccurs(rand_x, rand_y, 0, 0));
        calculateCellCoordinates(this_step, rand_x + rand_y * xdim);
    }

    void DispersalCoordinator::disperseDispersalMap(Step &this_step)
    {
        // Generate random number 0-1
        double random_no = NR->d01();
        // Now find the cell with that value
        // Now we get the cell reference
        unsigned long row_ref = calculateCellReference(this_step);
        auto begin = dispersal_prob_map.begin() + dispersal_prob_map.index(row_ref, 0);
        auto end = dispersal_prob_map.begin() + dispersal_prob_map.index(row_ref, dispersal_prob_map.getCols());
        unsigned long out_col = std::lower_bound(begin, end, random_no) - begin;

//        // Interval bisection on the cells to get the dispersal value
//        unsigned long min_col = 0;
//        unsigned long max_col = dispersal_prob_map.getCols() - 1;
//        while(max_col - min_col > 1)
//        {
//            auto to_check = static_cast<unsigned long>(floor(double(max_col - min_col) / 2.0) + min_col);
//            if(dispersal_prob_map.get(row_ref, to_check) < random_no)
//            {
//                min_col = to_check;
//            }
//            else
//            {
//                max_col = to_check;
//            }
//        }
        // Now get the coordinates of our cell reference
        calculateCellCoordinates(this_step, out_col);
#ifdef DEBUG
        if(landscape->getVal(this_step.x, this_step.y, this_step.xwrap, this_step.ywrap, *generation) < 1.0)
        {
            stringstream ss;
            ss << "Dispersal attempted to cell of zero density " << this_step.x << ", " << this_step.y;
            ss << " (" << this_step.xwrap << ", " << this_step.ywrap << ")" << endl;
            throw FatalException(ss.str());
        }
#endif // DEBUG
    }

    void DispersalCoordinator::calculateCellCoordinates(Step &this_step, const unsigned long &col_ref)
    {
        this_step.x = long(floor(fmod(double(col_ref), xdim)));
        this_step.y = long(floor(double(col_ref) / xdim));
        this_step.xwrap = 0;
        this_step.ywrap = 0;
        // Convert back to sample map
        landscape->convertFineToSample(this_step.x, this_step.xwrap, this_step.y, this_step.ywrap);

    }

    unsigned long DispersalCoordinator::calculateCellReference(Step &this_step) const
    {
        Cell cell;
        cell.x = landscape->convertSampleXToFineX(this_step.x, this_step.xwrap);
        cell.y = landscape->convertSampleYToFineY(this_step.y, this_step.ywrap);
        return calculateCellIndex(cell);
    }

    unsigned long DispersalCoordinator::calculateCellIndex(const Cell &cell) const
    {
        return cell.x + (cell.y * xdim);
    }

    void DispersalCoordinator::disperseDensityMap(Step &this_step)
    {
        bool fail;
        fail = true;
        // Store the starting positions
        long startx, starty, startxwrap, startywrap;
        startx = this_step.x;
        starty = this_step.y;
        startxwrap = this_step.xwrap;
        startywrap = this_step.ywrap;
        // keep looping until we reach a viable place to move from.
        // Store the density in the end location.
        unsigned long density;
        double dist, angle;
        // First check to see if the source cell is non-habitat.
        // If this is the case, then find the nearest neighbouring habitat pixel and only pick dispersal distances
        // greater than or equal to that distance.
        if(!landscape->getVal(this_step.x, this_step.y, this_step.xwrap, this_step.ywrap, *generation))
        {
            auto min_distance = landscape->distanceToNearestHabitat(this_step.x,
                                                                    this_step.y,
                                                                    this_step.xwrap,
                                                                    this_step.ywrap,
                                                                    *generation);
#ifdef DEBUG
            if(!landscape->isOnMap(this_step.x, this_step.y + min_distance, this_step.xwrap, this_step.ywrap)
               && !landscape->isOnMap(this_step.x + min_distance, this_step.y, this_step.xwrap, this_step.ywrap)
               && !landscape->isOnMap(this_step.x - min_distance, this_step.y, this_step.xwrap, this_step.ywrap)
               && !landscape->isOnMap(this_step.x, this_step.y - min_distance, this_step.xwrap, this_step.ywrap))
            {
                stringstream ss;
                ss << "Minimum distance calculated of " << min_distance << " for cell at x, y (" << this_step.x;
                ss << ", " << this_step.y << ") and wrap (" << this_step.xwrap << ", " << this_step.ywrap << ")";
                ss << " is outside of bounds of map at generation " << *generation << endl;
                ss << "Please report this bug." << endl;
                throw FatalException(ss.str());
            }
#endif // DEBUG
            unsigned long counter = 0;
            while(fail)
            {
                counter++;
                dist = NR->dispersalMinDistance(min_distance);
#ifdef DEBUG
                if(dist < min_distance)
                {
                    throw FatalException("Distance is less than minimum distance: please report this bug.");
                }
#endif // DEBUG
                angle = NR->direction();
                density = landscape->runDispersal(dist,
                                                  angle,
                                                  this_step.x,
                                                  this_step.y,
                                                  this_step.xwrap,
                                                  this_step.ywrap,
                                                  fail,
                                                  *generation);
                if(!fail)
                {
                    fail = !checkEndPoint(density,
                                          this_step.x,
                                          this_step.y,
                                          this_step.xwrap,
                                          this_step.ywrap,
                                          startx,
                                          starty,
                                          startxwrap,
                                          startywrap);
                }
                // This is a hack for those scenarios where habitat disappears and there is no easy replacement - then
                // the parent just comes from a nearest habitat cell that exists.
                if(counter > 10000000)
                {
#ifdef DEBUG
                    stringstream ss;
                    ss << "No possible parent found for cell at x, y (" << this_step.x;
                    ss << ", " << this_step.y << ") and wrap (" << this_step.xwrap << ", " << this_step.ywrap << ")";
                    ss << " at generation " << generation << " and with minimum distance of " << min_distance;
                    ss << ". Moving to nearest habitat cell." << endl;
#endif // DEBUG
                    disperseNearestHabitat(this_step);
                    fail = false;
                }

            }
        }
        else
        {
            while(fail)
            {
                angle = NR->direction();
                dist = NR->dispersal();
                density = landscape->runDispersal(dist,
                                                  angle,
                                                  this_step.x,
                                                  this_step.y,
                                                  this_step.xwrap,
                                                  this_step.ywrap,
                                                  fail,
                                                  *generation);
                // Discard the dispersal event a percentage of the time, based on the maximum value of the habitat map.
                // This is to correctly mimic less-dense cells having a lower likelihood of being the parent to the cell.
                if(!fail)
                {
                    fail = !checkEndPoint(density,
                                          this_step.x,
                                          this_step.y,
                                          this_step.xwrap,
                                          this_step.ywrap,
                                          startx,
                                          starty,
                                          startxwrap,
                                          startywrap);
                }
            }

        }
#ifdef DEBUG
        if(landscape->getVal(this_step.x, this_step.y, this_step.xwrap, this_step.ywrap, *generation) == 0 && !fail)
        {
            stringstream ss;
            ss << "x,y: " << this_step.x << "," << this_step.y;
            ss << " x,y wrap: " << this_step.xwrap << "," << this_step.ywrap << "Habitat cover: ";
            ss << landscape->getVal(this_step.x, this_step.y, this_step.xwrap, this_step.ywrap, *generation) << endl;
            writeLog(50, ss);
            throw FatalException("ERROR_MOVE_007: Dispersal attempted to non-habitat. Check dispersal function.");
        }
#endif
    }

    void DispersalCoordinator::disperseNearestHabitat(Step &this_step)
    {
        double end_x = this_step.x + 0.5;
        double end_y = this_step.y + 0.5;
        long end_x_wrap = 0;
        long end_y_wrap = 0;
        landscape->findNearestHabitatCell(this_step.x,
                                          this_step.y,
                                          this_step.xwrap,
                                          this_step.ywrap,
                                          end_x,
                                          end_y,
                                          *generation);
        landscape->fixGridCoordinates(end_x, end_y, end_x_wrap, end_y_wrap);
        end_x_wrap += this_step.xwrap;
        end_y_wrap += this_step.ywrap;
        if(!landscape->checkMap(end_x, end_y, end_x_wrap, end_y_wrap, *generation))
        {
            stringstream ss;
            ss << "Attempted nearest habitat cell is not habitat! Please report this bug." << endl;
            ss << "Nearby habitat cell at " << end_x << ", " << end_y << " (" << end_x_wrap << ", " << end_y_wrap;
            ss << ") does not contain habitat. Initial cell was ";
            ss << this_step.x << ", " << this_step.y << " (" << this_step.xwrap << ", " << this_step.ywrap;
            ss << ". Density of new cell was "
               << landscape->checkMap(end_x, end_y, end_x_wrap, end_y_wrap, *generation);
            double tmpx, tmpy;
            landscape->findNearestHabitatCell(this_step.x,
                                              this_step.y,
                                              this_step.xwrap,
                                              this_step.ywrap,
                                              tmpx,
                                              tmpy,
                                              *generation);
            ss << "Coords of nearest habitat :" << tmpx << ", " << tmpy << endl;
            ss << endl;
            throw FatalException(ss.str());
        }
        this_step.x = static_cast<long>(end_x);
        this_step.y = static_cast<long>(end_y);
        this_step.xwrap = end_x_wrap;
        this_step.ywrap = end_y_wrap;
    }

    void DispersalCoordinator::setEndPointFptr(const bool &restrict_self)
    {
        if(restrict_self)
        {
            if(reproduction_map->isNull())
            {
                checkEndPointFptr = &DispersalCoordinator::checkEndPointRestricted;
            }
            else
            {
                checkEndPointFptr = &DispersalCoordinator::checkEndPointDensityRestrictedReproduction;
            }
        }
        else
        {
            if(reproduction_map->isNull())
            {
                checkEndPointFptr = &DispersalCoordinator::checkEndPointDensity;
            }
            else
            {
                checkEndPointFptr = &DispersalCoordinator::checkEndPointDensityReproduction;
            }
        }
    }

    bool DispersalCoordinator::checkEndPoint(const unsigned long &density,
                                             long &x,
                                             long &y,
                                             long &xwrap,
                                             long &ywrap,
                                             const long &startx,
                                             const long &starty,
                                             const long &startxwrap,
                                             const long &startywrap)
    {
        return (this->*checkEndPointFptr)(density, x, y, xwrap, ywrap, startx, starty, startxwrap, startywrap);
    }

    bool DispersalCoordinator::checkEndPointDensity(const unsigned long &density,
                                                    long &x,
                                                    long &y,
                                                    long &xwrap,
                                                    long &ywrap,
                                                    const long &startx,
                                                    const long &starty,
                                                    const long &startxwrap,
                                                    const long &startywrap)
    {
        if((double(density) / double(landscape->getHabitatMax())) < NR->d01())
        {
            x = startx;
            y = starty;
            xwrap = startxwrap;
            ywrap = startywrap;
            return false;
        }
        return true;
    }

    bool DispersalCoordinator::checkEndPointRestricted(const unsigned long &density,
                                                       long &x,
                                                       long &y,
                                                       long &xwrap,
                                                       long &ywrap,
                                                       const long &startx,
                                                       const long &starty,
                                                       const long &startxwrap,
                                                       const long &startywrap)
    {
        if(startx == x && starty == y && startxwrap == xwrap && startywrap == ywrap)
        {
            return false;
        }
        return checkEndPointDensity(density, x, y, xwrap, ywrap, startx, starty, startxwrap, startywrap);
    }

    bool DispersalCoordinator::checkEndPointDensityReproduction(const unsigned long &density,
                                                                long &x,
                                                                long &y,
                                                                long &xwrap,
                                                                long &ywrap,
                                                                const long &startx,
                                                                const long &starty,
                                                                const long &startxwrap,
                                                                const long &startywrap)
    {
        if(checkEndPointDensity(density, x, y, xwrap, ywrap, startx, starty, startxwrap, startywrap))
        {
            if(!reproduction_map->actionOccurs(x, y, xwrap, ywrap))
            {
                x = startx;
                y = starty;
                xwrap = startxwrap;
                ywrap = startywrap;
                return false;
            }
            return true;
        }
        return false;

    }

    bool DispersalCoordinator::checkEndPointDensityRestrictedReproduction(const unsigned long &density,
                                                                          long &x,
                                                                          long &y,
                                                                          long &xwrap,
                                                                          long &ywrap,
                                                                          const long &startx,
                                                                          const long &starty,
                                                                          const long &startxwrap,
                                                                          const long &startywrap)
    {
        if(checkEndPointRestricted(density, x, y, xwrap, ywrap, startx, starty, startxwrap, startywrap))
        {
            if(!reproduction_map->actionOccurs(x, y, xwrap, ywrap))
            {
                x = startx;
                y = starty;
                xwrap = startxwrap;
                ywrap = startywrap;
                return false;
            }
            return true;
        }
        return false;

    }

    void DispersalCoordinator::disperse(Step &this_step)
    {
        (this->*doDispersal)(this_step);
    }

    double DispersalCoordinator::getSelfDispersalValue(const Cell &cell) const
    {
        if(!full_dispersal_map)
        {
            return 1.0;
        }
        unsigned long cell_index = calculateCellIndex(cell);
        if(cell_index >= raw_dispersal_prob_map.getCols())
        {
            stringstream ss;
            ss << "Index of " << cell_index << " for cell " << cell.x << ", " << cell.y
               << " is out of range of dispersal map with bounds " << raw_dispersal_prob_map.getCols() << ", "
               << raw_dispersal_prob_map.getRows() << endl;
            throw FatalException(ss.str());
        }
        return raw_dispersal_prob_map.getCopy(cell_index, cell_index);
    }

    double DispersalCoordinator::sumDispersalValues(const Cell &cell) const
    {
        if(!full_dispersal_map)
        {
            return raw_dispersal_prob_map.getCols();
        }
        unsigned long cell_index = calculateCellIndex(cell);
        double sum_total = 0.0;
        for(unsigned long x = 0; x < raw_dispersal_prob_map.getCols(); x++)
        {
            sum_total += raw_dispersal_prob_map.get(cell_index, x);
        }
        return sum_total;

    }

    void DispersalCoordinator::reimportRawDispersalMap()
    {
        if(raw_dispersal_prob_map.getCols() == 0 || raw_dispersal_prob_map.getRows() == 0)
        {
            raw_dispersal_prob_map.import(dispersal_prob_map.getFileName());
        }
    }

    void DispersalCoordinator::removeSelfDispersal()
    {
        reimportRawDispersalMap();
        Map<double> backup_dispersal_prob_map;
        backup_dispersal_prob_map = raw_dispersal_prob_map;
        for(unsigned long y = 0; y < raw_dispersal_prob_map.getRows(); y++)
        {
            raw_dispersal_prob_map.get(y, y) = 0.0;
        }
        dispersal_prob_map = raw_dispersal_prob_map;
        addDensity();
        addReproduction();
        fixDispersal();
        raw_dispersal_prob_map = backup_dispersal_prob_map;
#ifdef DEBUG
        validateNoSelfDispersalInDispersalMap();
#endif //DEBUG

    }

}