// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details
/**
 * @author Samuel Thompson
 * @file Landscape.cpp
 * @brief Contains the Landscape class implementation for easy referencing of the respective coarse and fine map within the
 * same coordinate system.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#define _USE_MATH_DEFINES

#include <cmath>
#include "Landscape.h"
#include "file_system.h"

uint32_t importToMapAndRound(string map_file, Map<uint32_t> &map_in, unsigned long map_x,
                             unsigned long map_y,
                             unsigned long scalar)
{
#ifndef SIZE_LIMIT
    if(map_x > 1000000 || map_y > 1000000)
    {
        throw runtime_error("Extremely large map sizes set for " + map_file + ": " + to_string(map_x) + ", " +
                            to_string(map_y) + "\n");
    }
#endif
    Map<float> temp_matrix;
    temp_matrix.setSize(map_y, map_x);
#ifdef DEBUG
    writeInfo("Calculating fine map");
#endif
    if(map_file == "null")
    {
        writeInfo("Setting null map.\n");
        for(unsigned long i = 0; i < map_y; i++)
        {
            for(unsigned long j = 0; j < map_x; j++)
            {
                temp_matrix[i][j] = 1.0;
            }
        }
    }
    else  // There is a map to read in.
    {
        temp_matrix.import(map_file);
        map_in.open(map_file);
        map_in.getRasterBand();
        map_in.getMetaData();
        map_in.close();
    }
#ifdef DEBUG
    writeInfo("import complete");
#endif
    uint32_t max_value = 0;
    map_in.setSize(temp_matrix.getRows(), temp_matrix.getCols());

    for(unsigned long i = 0; i < temp_matrix.getRows(); i++)
    {
        for(unsigned long j = 0; j < temp_matrix.getCols(); j++)
        {
            map_in[i][j] = (uint32_t) (max(round(temp_matrix[i][j] * (double) scalar), 0.0));
            if(map_in[i][j] > max_value)
            {
                max_value = map_in[i][j];
            }
        }
    }
    temp_matrix.close();
    return max_value;
}

unsigned long archimedesSpiralX(const double &centre_x, const double &centre_y, const double &radius,
                                const double &theta)
{
    return static_cast<unsigned long>(floor(radius * cos(theta) + centre_x));
}

unsigned long archimedesSpiralY(const double &centre_x, const double &centre_y, const double &radius,
                                const double &theta)
{
    return static_cast<unsigned long>(floor(radius * sin(theta) + centre_y));
}

void Landscape::setDims(shared_ptr<SimParameters> mapvarsin)
{
    if(!check_set_dim)  // checks to make sure it hasn't been run already.
    {
#ifdef DEBUG
        if(mapvarsin == nullptr)
        {
            throw FatalException("SimParameters pointer is nullptr. Please report this bug.");
        }
#endif // DEBUG
        mapvars = mapvarsin;
        mapvars->setHistorical(0);
        deme = mapvars->deme;
        x_dim = mapvars->grid_x_size;
        y_dim = mapvars->grid_y_size;
        scale = mapvars->coarse_map_scale;
        nUpdate = 0;
        check_set_dim = true;
        update_time = 0;
        updateMap(0);
        gen_since_historical = mapvars->gen_since_historical;
        if(gen_since_historical == 0)
        {
            gen_since_historical = 0.000000000000000001;
        }
        habitat_change_rate = mapvars->habitat_change_rate;
        landscape_type = mapvars->landscape_type;
    }
    else
    {
        writeError("ERROR_MAP_001: Dimensions have already been set");
    }
}

bool Landscape::checkMapExists()
{
    for(unsigned long i = 0; i < mapvars->configs.getSectionOptionsSize(); i++)
    {
        string tmppath = mapvars->configs[i].getOption("path");
        if(!doesExistNull(tmppath))
        {
            return false;
        }
    }
    return true;
}

void Landscape::calcFineMap()
{
    string fileinput = mapvars->fine_map_file;
    unsigned long mapxsize = mapvars->fine_map_x_size;
    unsigned long mapysize = mapvars->fine_map_y_size;
    if(!check_set_dim)  // checks that the dimensions have been set.
    {
        throw FatalException("ERROR_MAP_002: dimensions not set.");
    }
    // Note that the default "null" type is to have 100% forest cover in every cell.
    fine_max = importToMapAndRound(fileinput, fine_map, mapxsize, mapysize, deme);
}

void Landscape::calcHistoricalFineMap()
{
    string file_input = mapvars->historical_fine_map_file;
    unsigned long map_x_size = mapvars->fine_map_x_size;
    unsigned long map_y_size = mapvars->fine_map_y_size;
    if(!check_set_dim)  // checks that the dimensions have been set.
    {
        throw FatalException("ERROR_MAP_002: dimensions not set.");
    }
    has_historical = file_input != "none";
    historical_fine_max = 0;
    if(has_historical)
    {
        historical_fine_max = importToMapAndRound(file_input, historical_fine_map, map_x_size, map_y_size, deme);
    }
}

void Landscape::calcCoarseMap()
{
    string file_input = mapvars->coarse_map_file;
    unsigned long map_x_size = mapvars->coarse_map_x_size;
    unsigned long map_y_size = mapvars->coarse_map_y_size;
    if(!check_set_dim)  // checks that the dimensions have been set.
    {
        throw FatalException("ERROR_MAP_003: dimensions not set.");
    }
    has_coarse = file_input != "none";
    coarse_max = 0;
    if(has_coarse)
    {
        coarse_max = importToMapAndRound(file_input, coarse_map, map_x_size, map_y_size, deme);
    }
}

void Landscape::calcHistoricalCoarseMap()
{
    string file_input = mapvars->historical_coarse_map_file;
    unsigned long map_x_size = mapvars->coarse_map_x_size;
    unsigned long map_y_size = mapvars->coarse_map_y_size;
    if(!check_set_dim)  // checks that the dimensions have been set.
    {
        throw FatalException("ERROR_MAP_003: dimensions not set.");
    }
    historical_coarse_max = 0;
    if(has_coarse)
    {
        has_historical = file_input != "none";
        if(has_historical)
        {
            historical_coarse_max = importToMapAndRound(file_input, historical_coarse_map, map_x_size, map_y_size,
                                                        deme);
        }
    }
}

void Landscape::setTimeVars(double gen_since_historical_in, double habitat_change_rate_in)
{
    update_time = 0;
    gen_since_historical = gen_since_historical_in;
    habitat_change_rate = habitat_change_rate_in;
}

void Landscape::calcOffset()
{
    if(mapvars->times_file != "null")
    {
        mapvars->setHistorical(0);
    }
    if(fine_map.getCols() == 0 || fine_map.getRows() == 0)
    {
        throw FatalException("ERROR_MAP_004: fine map not set.");
    }
    if(coarse_map.getCols() == 0 || coarse_map.getRows() == 0)
    {
        if(has_coarse)
        {
            coarse_map.setSize(fine_map.getRows(), fine_map.getCols());
        }
    }
    if(checkAllDimensionsZero())
    {
        calculateOffsetsFromMaps();
    }
    else
    {
        calculateOffsetsFromParameters();
    }
    dispersal_relative_cost = mapvars->dispersal_relative_cost;
#ifdef DEBUG
    stringstream os;
    os << "\nfinex: " << fine_x_min << "," << fine_x_max << endl;
    os << "finey: " << fine_y_min << "," << fine_y_max << endl;
    os << "coarsex: " << coarse_x_min << "," << coarse_x_max << endl;
    os << "coarsey: " << coarse_y_min << "," << coarse_y_max << endl;
    os << "offsets: "
       << "(" << fine_x_offset << "," << fine_y_offset << ")(" << coarse_x_offset << "," << coarse_y_offset << ")"
       << endl;
    os << "historical fine file: " << historical_fine_map << endl;
    os << "historical coarse file: " << historical_coarse_map << endl;
    writeInfo(os.str());
#endif
    //		os << "fine variables: " << finexmin << "," << fine_x_max << endl;
    //		os << "coarse variabes: " << coarse_x_min << "," << coarse_x_max << endl;
    if(fine_x_min < coarse_x_min || fine_x_max > coarse_x_max || (fine_x_max - fine_x_min) < x_dim ||
       (fine_y_max - fine_y_min) < y_dim)
    {
        throw FatalException(
                "ERROR_MAP_006: FATAL - fine map extremes outside coarse map or sample grid larger than fine map");
    }
}

bool Landscape::checkAllDimensionsZero()
{
    return mapvars->fine_map_x_offset == 0 && mapvars->fine_map_y_offset == 0 && mapvars->coarse_map_x_offset == 0 &&
           mapvars->coarse_map_y_offset == 0 && mapvars->sample_x_offset == 0 && mapvars->sample_y_offset == 0 &&
           mapvars->fine_map_x_size == 0 && mapvars->fine_map_y_size == 0 && mapvars->coarse_map_x_size == 0 &&
           mapvars->coarse_map_y_size == 0;
}

void Landscape::calculateOffsetsFromMaps()
{
    long x_offset, y_offset;
    if(mapvars->sample_mask_file != "null" && mapvars->sample_mask_file != "none")
    {
        writeInfo("Calculating offsets from maps...\n");
        // Opens an empty map object for the sample mask file and then calculates the offsets.
        Map<uint32_t> tmp_sample_map;
        tmp_sample_map.open(mapvars->sample_mask_file);
        tmp_sample_map.getRasterBand();
        tmp_sample_map.getMetaData();
        tmp_sample_map.calculateOffset(fine_map, x_offset, y_offset);
        if(tmp_sample_map.roundedScale(fine_map) != 1)
        {
            writeInfo("Sample map resolution does not match fine map resolution.\n");
        }
        tmp_sample_map.close();
        if(x_offset < 0 || y_offset < 0)
        {
            stringstream ss;
            ss << "Fine map upper-left coordinates: " << fine_map.getUpperLeftX() << ", " << fine_map.getUpperLeftY();
            ss << endl << "Sample map upper-left coordinates: " << tmp_sample_map.getUpperLeftX() << ", ";
            ss << tmp_sample_map.getUpperLeftY() << endl;
            writeInfo(ss.str());
            ss.str("");
            ss << "Offsets of " << mapvars->fine_map_file << " from " << mapvars->sample_mask_file << " are negative (";
            ss << x_offset << ", " << y_offset << "): ";
            ss << "check map files are set correctly.\n" << endl;
            throw FatalException(ss.str());
        }
        mapvars->fine_map_x_offset = static_cast<unsigned long>(x_offset);
        mapvars->fine_map_y_offset = static_cast<unsigned long>(y_offset);
    }
    mapvars->coarse_map_x_size = coarse_map.getCols();
    mapvars->coarse_map_y_size = coarse_map.getRows();
    mapvars->fine_map_x_size = fine_map.getCols();
    mapvars->fine_map_y_size = fine_map.getRows();
    mapvars->sample_x_offset = 0;
    mapvars->sample_y_offset = 0;
    mapvars->sample_x_size = mapvars->fine_map_x_size;
    mapvars->sample_y_size = mapvars->fine_map_y_size;
    mapvars->grid_x_size = mapvars->fine_map_x_size;
    mapvars->grid_y_size = mapvars->fine_map_y_size;
    x_dim = mapvars->grid_x_size;
    y_dim = mapvars->grid_y_size;
    fine_map.calculateOffset(coarse_map, x_offset, y_offset);
    mapvars->coarse_map_x_offset = static_cast<unsigned long>(x_offset);
    mapvars->coarse_map_y_offset = static_cast<unsigned long>(y_offset);
    mapvars->coarse_map_scale = fine_map.roundedScale(coarse_map);
    scale = mapvars->coarse_map_scale;
    if(x_offset < 0 || y_offset < 0)
    {
        stringstream ss;
        ss << "Fine map upper-left coordinates: " << fine_map.getUpperLeftX() << ", " << fine_map.getUpperLeftY();
        ss << endl << "Coarse map upper-left coordinates: " << coarse_map.getUpperLeftX() << ", ";
        ss << fine_map.getUpperLeftY();
        writeInfo(ss.str());
        ss.str("");
        ss << "Offsets of " << mapvars->coarse_map_file << " from " << mapvars->fine_map_file << " are negative (";
        ss << x_offset << ", " << y_offset << "): ";
        ss << "check map files are set correctly." << endl;
        throw FatalException(ss.str());
    }
    stringstream ss;
    ss << "Dimensions detected as: " << endl;
    ss << "Fine map" << endl;
    ss << "-dimensions: " << fine_map.getCols() << ", " << fine_map.getRows() << endl;
    ss << "-offsets: " << mapvars->fine_map_x_offset << ", " << mapvars->fine_map_y_offset << endl;
    ss << "Coarse map" << endl;
    ss << "-dimensions: " << coarse_map.getCols() << ", " << coarse_map.getRows() << endl;
    ss << "-offsets: " << mapvars->coarse_map_x_offset << ", " << mapvars->coarse_map_y_offset << endl;
    ss << "-scale: " << mapvars->coarse_map_scale << endl;
    writeInfo(ss.str());
    calculateOffsetsFromParameters();
}

void Landscape::calculateOffsetsFromParameters()
{
    fine_x_offset = mapvars->fine_map_x_offset + mapvars->sample_x_offset;
    fine_y_offset = mapvars->fine_map_y_offset + mapvars->sample_y_offset;
    coarse_x_offset = mapvars->coarse_map_x_offset;
    coarse_y_offset = mapvars->coarse_map_y_offset;
    scale = mapvars->coarse_map_scale;
    // this is the location of the top left (or north west) corner of the respective map
    // and the x and y distance from the top left of the grid object that contains the initial lineages.
    fine_x_min = -fine_x_offset;
    fine_y_min = -fine_y_offset;
    fine_x_max = fine_x_min + (fine_map.getCols());
    fine_y_max = fine_y_min + (fine_map.getRows());
    if(has_coarse) // Check if there is a coarse map
    {
        coarse_x_min = -coarse_x_offset - fine_x_offset;
        coarse_y_min = -coarse_y_offset - fine_y_offset;
        coarse_x_max = coarse_x_min + scale * (coarse_map.getCols());
        coarse_y_max = coarse_y_min + scale * (coarse_map.getRows());
    }
    else // Just set the offsets to the same as the fine map
    {
        coarse_x_min = fine_x_min;
        coarse_y_min = fine_y_min;
        coarse_x_max = fine_x_max;
        coarse_y_max = fine_y_max;
        scale = 1;
    }
}

void Landscape::validateMaps()
{
    stringstream os;
    os << "\rValidating maps..." << flush;
#ifdef historical_mode
    double dTotal = fine_map.getCols() + coarse_map.getCols();
    unsigned long iCounter = 0;
#endif // historical_mode
    if(has_historical)
    {
        if(fine_map.getCols() == historical_fine_map.getCols() && fine_map.getRows() == historical_fine_map.getRows() &&
           coarse_map.getCols() == historical_coarse_map.getCols() && coarse_map.getRows() ==
                                                                      historical_coarse_map.getRows())
        {
            os << "\rValidating maps...map sizes okay" << flush;
            writeInfo(os.str());
        }
        else
        {
            throw FatalException(
                    "ERROR_MAP_009: Landscape validation failed - modern and historical maps are not the same dimensions.");
        }
#ifdef historical_mode
        for(unsigned long i = 0; i < fine_map.getCols(); i++)
        {
            for(unsigned long j = 0; j < fine_map.getRows(); j++)
            {
                if(fine_map[j][i] > historical_fine_map[j][i])
                {
#ifdef DEBUG
                    stringstream ss;
                    ss << "fine map: " << fine_map[j][i] << " historical map: " << historical_fine_map[j][i];
                    ss << " x,y: " << i << "," << j << endl;
                    writeLog(50, ss);
#endif //DEBUG
                    throw FatalException("ERROR_MAP_007: Landscape validation failed - fine map value larger "
                                              "than historical fine map value.");
                }
            }
            double dPercentComplete = 100 * ((double)(i + iCounter) / dTotal);
            if(i % 1000 == 0)
            {
                os.str("");
                os << "\rValidating maps..." << dPercentComplete << "%                " << flush;
                writeInfo(os.str());
            }
        }
#endif // historical mode
    }
#ifdef historical_mode
    iCounter = fine_map.getCols();
    stringstream ss;
    if(has_historical)
    {
        for(unsigned long i = 0; i < coarse_map.getCols(); i++)
        {
            for(unsigned long j = 0; j < coarse_map.getRows(); j++)
            {
                if(coarse_map[j][i] > historical_coarse_map[j][i])
                {
                    ss << "coarse map: " << coarse_map[j][i] << " historical map: " << historical_coarse_map[j][i];
                    ss << " coarse map x+1: " << coarse_map[j][i + 1]
                       << " historical map: " << historical_coarse_map[j][i + 1];
                    ss << " x,y: " << i << "," << j;
                    writeLog(50, ss);
                    throw FatalException("ERROR_MAP_008: Landscape validation failed - coarse map value larger "
                                         "than historical coarse map value.");
                }
            }
            double dPercentComplete = 100 * ((double) (i + iCounter) / dTotal);
            if(i % 1000 == 0)
            {
                os.str("");
                os << "\rValidating maps..." << dPercentComplete << "%                " << flush;
                writeInfo(os.str());
            }
        }

    }
#endif // historical_mode
    os.str("");
    os << "\rValidating maps complete                                       " << endl;
    writeInfo(os.str());
}

bool Landscape::updateMap(double generation)
{
    // only update the map if the historical state has not been reached.
    if(!mapvars->is_historical && has_historical)
    {
        if(mapvars->gen_since_historical < generation)
        {
            // Only update the map if the maps have actually changed
            if(mapvars->setHistorical(nUpdate + 1))
            {
                stringstream ss;
                ss << "\nUpdating historical maps at " << generation << "...\n";
                writeInfo(ss.str());
                fine_max = historical_fine_max;
                fine_map = historical_fine_map;
                coarse_max = historical_coarse_max;
                coarse_map = historical_coarse_map;
                doUpdate();
                return true;
            }
        }
    }
    return false;
}

void Landscape::doUpdate()
{
    nUpdate++;
    // historical_fine_map = mapvars->historical_fine_map_file;
    // historical_coarse_map = mapvars->historical_coarse_map_file;
    current_map_time = gen_since_historical;
    gen_since_historical = mapvars->gen_since_historical;
    if(gen_since_historical == 0)
    {
        gen_since_historical = 0.000000000000000001;
    }
    habitat_change_rate = mapvars->habitat_change_rate;
    calcHistoricalFineMap();
    calcHistoricalCoarseMap();
    if(has_historical)
    {
        is_historical = mapvars->is_historical;
    }
    recalculateHabitatMax();
}

void Landscape::resetHistorical()
{
    nUpdate = 0;
    doUpdate();
}

void Landscape::setLandscape(string landscape_type)
{
    infinite_boundaries = true;
    if(landscape_type == "infinite")
    {
        writeInfo("Setting infinite landscape.\n");
        getValFunc = &Landscape::getValInfinite;
    }
    else if(landscape_type == "tiled_coarse")
    {
        writeInfo("Setting tiled coarse infinite landscape.\n");
        getValFunc = &Landscape::getValCoarseTiled;
    }
    else if(landscape_type == "tiled_fine")
    {
        writeInfo("Setting tiled fine infinite landscape.\n");
        getValFunc = &Landscape::getValFineTiled;
    }
    else if(landscape_type == "closed")
    {
        infinite_boundaries = false;
        getValFunc = &Landscape::getValFinite;
    }
    else
    {
        throw FatalException("Provided landscape type is not a valid option: " + landscape_type);
    }
}

unsigned long Landscape::getVal(const double &x, const double &y,
                                const long &xwrap, const long &ywrap, const double &current_generation)
{
    return (this->*getValFunc)(x, y, xwrap, ywrap, current_generation);
}

unsigned long Landscape::getValInfinite(
        const double &x, const double &y, const long &xwrap, const long &ywrap, const double &current_generation)
{
    double xval, yval;
    xval = x + (x_dim * xwrap);  //
    yval = y + (y_dim * ywrap);
    //		// return 0 if the requested coordinate is completely outside the map
    if(xval < coarse_x_min || xval >= coarse_x_max || yval < coarse_y_min || yval >= coarse_y_max)
    {
        return deme;
    }
    return getValFinite(x, y, xwrap, ywrap, current_generation);
}

unsigned long Landscape::getValCoarseTiled(
        const double &x, const double &y, const long &xwrap, const long &ywrap, const double &current_generation)
{
    double newx = fmod(x + (xwrap * x_dim) + fine_x_offset + coarse_x_offset, coarse_map.getCols());
    double newy = fmod(y + (ywrap * y_dim) + fine_x_offset + coarse_x_offset, coarse_map.getRows());
    if(newx < 0)
    {
        newx += coarse_map.getCols();
    }
    if(newy < 0)
    {
        newy += coarse_map.getRows();
    }
    return getValCoarse(newx, newy, current_generation);
}

unsigned long Landscape::getValFineTiled(
        const double &x, const double &y, const long &xwrap, const long &ywrap, const double &current_generation)
{

    double newx = fmod(x + (xwrap * x_dim) + fine_x_offset, fine_map.getCols());
    double newy = fmod(y + (ywrap * y_dim) + fine_y_offset, fine_map.getRows());
    // Now adjust for incorrect wrapping behaviour of fmod
    if(newx < 0)
    {
        newx += fine_map.getCols();
    }
    if(newy < 0)
    {
        newy += fine_map.getRows();
    }
#ifdef DEBUG
    if(newx >= fine_map.getCols() || newx < 0 || newy >= fine_map.getRows() || newy < 0)
    {
        stringstream ss;
        ss << "Fine map indexing out of range of fine map." << endl;
        ss << "x, y: " << newx << ", " << newy << endl;
        ss << "cols, rows: " << fine_map.getCols() << ", " << fine_map.getRows() << endl;
        throw out_of_range(ss.str());
    }
#endif
    return getValFine(newx, newy, current_generation);
}

unsigned long Landscape::getValCoarse(const double &xval, const double &yval, const double &current_generation)
{
    unsigned long retval = 0;
    if(has_historical)
    {
        if(is_historical || historical_coarse_map[yval][xval] == coarse_map[yval][xval])
        {
            return historical_coarse_map[yval][xval];
        }
        else
        {
            double currentTime = current_generation - current_map_time;
            retval = (unsigned long) floor(coarse_map[yval][xval] +
                                           (habitat_change_rate *
                                            ((historical_coarse_map[yval][xval] - coarse_map[yval][xval]) /
                                             (gen_since_historical - current_map_time)) * currentTime));
        }
    }
    else
    {
        return coarse_map[yval][xval];
    }
#ifdef historical_mode
    if(retval > historical_coarse_map[yval][xval])
        {
            string ec =
                "Returned value greater than historical value. Check file input. (or disable this error before "
                "compilation.\n";
            ec += "historical value: " + to_string((long long)historical_coarse_map[yval][xval]) +
                  " returned value: " + to_string((long long)retval);
            throw FatalException(ec);
        }
// Note that debug mode will throw an exception if the returned value is less than the historical state

#endif
    return retval;
}

unsigned long Landscape::getValFine(const double &xval, const double &yval, const double &current_generation)
{
    unsigned long retval = 0;
    if(has_historical)
    {
        if(is_historical || historical_fine_map[yval][xval] == fine_map[yval][xval])
        {
            retval = historical_fine_map[yval][xval];
        }
        else
        {
            double currentTime = current_generation - current_map_time;
#ifdef historical_mode
            retval = (unsigned long)floor(fine_map[yval][xval] +
                                           (habitat_change_rate * ((historical_fine_map[yval][xval] - fine_map[yval][xval]) /
                                                   (gen_since_historical-current_map_time)) * currentTime));
#else
            retval = (unsigned long) floor(fine_map[yval][xval] +
                                           (habitat_change_rate *
                                            ((static_cast<double>(historical_fine_map[yval][xval]) -
                                              static_cast<double>(fine_map[yval][xval])) /
                                             (gen_since_historical - current_map_time)) * currentTime));
#endif
        }
    }
    else
    {
        return fine_map[yval][xval];
    }
// os <<fine_map[yval][xval] << "-"<< retval << endl;
// Note that debug mode will throw an exception if the returned value is less than the historical state
#ifdef historical_mode
    if(has_historical)
    {
        if(retval > historical_fine_map[yval][xval])
        {
            throw FatalException("Returned value greater than historical value. Check file input. (or disable this "
                                      "error before compilation.");
        }
    }
#endif
    return retval;
}

unsigned long Landscape::getValFinite(
        const double &x, const double &y, const long &xwrap, const long &ywrap, const double &current_generation)
{

    double xval, yval;
    xval = x + (x_dim * xwrap);  //
    yval = y + (y_dim * ywrap);
    //		// return 0 if the requested coordinate is completely outside the map
    if(xval < coarse_x_min || xval >= coarse_x_max || yval < coarse_y_min || yval >= coarse_y_max)
    {
        return 0;
    }
    if((xval < fine_x_min || xval >= fine_x_max || yval < fine_y_min ||
        yval >= fine_y_max) && has_coarse)  // check if the coordinate comes from the coarse resolution map.
    {
        // take in to account the fine map offsetting
        xval += fine_x_offset;
        yval += fine_y_offset;
        // take in to account the coarse map offsetting and the increased scale of the larger map.
        xval = floor((xval + coarse_x_offset) / scale);
        yval = floor((yval + coarse_y_offset) / scale);
        return getValCoarse(xval, yval, current_generation);
    }
    // take in to account the fine map offsetting
    // this is done twice to avoid having all the comparisons involve additions.
    xval += fine_x_offset;
    yval += fine_y_offset;
    return getValFine(xval, yval, current_generation);

}

unsigned long Landscape::convertSampleXToFineX(const unsigned long &x, const long &xwrap)
{
    return x + fine_x_offset + (xwrap * x_dim);
}

unsigned long Landscape::convertSampleYToFineY(const unsigned long &y, const long &ywrap)
{
    return y + fine_y_offset + (ywrap * y_dim);
}

void Landscape::convertFineToSample(long &x, long &xwrap, long &y, long &ywrap)
{
    auto tmpx = double(x - fine_x_offset);
    auto tmpy = double(y - fine_y_offset);
    fixGridCoordinates(tmpx, tmpy, xwrap, ywrap);
    x = static_cast<long>(floor(tmpx));
    y = static_cast<long>(floor(tmpy));
}

unsigned long Landscape::getInitialCount(double dSample, DataMask &samplemask)
{
    unsigned long toret;
    toret = 0;
    long x, y;
    long xwrap, ywrap;
    unsigned long max_x, max_y;
    if(samplemask.isNull())
    {
        max_x = fine_map.getCols();
        max_y = fine_map.getRows();
    }
    else
    {
        max_x = samplemask.sample_mask.getCols();
        max_y = samplemask.sample_mask.getRows();
    }
    for(unsigned long i = 0; i < max_x; i++)
    {
        for(unsigned long j = 0; j < max_y; j++)
        {
            x = i;
            y = j;
            xwrap = 0;
            ywrap = 0;
            samplemask.recalculateCoordinates(x, y, xwrap, ywrap);
            toret += (unsigned long) (max(floor(dSample * (getVal(x, y, xwrap, ywrap, 0)) *
                                                samplemask.getExactValue(x, y, xwrap, ywrap)), 0.0));
        }
    }
    return toret;
}

shared_ptr<SimParameters> Landscape::getSimParameters()
{
    if(!mapvars)
    {
        throw FatalException("Simulation current_metacommunity_parameters have not yet been set.");
    }
    return mapvars;
}

bool
Landscape::checkMap(const double &x, const double &y, const long &xwrap, const long &ywrap, const double generation)
{
    return getVal(x, y, xwrap, ywrap, generation) != 0;
}

bool Landscape::isOnFine(const double &x, const double &y, const long &xwrap, const long &ywrap)
{
    double tmpx, tmpy;
    tmpx = x + xwrap * x_dim;
    tmpy = y + ywrap * y_dim;
    return !(tmpx < fine_x_min || tmpx >= fine_x_max || tmpy < fine_y_min || tmpy >= fine_y_max);
}

bool Landscape::isOnCoarse(const double &x, const double &y, const long &xwrap, const long &ywrap)
{
    double xval, yval;
    xval = x + xwrap * x_dim;
    yval = y + ywrap * y_dim;
    return xval < coarse_x_min || xval >= coarse_x_max || yval < coarse_y_min || yval >= coarse_y_max;
}

bool Landscape::isOnMap(const double &x, const double &y, const long &xwrap, const long &ywrap)
{
    if(infinite_boundaries)
    {
        return true;
    }
    if(has_coarse && isOnCoarse(x, y, xwrap, ywrap))
    {
        return true;
    }
    return isOnFine(x, y, xwrap, ywrap);
}

void Landscape::fixGridCoordinates(double &x, double &y, long &xwrap, long &ywrap)
{
    xwrap += floor(x / x_dim);
    ywrap += floor(y / y_dim);
    x = x - xwrap * x_dim;
    y = y - ywrap * y_dim;
}

unsigned long Landscape::runDispersal(const double &dist,
                                      const double &angle,
                                      long &startx,
                                      long &starty,
                                      long &startxwrap,
                                      long &startywrap,
                                      bool &disp_comp,
                                      const double &generation)
{
// Checks that the start point is not out of matrix - this might have to be disabled to ensure that when updating the
// map, it doesn't cause problems.
#ifdef historical_mode
    if(!checkMap(startx, starty, startxwrap, startywrap, generation))
    {
        disp_comp = true;
        return;
    }
#endif

    // Different calculations for each quadrant to ensure that the dispersal reads the probabilities correctly.
    double newx, newy;
    newx = startx + (x_dim * startxwrap) + 0.5;
    newy = starty + (y_dim * startywrap) + 0.5;
    if(dispersal_relative_cost == 1)
    {
        // then nothing complicated is required and we can jump straight to the final point.
        newx += dist * cos(angle);
        newy += dist * sin(angle);
    }
    else  // we need to see which deforested patches we pass over
    {
        long boost;
        boost = 1;
        double cur_dist, tot_dist, l;
        cur_dist = 0;
        tot_dist = 0;
        // Four different calculations for the different quadrants.
        if(angle > 7 * M_PI_4 || angle <= M_PI_4)
        {
            // Continue while the dist travelled is less than the dist energy
            while(cur_dist < dist)
            {
                // Check if the starting position of the loop is in the fine map or not.
                if(isOnFine(newx, newy, 0, 0))
                {
                    // Keep the standard movement rate
                    boost = 1;
                }
                else
                {
                    // Accellerate the travel speed if the point is outside the fine grid.
                    // Note this means that lineages travelling from outside the fine grid to within the
                    // fine grid may
                    // see 1 grid's worth of approximation, rather than exact values.
                    // This is an acceptable approximation!
                    boost = deme;
                }

                // Add the value to the new x and y values.
                newx = newx + boost;
                newy = newy + boost * tan(angle);
                // Check if the new point is within forest.
                if(checkMap(newx, newy, 0, 0, generation))
                {
                    l = 1;
                }
                else
                {
                    l = dispersal_relative_cost;
                }
                // Move forward different dists based on the difficulty of moving through forest.
                cur_dist = cur_dist + l * boost * (1 / cos(angle));
                tot_dist = tot_dist + boost * (1 / cos(angle));
            }
        }
        else if(angle > 3 * M_PI_4 && angle <= 5 * M_PI_4)
        {
            while(cur_dist < dist)
            {
                if(isOnFine(newx, newy, 0, 0))
                {
                    boost = 1;
                }
                else
                {
                    boost = deme;
                }
                // Add the change to the new x and y values.
                newx = newx - boost;
                newy = newy + boost * tan(M_PI - angle);
                if(checkMap(newx, newy, 0, 0, generation))
                {
                    l = 1;
                }
                else
                {
                    l = dispersal_relative_cost;
                }
                cur_dist = cur_dist + boost * l * (1 / cos(M_PI - angle));
                tot_dist = tot_dist + boost * (1 / cos(M_PI - angle));
            }
        }
        else if(angle > M_PI_4 && angle <= 3 * M_PI_4)
        {
            while(cur_dist < dist)
            {
                if(isOnFine(newx, newy, 0, 0))
                {
                    boost = 1;
                }
                else
                {
                    boost = deme;
                }
                // Add the change to the new x and y values.
                newx = newx + boost * tan(angle - M_PI_2);
                newy = newy + boost;
                if(checkMap(newx, newy, 0, 0, generation))
                {
                    l = 1;
                }
                else
                {
                    l = dispersal_relative_cost;
                }
                cur_dist = cur_dist + l * boost / cos(angle - M_PI_2);
                tot_dist = tot_dist + boost / cos(angle - M_PI_2);
            }
        }
        else if(angle > 5 * M_PI_4 && angle <= 7 * M_PI_4)
        {
            //				os << "...ang4..." <<  flush;
            while(cur_dist < dist)
            {
                if(isOnFine(newx, newy, 0, 0))
                {
                    boost = 1;
                }
                else
                {
                    boost = deme;
                }
                newx = newx + boost * tan(3 * M_PI_2 - angle);
                newy = newy - boost;
                if(checkMap(newx, newy, 0, 0, generation))
                {
                    l = 1;
                }
                else
                {
                    l = dispersal_relative_cost;
                }
                cur_dist = cur_dist + l * boost / cos(3 * M_PI_2 - angle);
                tot_dist = tot_dist + boost / cos(3 * M_PI_2 - angle);
            }
        }
        // Move the point back to get the exact placement
        if(checkMap(newx, newy, 0, 0, generation))
        {
            tot_dist = tot_dist - min(cur_dist - dist, (double(boost) - 0.001));
        }
        else
        {
            disp_comp = true;
        }
        newx = startx + 0.5 + tot_dist * cos(angle);
        newy = starty + 0.5 + tot_dist * sin(angle);
    }
    unsigned long ret = getVal(newx, newy, 0, 0, generation);
    if(ret > 0)
    {
        long newxwrap, newywrap;
        newxwrap = 0;
        newywrap = 0;
        fixGridCoordinates(newx, newy, newxwrap, newywrap);
#ifdef DEBUG
        if(!checkMap(newx, newy, newxwrap, newywrap, generation))
        {
            throw FatalException(string(
                    "ERROR_MOVE_007: Dispersal attempted to non-forest. Check dispersal function. Forest cover: " +
                    to_string((long long) getVal(newx, newy, newxwrap, newywrap, generation))));
        }
#endif
        startx = newx;
        starty = newy;
        startxwrap = newxwrap;
        startywrap = newywrap;
        disp_comp = false;
    }
    return ret;
};

double Landscape::distanceToNearestHabitat(const long &start_x, const long &start_y, const long &start_x_wrap,
                                           const long &start_y_wrap, const double &generation)
{
    long end_x = start_x;
    long end_y = start_y;
    double theta = 0;
    double radius = 1.0;
    if(!getVal(end_x, end_y, start_x_wrap, start_y_wrap, generation))
    {
        while(true)
        {
            theta += 0.5 * M_PI / (2.0 * max(radius, 1.0));
            radius = theta / (2 * M_PI);
            end_x = archimedesSpiralX(start_x, start_y, radius, theta);
            end_y = archimedesSpiralY(start_x, start_y, radius, theta);

            // Double check that the distance is not greater than the map size
            // This acts as a fail-safe in case someone presents a historical map with no habitat cells on
            if(!isOnMap(end_x, end_y, start_x_wrap, start_y_wrap))
            {
                if(radius > fine_map.getCols() && radius > fine_map.getRows() &&
                   radius > coarse_map.getCols() * scale && radius > coarse_map.getRows() * scale)
                    throw FatalException(
                            "Could not find a habitat cell for parent. Check that your map files always have a "
                            "place for lineages to disperse from.");
            }
            else
            {
                if(getVal(end_x, end_y, start_x_wrap, start_y_wrap, generation))
                {
                    break;
                }
            }
        }
    }
    return pow(pow((start_x - end_x), 2) + pow((start_y - end_y), 2), 0.5);
}

void Landscape::clearMap()
{
    current_map_time = 0;
    check_set_dim = false;
    is_historical = false;
}

string Landscape::printVars()
{
    stringstream os;
    os << "fine x limits: " << fine_x_min << " , " << fine_x_max << endl;
    os << "fine y limits: " << fine_y_min << " , " << fine_y_max << endl;
    os << "fine map offset: " << fine_x_offset << " , " << fine_y_offset << endl;
    os << "coarse x limits: " << coarse_x_min << " , " << coarse_x_max << endl;
    os << "coarse y limits: " << coarse_y_min << " , " << coarse_y_max << endl;
    os << "x,y dims: " << x_dim << " , " << y_dim << endl;
    return os.str();
}

unsigned long Landscape::getHabitatMax()
{
    return habitat_max;
}

bool Landscape::hasHistorical()
{
    return has_historical;
}

void Landscape::recalculateHabitatMax()
{
    habitat_max = 0;
    if(is_historical && has_historical)
    {
        if(habitat_max < historical_fine_max)
        {
            habitat_max = historical_fine_max;
        }
        if(habitat_max < historical_coarse_max)
        {
            habitat_max = historical_coarse_max;
        }
    }
    else
    {
        if(habitat_max < fine_max)
        {
            habitat_max = fine_max;
        }
        if(habitat_max < coarse_max)
        {
            habitat_max = coarse_max;
        }
        if(habitat_max < historical_fine_max)
        {
            habitat_max = historical_fine_max;
        }
        if(habitat_max < historical_coarse_max)
        {
            habitat_max = historical_coarse_max;
        }
    }
#ifdef DEBUG
    if(habitat_max > 10000)
    {
        stringstream ss;
        writeLog(10, "habitat_max may be unreasonably large: " + to_string(habitat_max));
        ss << "fine, coarse, pfine, pcoarse: " << fine_max << ", " << coarse_max;
        ss << ", " << historical_fine_max << ", " << historical_coarse_max << endl;
    }
#endif
}

