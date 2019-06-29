
/**
 * @author Samuel Thompson
 * @file DataMask.cpp
 * @brief  Contains the DataMask class for describing the spatial sampling pattern on a landscape.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include "DataMask.h"
#include "Landscape.h"
#include "Logging.h"

DataMask::DataMask() : inputfile(""), isNullSample(true), isGridOffset(false), x_offset(0), y_offset(0), x_dim(0),
                       y_dim(0), mask_x_dim(0), mask_y_dim(0), getProportionfptr(nullptr), sample_mask(),
                       sample_mask_exact()
{
    getProportionfptr = &DataMask::getBoolProportion;
}

bool DataMask::isNull()
{
    return isNullSample;
}

void DataMask::setup(const shared_ptr<SimParameters> sim_parameters)
{
#ifdef DEBUG
    if((sim_parameters->grid_x_size > sim_parameters->sample_x_size ||
        sim_parameters->grid_y_size > sim_parameters->sample_y_size) && !isNullSample)
    {
        writeLog(50, "Grid size: " + to_string(sim_parameters->grid_x_size) + ", " +
                     to_string(sim_parameters->grid_y_size));
        writeLog(50, "Sample mask size: " + to_string(sim_parameters->sample_x_size) + ", " +
                     to_string(sim_parameters->sample_y_size));
        throw FatalException("Datamask dimensions do not make sense");
    }
#endif // DEBUG
    inputfile = sim_parameters->sample_mask_file;
    x_dim = sim_parameters->grid_x_size;
    y_dim = sim_parameters->grid_y_size;
    mask_x_dim = sim_parameters->sample_x_size;
    mask_y_dim = sim_parameters->sample_y_size;
    x_offset = sim_parameters->sample_x_offset;
    y_offset = sim_parameters->sample_y_offset;
    if(x_dim != mask_x_dim || y_dim != mask_y_dim)
    {
        isGridOffset = true;
    }
    else
    {
        isGridOffset = false;
    }
}

bool DataMask::checkCanUseDefault(const shared_ptr<SimParameters> sim_parameters)
{
    if(sim_parameters->sample_mask_file == "null")
    {
        if(sim_parameters->fine_map_x_size == sim_parameters->sample_x_size &&
           sim_parameters->fine_map_y_size == sim_parameters->sample_y_size &&
           sim_parameters->fine_map_x_offset == 0 && sim_parameters->fine_map_y_offset == 0)
        {
            isNullSample = true;
        }
        else
        {
            isNullSample = false;
        }
    }
    else
    {
        isNullSample = sim_parameters->sample_mask_file == "none";
    }
    return isNullSample;
}

void DataMask::importBooleanMask(unsigned long xdim, unsigned long ydim, unsigned long mask_xdim,
                                 unsigned long mask_ydim,
                                 unsigned long xoffset, unsigned long yoffset, string inputfile_in)
{
    shared_ptr<SimParameters> tmp_sim_parameters = make_shared<SimParameters>();
    tmp_sim_parameters->sample_mask_file = inputfile_in;
    tmp_sim_parameters->grid_x_size = xdim;
    tmp_sim_parameters->grid_y_size = ydim;
    tmp_sim_parameters->sample_x_size = mask_xdim;
    tmp_sim_parameters->sample_y_size = mask_ydim;
    tmp_sim_parameters->sample_x_offset = xoffset;
    tmp_sim_parameters->sample_y_offset = yoffset;
    setup(tmp_sim_parameters);
    isNullSample = inputfile_in == "null" || inputfile_in == "none";
    if(!isNullSample)
    {
        doImport();
    }
}

void DataMask::doImport()
{
    sample_mask.setSize(mask_y_dim, mask_x_dim);
    sample_mask.import(inputfile);
    sample_mask.close();
    completeBoolImport();
}

void DataMask::completeBoolImport()
{
    mask_x_dim = sample_mask.getCols();
    mask_y_dim = sample_mask.getRows();
    getProportionfptr = &DataMask::getBoolProportion;
}

void DataMask::setupNull(const shared_ptr<SimParameters> mapvarin)
{
    sample_mask.setSize(mapvarin->fine_map_y_size, mapvarin->fine_map_x_size);
    for(unsigned long i = 0; i < sample_mask.getRows(); i++)
    {
        for(unsigned long j = 0; j < sample_mask.getCols(); j++)
        {
            sample_mask.setValue(i, j, i + y_offset < mask_y_dim && j + x_offset < mask_x_dim);
        }
    }
    completeBoolImport();
}

void DataMask::importSampleMask(const shared_ptr<SimParameters> mapvarin)
{
    setup(mapvarin);
    if(!checkCanUseDefault(mapvarin))
    {
        if(inputfile == "null")
        {
            setupNull(mapvarin);
        }
        else if(mapvarin->uses_spatial_sampling)
        {
#ifdef DEBUG
            writeLog(10, "Using spatial sampling.");
            writeLog(10, "Mask dimensions: " + to_string(mask_x_dim) + ", " + to_string(mask_y_dim));
#endif // DEBUG
            sample_mask_exact.setSize(mask_y_dim, mask_x_dim);
            sample_mask_exact.import(inputfile);
            sample_mask_exact.close();
            mask_x_dim = sample_mask_exact.getCols();
            mask_y_dim = sample_mask_exact.getRows();
            getProportionfptr = &DataMask::getSampleProportion;
        }
        else
        {
            doImport();
        }
    }
    else
    {
        if(mapvarin->uses_spatial_sampling)
        {
            // This could perhaps be a warning, but I'd prefer to have the warning/prohibit potential in python
            // and throw a full exception here.
            throw FatalException("Cannot use a spatial sampling routine when the map file is null.");
        }
        getProportionfptr = &DataMask::getNullProportion;
    }
}

bool DataMask::getVal(const long &x, const long &y, const long &xwrap, const long &ywrap) const
{
    long xval = x + (xwrap * x_dim) + x_offset;
    long yval = y + (ywrap * y_dim) + y_offset;
    if(isNullSample)
    {
        return true;
    }
#ifdef DEBUG
    if(xval < 0 || xval >= (long) mask_x_dim || yval < 0 || yval >= (long) mask_y_dim)
    {
        stringstream ss;
        ss << "Get value on samplemask requested for non index." << endl;
        ss << "x, y: " << x << ", " << y << endl;
        ss << "dimensions x,y: " << mask_x_dim << ", " << mask_y_dim << endl;
        ss << "x, y wrap: " << xwrap << ", " << ywrap << endl;
        ss << "xval, yval: " << xval << ", " << yval << endl;
        ss << "offsets x, y: " << x_offset << ", " << y_offset << endl;
        writeLog(50, ss);
        ss.str("Get value on samplemask requested for non index.");
        throw out_of_range(ss.str());
    }
#endif
    return sample_mask.getCopy(yval, xval);
}

double DataMask::getNullProportion(const long &x, const long &y, const long &xwrap, const long &ywrap) const
{
    return 1.0;
}

double DataMask::getBoolProportion(const long &x, const long &y, const long &xwrap, const long &ywrap) const
{

    if(getVal(x, y, xwrap, ywrap))
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

double DataMask::getSampleProportion(const long &x, const long &y, const long &xwrap, const long &ywrap) const
{
#ifdef DEBUG
    if(isNullSample || sample_mask_exact.getCols() == 0)
    {
        throw out_of_range("Cannot get the exact value from a samplemask if we are using a null mask, or the "
                           "exact samplemask has not been properly imported.");
    }
#endif // DEBUG
    long xval = x + (xwrap * x_dim) + x_offset;
    long yval = y + (ywrap * y_dim) + y_offset;
    return sample_mask_exact.getCopy(yval, xval);
}

double DataMask::getExactValue(const long &x, const long &y, const long &xwrap, const long &ywrap) const
{
    return (this->*getProportionfptr)(x, y, xwrap, ywrap);
}

void DataMask::convertBoolean(shared_ptr<Landscape> map1, const double &deme_sampling, const double &generation)
{
    // Clear the old boolean object and set the new size
    sample_mask.setSize(y_dim, x_dim);
    for(unsigned long y = 0; y < y_dim; y++)
    {
        for(unsigned long x = 0; x < x_dim; x++)
        {
            long tmp_x = x;
            long tmp_y = y;
            long tmp_xwrap = 0;
            long tmp_ywrap = 0;
            recalculateCoordinates(tmp_x, tmp_y, tmp_xwrap, tmp_ywrap);
            double density = map1->getVal(tmp_x, tmp_y, tmp_xwrap, tmp_ywrap, generation) * deme_sampling;
            sample_mask.setValue(y, x, density >= 1.0);
        }
    }
}

void DataMask::clearSpatialMask()
{
    sample_mask_exact.setSize(0, 0);
}

void DataMask::recalculateCoordinates(long &x, long &y, long &x_wrap, long &y_wrap)
{
    if(isGridOffset)
    {
        x_wrap = (long) ((floor((x - (double) x_offset) / (double) x_dim)));
        y_wrap = (long) ((floor((y - (double) y_offset) / (double) y_dim)));
        x += -x_offset - (x_wrap * x_dim);
        y += -y_offset - (y_wrap * y_dim);
    }
}
