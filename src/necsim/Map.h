// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details
/**
 * @author Samuel Thompson
 * @file Map.h
 * @brief Contains Map for importing .tif files and obtaining a variety of information from them.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef MAP_H
#define MAP_H
#ifdef with_gdal

#include <string>
#include <cstring>
#include <sstream>
#include <gdal_priv.h>
#include <cpl_conv.h> // for CPLMalloc()
#include <sstream>
#include "Logging.h"
#include "Matrix.h"
#include "custom_exceptions.h"
#include "cpl_custom_handler.h"

using namespace std;
#ifdef DEBUG

#include "custom_exceptions.h"

#endif // DEBUG

/**
 * @brief Read a a tif file to a matrix and obtain spatial metadata.
 * @tparam T The type of the Matrix to create.
 */
template<class T>
class Map : public virtual Matrix<T>
{
protected:
    GDALDataset* po_dataset;
    GDALRasterBand* po_band;
    unsigned long block_x_size, block_y_size;
    double no_data_value;
    string file_name;
    GDALDataType gdal_data_type;
    // Contains the error object to check for any gdal issues
    CPLErr cpl_error;
    double upper_left_x, upper_left_y, x_res, y_res;
    using Matrix<T>::matrix;
    using Matrix<T>::num_cols;
    using Matrix<T>::num_rows;
    bool cpl_error_set;
public:
    using Matrix<T>::setSize;
    using Matrix<T>::getCols;
    using Matrix<T>::getRows;
    using Matrix<T>::get;
    using Matrix<T>::operator*;
    using Matrix<T>::operator/;
    using Matrix<T>::operator+;
    using Matrix<T>::operator-;
    using Matrix<T>::operator-=;
    using Matrix<T>::operator+=;
    using Matrix<T>::operator*=;
    using Matrix<T>::operator/=;

    Map() : Matrix<T>(0, 0), po_dataset(nullptr), po_band(nullptr), block_x_size(0), block_y_size(0),
            no_data_value(0.0), file_name(""), gdal_data_type(GDT_Unknown), cpl_error(CE_None), upper_left_x(0.0),
            upper_left_y(0.0), x_res(0.0), y_res(0.0), cpl_error_set(false)
    {
        setCPLErrorHandler();
        GDALAllRegister();
        removeCPLErrorHandler();
    }

    ~Map() override
    {
        close();
        removeCPLErrorHandler();
    }

    void setCPLErrorHandler()
    {
        if(!cpl_error_set)
        {
            cpl_error_set = true;
            CPLPushErrorHandler(cplNecsimCustomErrorHandler);
        }
    }

    void removeCPLErrorHandler()
    {
        if(cpl_error_set)
        {
            cpl_error_set = false;
            CPLPopErrorHandler();
        }
    }

    /**
     * @brief Opens the provided filename to the poDataset object.
     * @param filename file to open in read-only mode.
     */
    void open(const string &filename_in)
    {
        if(!po_dataset)
        {
            file_name = filename_in;
            setCPLErrorHandler();
            po_dataset = (GDALDataset*) GDALOpen(file_name.c_str(), GA_ReadOnly);
            removeCPLErrorHandler();
        }
        else
        {
            throw FatalException("File already open at " + file_name);
        }
        if(!po_dataset)
        {
            string s = "File " + file_name + " not found.";
            throw runtime_error(s);
        }
    }

    /**
     * @brief Overloaded open for using the preset file name.
     */
    void open()
    {
        open(file_name);
    }

    /**
     * @brief Checks if the connection to the map file has already been opened.
     *
     * All this does is check if poDataset is a null pointer.
     * @return true if poDataset is a null pointer.
     */
    bool isOpen()
    {
        return po_dataset != nullptr;
    }

    /**
     * @brief Destroys the connection to the dataset.
     */
    void close()
    {
        if(po_dataset)
        {
            setCPLErrorHandler();
            GDALClose(po_dataset);
            removeCPLErrorHandler();
            po_dataset = nullptr;
            po_band = nullptr;
        }
    }

    /**
     * @brief Sets the raster band to the first raster.
     */
    void getRasterBand()
    {
        setCPLErrorHandler();
        if(po_dataset->GetRasterCount() < 1)
        {
            stringstream ss;
            ss << "No raster band detected in " << file_name << endl;
            throw FatalException(ss.str());
        }
        po_band = po_dataset->GetRasterBand(1);
        removeCPLErrorHandler();
        if(po_band == nullptr)
        {
            stringstream ss;
            ss << "Could not open po_band in " << file_name << ": returned a nullptr" << endl;
            throw FatalException(ss.str());
        }
    }

    /**
     * @brief Obtains the x and y dimensions from the tif file for reading in blocks.
     */
    void getBlockSizes()
    {
        setCPLErrorHandler();
        block_x_size = static_cast<unsigned long>(po_dataset->GetRasterXSize());
        block_y_size = static_cast<unsigned long>(po_dataset->GetRasterYSize());
        removeCPLErrorHandler();
    }

    /**
     * @brief Sets the no data, data type and data type name values from the tif file.
     */
    void getMetaData()
    {
#ifdef DEBUG
        if(po_dataset == nullptr)
        {
            throw FatalException("po_dataset is nullptr. This is likely a problem with gdal. Please report this bug.");
        }
        if(po_band == nullptr)
        {
            throw FatalException("po_band is nullptr. This is likely a problem with gdal. Please report this bug.");
        }
#endif // DEBUG
        setCPLErrorHandler();
        try
        {
            int pbSuccess;
            no_data_value = po_band->GetNoDataValue(&pbSuccess);
            if(!pbSuccess)
            {
                no_data_value = 0.0;
            }
        }
        catch(out_of_range &out_of_range1)
        {
            no_data_value = 0.0;
        }
        stringstream ss;
        ss << "No data value is: " << no_data_value << endl;
        writeInfo(ss.str());
        // Check sizes match
        gdal_data_type = po_band->GetRasterDataType();
        const static double default_values[6] = {0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        double geoTransform[6] = {0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        writeInfo("Getting geo transform...");
        cpl_error = po_dataset->GetGeoTransform(geoTransform);
        if(cpl_error >= CE_Warning)
        {
            ss.str("");
            ss << "cpl error detected greater than or equal to warning level." << endl;
            writeInfo(ss.str());
            string message = "No transform present in dataset for " + file_name;
            cplNecsimCustomErrorHandler(cpl_error, 6, message.c_str());
            CPLErrorReset();
            cpl_error = CE_None;
            std::copy(std::begin(default_values), std::end(default_values), std::begin(geoTransform));
        }
        writeInfo("done.\n");
        ss.str("");
        ss << "Affine transform is " << geoTransform[0] << ", " << geoTransform[1] << ", " << geoTransform[2] << ", ";
        ss << geoTransform[3] << ", " << geoTransform[4] << ", " << geoTransform[5] << endl;
        writeInfo(ss.str());
        upper_left_x = geoTransform[0];
        upper_left_y = geoTransform[3];
        x_res = geoTransform[1];
        y_res = -geoTransform[5];
        removeCPLErrorHandler();
        //		checkTifImportFailure();
#ifdef DEBUG
        printMetaData();
#endif // DEBUG
    }

#ifdef DEBUG

    void printMetaData()
    {
        setCPLErrorHandler();
        stringstream ss;
        const char* dt_name = GDALGetDataTypeName(gdal_data_type);
        ss << "Filename: " << file_name << endl;
        writeLog(10, ss.str());
        ss.str("");
        ss << "data type: " << gdal_data_type << "(" << dt_name << ")" << endl;
        writeLog(10, ss.str());
        ss.str("");
        ss << "Geo-transform (ulx, uly, x res, y res): " << upper_left_x << ", " << upper_left_y << ", ";
        ss << x_res << ", " << y_res << ", " << endl;
        writeLog(10, ss.str());
        ss.str("");
        ss << "No data value: " << no_data_value << endl;
        writeLog(10, ss.str());
        removeCPLErrorHandler();

    }

#endif //DEBUG

    /**
     * @brief Gets the upper left x (longitude) coordinate
     * @return upper left x of the map
     */
    double getUpperLeftX()
    {
        return upper_left_x;
    }

    /**
     * @brief Gets the upper left y (latitude) coordinate
     * @return upper left y of the map
     */
    double getUpperLeftY()
    {
        return upper_left_y;
    }

    /**
     * @brief Imports the matrix from a csv file.
     *
     * @throws runtime_error: if type detection for the filename fails.
     * @param filename the file to import.
     */
    void import(const string &filename) override
    {
        if(!importTif(filename))
        {
            Matrix<T>::import(filename);
        }
    }

    /**
     * @brief Imports the matrix from a tif file using the gdal library functions.
     * @note Opens a connection to the file object, which should be closed.
     * @param filename the path to the file to import.
     */
    bool importTif(const string &filename)
    {

        if(filename.find(".tif") != string::npos)
        {
            setCPLErrorHandler();
            stringstream ss;
            ss << "Importing " << filename << " " << endl;
            writeInfo(ss.str());
            open(filename);
            getRasterBand();
            getBlockSizes();
            getMetaData();
            // If the sizes are 0 then use the raster sizes
            if(num_cols == 0 || num_rows == 0)
            {
                setSize(block_y_size, block_x_size);
            }
            // Check sizes
            if((num_cols != block_x_size || num_rows != block_y_size) || num_cols == 0 ||
               num_rows == 0)
            {
                stringstream stringstream1;
                stringstream1 << "Raster data size does not match inputted dimensions for " << filename
                              << ". Using raster sizes."
                              << endl;
                stringstream1 << "Old dimensions: " << num_cols << ", " << num_rows << endl;
                stringstream1 << "New dimensions: " << block_x_size << ", " << block_y_size << endl;
                writeWarning(stringstream1.str());
                setSize(block_y_size, block_x_size);
            }
            // Check the data types are support
            const char* dt_name = GDALGetDataTypeName(gdal_data_type);
            if(gdal_data_type == 0 || gdal_data_type > 7)
            {
                removeCPLErrorHandler();
                throw FatalException("Data type of " + string(dt_name) + " is not supported.");
            }
#ifdef DEBUG
            if(sizeof(T) * 8 != gdal_data_sizes[gdal_data_type])
            {
                stringstream ss2;
                ss2 << "Object data size: " << sizeof(T) * 8 << endl;
                ss2 << "Tif data type: " << dt_name << ": " << gdal_data_sizes[gdal_data_type] << " bytes" << endl;
                ss2 << "Tif data type does not match object data size in " << file_name << endl;
                writeWarning(ss2.str());
            }
#endif
            // Use the overloaded method for importing between types
            internalImport();
            writeInfo("done.\n");
            removeCPLErrorHandler();
            return true;
        }
        return false;
    }

    /**
     * @brief Opens the offset map and fetches the metadata.
     * @param offset_map the offset map to open (should be the larger map).
     * @return true if the offset map is opened within this function
     */
    bool openOffsetMap(Map &offset_map)
    {
        bool opened_here = false;
        if(!offset_map.isOpen())
        {
            opened_here = true;
            offset_map.open();
        }
        offset_map.getRasterBand();
        offset_map.getMetaData();
        return opened_here;
    }

    void closeOffsetMap(Map &offset_map, const bool &opened_here)
    {
        if(opened_here)
        {
            offset_map.close();
        }
    }

    /**
     * @brief Calculates the offset between the two maps.
     *
     * The offset_map should be larger and contain this map, otherwise returned values will be negative
     *
     * @note Opens a connection to the tif file (if it has not already been opened), which is then closed. If the
     * 		 connection is already open, then it will not be closed and it is assumed logic elsewhere achieves this.
     *
     * @note Offsets are returned as rounded integers at the resolution of the smaller map.
     *
     * @param offset_map the offset map to read from
     * @param offset_x the x offset variable to fill
     * @param offset_y the y offset variable to fill
     */
    void calculateOffset(Map &offset_map, long &offset_x, long &offset_y)
    {
        auto opened_here = openOffsetMap(offset_map);
        offset_x = static_cast<long>(round((upper_left_x - offset_map.upper_left_x) / x_res));
        offset_y = static_cast<long>(round((offset_map.upper_left_y - upper_left_y) / y_res));
        closeOffsetMap(offset_map, opened_here);
    }

    /**
     * @brief Calculates the relative scale of this map compared to the offset map.
     *
     * The offset map should be larger and contain this map.
     *
     * @note Only the x resolution is checked, it is assumed the x and y resolutions of both maps is the same (i.e. each
     * 		 cell on the map is a square.
     *
     * @param offset_map the offset map object to read from
     * @return the relative scale of the offset map
     */
    unsigned long roundedScale(Map &offset_map)
    {
        auto opened_here = openOffsetMap(offset_map);
        closeOffsetMap(offset_map, opened_here);
        return static_cast<unsigned long>(floor(offset_map.x_res / x_res));
    }

    /**
     * @brief Default importer when we rely on the default gdal method of converting between values.
     * Note that importing doubles to ints results in the values being rounded down.
     * @return true if a tif file exists and can be imported, false otherwise.
     */
    void internalImport()
    {
        writeWarning(
                "\nNo type detected for Map type. Attempting default importing (potentially undefined behaviour).");
        defaultImport();
    }

    /**
     * @brief Default import routine for any type. Provided as a separate function so implementation can be called from
     * any template class type.
     */
    void defaultImport()
    {
        setCPLErrorHandler();
        unsigned int number_printed = 0;
        for(uint32_t j = 0; j < num_rows; j++)
        {
            printNumberComplete(j, number_printed);
            cpl_error = po_band->RasterIO(GF_Read, 0, j, static_cast<int>(block_x_size), 1, &matrix[j*num_cols],
                                          static_cast<int>(block_x_size), 1, gdal_data_type, 0, 0);
            checkTifImportFailure();
            // Now convert the no data values to 0
            for(uint32_t i = 0; i < num_cols; i++)
            {
                if(get(j, i) == no_data_value)
                {
                    get(j, i) = 0;
                }
            }
        }
        removeCPLErrorHandler();
    }

    /**
     * @brief Imports from the supplied filename into the GeoTiff object, converting doubles to booleans.
     * The threshold for conversion is x>0.5 -> true, false otherwise.
     */
    void importFromDoubleAndMakeBool()
    {
        setCPLErrorHandler();
        writeInfo("\nConverting from double to boolean.\n");
        unsigned int number_printed = 0;
        // create an empty row of type float
        double* t1;
        t1 = (double*) CPLMalloc(sizeof(double) * num_cols);
        // import the data a row at a time, using our template row.
        for(uint32_t j = 0; j < num_rows; j++)
        {
            printNumberComplete(j, number_printed);
            cpl_error = po_band->RasterIO(GF_Read, 0, j, static_cast<int>(block_x_size), 1, &t1[0],
                                          static_cast<int>(block_x_size), 1, GDT_Float64, 0, 0);
            checkTifImportFailure();
            // now copy the data to our Map, converting float to int. Round or floor...? hmm, floor?
            for(unsigned long i = 0; i < num_cols; i++)
            {
                if(t1[i] == no_data_value)
                {
                    this->setValue(j, i, false);
                }
                else
                {
                    this->setValue(j, i, t1[i] >= 0.5);
                }
            }
        }
        CPLFree(t1);
        removeCPLErrorHandler();
    }

    /**
     * @brief Imports from the supplied filename into the GeoTiff object, converting doubles to booleans.
     * The threshold for conversion is x>0.5 -> true, false otherwise.
     *
     * @param dt_buff: the buffer type for the data
     * @tparam T2 the template type for data reading.
     */
    template<typename T2>
    void importUsingBuffer(GDALDataType dt_buff)
    {
        setCPLErrorHandler();
        stringstream ss;
        ss << "\nUsing buffer of type " << GDALGetDataTypeName(dt_buff) << " into array of dimensions ";
        ss << num_cols << " by " << num_rows << " with block size " << block_x_size << ", " << block_y_size << endl;
        writeInfo(ss.str());
        unsigned int number_printed = 0;
        // create an empty row of type float
        T2* t1;
        t1 = (T2*) CPLMalloc(sizeof(T2) * num_cols);
        if(CPLGetLastErrorNo() != CE_None)
        {
            stringstream os;
            os << "Error thrown by CPLMalloc: " << CPLGetLastErrorType() << ": " << CPLGetLastErrorMsg() << endl;
            throw FatalException(os.str());
        }
        // import the data a row at a time, using our template row.
        auto int_block_x_size = static_cast<int>(block_x_size);
#ifdef DEBUG
        if(sizeof(T2) * 8 != GDALGetDataTypeSize(dt_buff))
        {
            stringstream ss0;
            ss0 << "Size of template (" << sizeof(T2) << ") does not equal size of gdal buffer (";
            ss0 << GDALGetDataTypeSize(dt_buff) << "). Please report this bug." << endl;
            throw FatalException(ss0.str());
        }
        if(t1 == nullptr)
        {
            stringstream ss1;
            ss1 << "CPL malloc could not acquire " << sizeof(T2) * num_cols << " of space and returned a nullptr.";
            ss1 << " Please check that your install of gdal is fully functional, and otherwise report this bug."
                << endl;
            throw FatalException(ss1.str());
        }
        if(static_cast<unsigned long>(int_block_x_size) != block_x_size)
        {
            stringstream ss2;
            ss2 << "Error in integer conversion - please report this bug: " << int_block_x_size << " != ";
            ss2 << block_x_size << endl;
            throw FatalException(ss2.str());
        }
        if(po_band == nullptr)
        {
            CPLFree(t1);
            throw FatalException("po_band is nullptr during import using buffer - please report this bug.");
        }
#endif // DEBUG
        for(uint32_t j = 0; j < num_rows; j++)
        {
            printNumberComplete(j, number_printed);
            cpl_error = po_band->RasterIO(GF_Read, 0, j, int_block_x_size, 1, t1,
                                          int_block_x_size, 1, dt_buff, 0, 0);
            checkTifImportFailure();
            for(unsigned long i = 0; i < num_cols; i++)
            {
                if(t1[i] == no_data_value)
                {
                    this->setValue(j, i, static_cast<T>(0));
                }
                else
                {
                    this->setValue(j, i, static_cast<T>(t1[i]));
                }
            }
        }
        CPLFree(t1);
        removeCPLErrorHandler();
    }

    /**
     * @brief Print the percentage complete during import
     * @param j the reference for the counter
     * @param number_printed the number of previously printed lines
     */
    void printNumberComplete(const uint32_t &j, unsigned int &number_printed)
    {
        double dComplete = ((double) j / (double) num_rows) * 20;
        if(number_printed < dComplete)
        {
            stringstream os;
            os << "\rImporting " << file_name << " ";
            number_printed = 0;
            while(number_printed < dComplete)
            {
                os << ".";
                number_printed++;
            }
            os << flush;
            writeInfo(os.str());
        }
    }

    /**
     * @brief Checks the error code of the CPLErr object and formats the error
     */
    void checkTifImportFailure()
    {
        if(cpl_error >= CE_Warning && !cpl_error_set)
        {
            stringstream ss;
            ss << "\nCPL error thrown during import of " << file_name << endl;
            cplNecsimCustomErrorHandler(cpl_error, 3, ss.str().c_str());
            CPLErrorReset();
            cpl_error = CE_None;
        }
    }

    /**
     * @brief Output operator
     * @param os the output stream
     * @param m the object to write out
     * @return the modified output stream
     */
    friend ostream &operator>>(ostream &os, const Map &m)
    {
        return Matrix<T>::writeOut(os, m);
    }

    /**
     * @brief Input operator
     * @param is the input stream
     * @param m the object to write in
     * @return the modified input stream
     */
    friend istream &operator<<(istream &is, Map &m)
    {
        return Matrix<T>::readIn(is, m);
    }

    /**
     * @brief Equality operator
     * @param m the Map object to copy from
     * @return the self Map object
     */
    Map &operator=(const Map &m)
    {
        Matrix<T>::operator=(m);
        this->po_dataset = m.po_dataset;
        this->po_band = m.po_band;
        this->block_x_size = m.block_x_size;
        this->block_y_size = m.block_y_size;
        this->no_data_value = m.no_data_value;
        this->file_name = m.file_name;
        this->gdal_data_type = m.gdal_data_type;
        this->cpl_error = m.cpl_error;
        this->upper_left_x = m.upper_left_x;
        this->upper_left_y = m.upper_left_y;
        this->x_res = m.x_res;
        this->y_res = m.y_res;
        return *this;
    }

};

/**
 * @brief Overloaded imported for handling conversion of types to boolean. This function should only be once
 * elsewhere, so inlining is fine, allowing this file to remain header only.
 */
template<>
inline void Map<bool>::internalImport()
{
    if(gdal_data_type <= 5)
    {
        // Then the tif file type is an int/byte
        // we can just import as it is
        importUsingBuffer<uint8_t>(GDT_Byte);
    }
    else
    {
        // Conversion from double to boolean
        importFromDoubleAndMakeBool();
    }
}

/**
 * @brief Overloaded functions for importing from tifs and matching between gdal and C types.
 */
template<>
inline void Map<int8_t>::internalImport()
{
    importUsingBuffer<int16_t>(GDT_Int16);
}

template<>
inline void Map<uint8_t>::internalImport()
{
    gdal_data_type = GDT_Byte;
    defaultImport();
}

template<>
inline void Map<int16_t>::internalImport()
{
    gdal_data_type = GDT_Int16;
    defaultImport();
}

template<>
inline void Map<uint16_t>::internalImport()
{
    gdal_data_type = GDT_UInt16;
    defaultImport();
}

template<>
inline void Map<int32_t>::internalImport()
{
    gdal_data_type = GDT_Int32;
    defaultImport();
}

template<>
inline void Map<uint32_t>::internalImport()
{
    gdal_data_type = GDT_UInt32;
    defaultImport();
}

template<>
inline void Map<float>::internalImport()
{
    gdal_data_type = GDT_Float32;
    defaultImport();
}

template<>
inline void Map<double>::internalImport()
{
    gdal_data_type = GDT_Float64;
    defaultImport();
}

#endif // with_gdal
#endif //MAP_H
