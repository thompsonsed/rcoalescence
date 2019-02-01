//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file Matrix.h
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Contains a template for a matrix with all the basic matrix operations overloaded.
 *
 * @details Code supplied by James Rosindell with large usage of
 * <a> href = "http://www.devarticles.com/c/a/Cplusplus/Operator-Overloading-in-C-plus/1"> this website </a>, and
 * modified and updated by Samuel Thompson.
 * There are two distinct classes, Row and Matrix.
 * Most operations are low-level, but some higher level functions remain, such as importCsv().
 *
 * Contact: thompsonsed@gmail.com
 */


#ifndef MATRIX
#define MATRIX
#define null 0

#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#ifdef use_csv
#include<cmath>
#include <stdexcept>
#include "fast-cpp-csv-parser/csv.h"
#endif

#include <cstdint>
#include "Logging.h"

using namespace std;

// Array of data sizes for importing tif files.
const int gdal_data_sizes[] = {0, 8, 16, 16, 32, 32, 32, 64};

/**
 * @brief Contains a template Row class and basic operations.
 * Uses an array to store the row.
 * @tparam T the type of the values in the row
 */
template<class T>
class Row
{
private:
    // stores the number of columns in the row
    unsigned long numCols{};
    // an array to store the row
    T *row;
public:

    /**
     * @brief Standard constructor.
     * @param cols optionally provide the number of rows to initiate with.
     */
    explicit Row(unsigned long cols = 0) : row(nullptr)
    {
        setSize(cols);
    }

    /**
     * @brief Standard destructor
     */
    ~Row()
    {
        if(row)
        {
            delete[] row;
        }
    }

    /**
     * @brief Copy constructor.
     * @param r the Row object to copy from.
     */
    Row(const Row &r) : row(0)
    {
        setSize(r.numCols);
        // Change to standard copy
        copy(&r.row[0], &r.row[numCols], row);
    }

    /**
     * @brief Setter for the row size.
     * @param n the number of rows to initiate with.
     * SetRowSize() deletes any old data, and allocates space for new data, unless we set the number of columns to 0, in which case it merely deletes the data.
     * This lets us use this function for construction, destruction, and dynamic modification in one method.
     */
    void setSize(unsigned long n)
    {
        if(row)
        {
            delete[] row;
        }
        if(n > 0)
        {
            row = new T[n];
        }
        else
        {
            row = 0;
        }
        numCols = n;
    }

    /**
     * @brief Changes the size of the array.
     * @param n the new size to change to.
     * Note that no checks are performed that the new row size is larger than the old row size.
     * Thus is this function is used to shrink the row size, a bad_alloc error will likely be thrown.
     */
    void resize(unsigned long n)
    {
        auto t = new T[n];
        for(unsigned long i = 0; i < numCols; i++)
        {
            t[i] = row[i];
        }
        delete[] row;
        row = move(t);
        numCols = n;
    }

    /**
     * @brief Getter for the size of the array.
     * @return the number of columns.
     */
    unsigned long size()
    {
        return numCols;
    }

    /**
     * @brief Overloading the [] operator to allow for simple referencing.
     * @param column the column to get the value from.
     * @return the value in the specified column.
     * Note that different versions deal with values outside of (0,numCols) in different ways.
     * @note updated to throw an out_of_range exception if the column is out of the row range.
     */
    T &operator[](unsigned long column)
    {
        // assert(column<num_cols);
        // check we are within bounds
#ifdef DEBUG
        if(column < 0 || column >= numCols)
        {
            string err =
                    "ERROR_MAIN_013b: Tried to call an indices that was out of range of the row. Check row size definition. num_cols: " +
                    to_string((long long) numCols) + " index: " + to_string((long long) column);
            throw out_of_range(err);
        }
#endif
        column = column % numCols;
        return row[column];
    }

    /**
     * @brief Overloading the = operator to allow for copying data across.
     * @param r the Row object to copy data from.
     */
    Row &operator=(const Row &r)
    {
        setSize(r.numCols);
        for(unsigned long i = 0; i < numCols; i++)
        {
            row[i] = r.row[i];
        }
        return *this;
    }

    /**
     * @brief Overloading the << operator for outputting to the output stream.
     * @param os the output stream.
     * @param r the Row object to output from.
     * @return os the output stream.
     */
    friend ostream &operator<<(ostream &os, const Row &r)
    {
        os << r.numCols << ",";
        for(unsigned long c = 0; c < r.numCols; c++)
        {
            os << r.row[c] << ",";
        }
        return os;
    }

    /**
     * @brief Overloading the << operator for inputting from an input stream.
     * @param is the input stream.
     * @param r the Row object to input to.
     * @return the input stream.
     */
    friend istream &operator>>(istream &is, Row &r)
    {
        char delim;
        int n;
        is >> n;
        r.setSize(n);
        is >> delim;
        for(unsigned long c = 0; c < r.numCols; c++)
        {
            is >> r.row[c];
            is >> delim;
        }
        return is;
    }
};

/**
 * @brief A class containing the Matrix object, set up as an array of Row objects.
 * Includes basic operations, as well as the importCsv() function for more advanced reading from file.
 * @tparam T the type of the values in the matrix
 */
template<class T>
class Matrix
{

protected:

    // number of rows and columns
    unsigned long num_cols{};
    unsigned long num_rows{};
    // a matrix is an array of rows
    Row<T> *matrix;
public:

    /**
     * @brief The standard constructor
     * @param rows optionally provide the number of rows.
     * @param cols optionally provide the number of columns.
     */
    explicit Matrix(unsigned long rows = 0, unsigned long cols = 0) : matrix(nullptr)
    {
        setSize(rows, cols);
    }

    /**
     * @brief The copy constructor.
     * @param m a Matrix object to copy from.
     */
    Matrix(const Matrix &m) : matrix(nullptr)
    {
        setSize(m.num_rows, m.num_cols);
        copy(&m.matrix[0][0], &m.matrix[num_rows][num_cols], matrix);
    }

    /**
    * @brief The destructor.
    */
    virtual ~Matrix()
    {
        if(matrix)
        {
            delete[] matrix;
        }
    }

    /**
     * @brief Sets the matrix size.
     * Similar concept to that for Rows.
     * @param rows the number of rows.
     * @param cols the number of columns.
     */
    void setSize(unsigned long rows, unsigned long cols)
    {
        if(matrix)
        {
            delete[]matrix;
        }
        if(cols > 0 && rows > 0)
        {
            matrix = new Row<T>[rows];
            for(unsigned long i = 0; i < rows; i++)
            {
                matrix[i].setSize(cols);
            }
        }
        else
        {
            matrix = nullptr;
        }
        num_cols = cols;
        num_rows = rows;
    }

    /**
     * @brief Getter for the number of columns.
     * @return the number of columns.
     */
    unsigned long getCols() const
    {
        return num_cols;
    }

    /**
     * @brief Getter for the number of rows.
     * @return the number of rows.
     */
    unsigned long getRows() const
    {
        return num_rows;
    }

    /**
     * @brief Overoads the [] operator for Matrix.
     * Allows referencing of a value i,j using Matrix[i][j].
     * Includes error checking for if the indices are out of range of the matrix.
     * Note that this functionality has been altered since the original file generation.
     * @param index the row number to get the value from.
     * @return the matrix row object.
     */
    Row<T> &operator[](unsigned long index)
    {
#ifdef DEBUG
        if(index < 0 || index >= num_rows)
        {
            string err =
                    "ERROR_MAIN_013: Tried to call an indices that was out of range of the matrix. Check matrix size definition. num_rows: " +
                    to_string((long long) num_rows) + " index: " + to_string((long long) index);
            throw out_of_range(err);
        }
#endif
        index = index % num_rows;
        return matrix[index];
    }

    /**
     * @brief Overloading the = operator.
     * @param m the matrix to copy from.
     */
    Matrix &operator=(const Matrix &m)
    {
        setSize(m.num_rows, m.num_cols);
        for(unsigned long r = 0; r < num_rows; r++)
        {
            matrix[r] = Row<T>(m.matrix[r]);
        }
        return *this;
    }

    /**
     * @brief Overloading the + operator.
     * @param m the matrix to add to this matrix.
     * @return the matrix object which is the sum of the two matrices.
     */
    Matrix operator+(const Matrix &m) const
    {
        //Since addition creates a new matrix, we don't want to return a reference, but an actual matrix object.
        unsigned long newnumcols, newnumrows;
        if(num_cols > m.num_cols)
        {
            newnumcols = m.num_cols;
        }
        else
        {
            newnumcols = num_cols;
        }
        if(num_rows > m.num_rows)
        {
            newnumrows = m.num_rows;
        }
        else
        {
            newnumrows = num_rows;
        }

        Matrix result(newnumrows, newnumcols);
        for(unsigned long r = 0; r < newnumrows; r++)
        {
            for(unsigned long c = 0; c < newnumcols; c++)
            {
                result[r][c] = matrix[r][c] + m.matrix[r][c];
            }
        }
        return result;
    }

    /**
    * @brief Overloading the - operator.
     * @param m the matrix to subtract from this matrix.
     * @return the matrix object which is the subtraction of the two matrices.
     * */
    Matrix operator-(const Matrix &m) const
    {
        unsigned long newnumcols, newnumrows;
        if(num_cols > m.num_cols)
        {
            newnumcols = m.num_cols;
        }
        else
        {
            newnumcols = num_cols;
        }
        if(num_rows > m.num_rows)
        {
            newnumrows = m.num_rows;
        }
        else
        {
            newnumrows = num_rows;
        }
        Matrix result(newnumrows, newnumcols);
        for(unsigned long r = 0; r < newnumrows; r++)
        {
            for(unsigned long c = 0; c < newnumcols; c++)
            {
                result[r][c] = matrix[r][c] - m.matrix[r][c];
            }
        }
        return result;
    }

    /**
     * @brief Overloading the += operator so that the new object is written to the current object.
     * @param m the Matrix object to add to this matrix.
     */
    Matrix &operator+=(const Matrix &m)
    {
        unsigned long newnumcols, newnumrows;
        if(num_cols > m.num_cols)
        {
            newnumcols = m.num_cols;
        }
        else
        {
            newnumcols = num_cols;
        }
        if(num_rows > m.num_rows)
        {
            newnumrows = m.num_rows;
        }
        else
        {
            newnumrows = num_rows;
        }
        for(unsigned long r = 0; r < newnumrows; r++)
        {
            for(unsigned long c = 0; c < newnumcols; c++)
            {
                matrix[r][c] += m.matrix[r][c];
            }
        }
        return *this;
    }

    /**
     * @brief Overloading the -= operator so that the new object is written to the current object.
     * @param m the Matrix object to subtract from this matrix.
     */
    Matrix &operator-=(const Matrix &m)
    {
        unsigned long newnumcols, newnumrows;
        if(num_cols > m.num_cols)
        {
            newnumcols = m.num_cols;
        }
        else
        {
            newnumcols = num_cols;
        }
        if(num_rows > m.num_rows)
        {
            newnumrows = m.num_rows;
        }
        else
        {
            newnumrows = num_rows;
        }
        for(unsigned long r = 0; r < newnumrows; r++)
        {
            for(unsigned long c = 0; c < newnumcols; c++)
            {
                matrix[r][c] -= m.matrix[r][c];
            }
        }
        return *this;
    }

    /**
     * @brief Overloading the * operator for scaling.
     * @param s the constant to scale the matrix by.
     * @return the scaled matrix.
     */
    Matrix operator*(const double s) const
    {
        Matrix result(num_rows, num_cols);
        for(unsigned long r = 0; r < num_rows; r++)
        {
            for(unsigned long c = 0; c < num_cols; c++)
            {
                result[r][c] = matrix[r][c] * s;
            }
        }
        return result;
    }

    /**
     * @brief Overloading the * operator for matrix multiplication.
     * Multiplies each value in the matrix with its corresponding value in the other matrix.
     * @param m the matrix to multiply with
     * @return the product of each ith,jth value of the matrix.
     */
    Matrix operator*(Matrix &m) const
    {
        unsigned long newnumcols;
        if(num_cols > m.num_rows)
        {
            newnumcols = m.num_rows;
        }
        else
        {
            newnumcols = num_cols;
        }

        Matrix result(num_rows, m.num_cols);
        for(unsigned long r = 0; r < num_rows; r++)
        {
            for(unsigned long c = 0; c < m.num_cols; c++)
            {
                for(unsigned long i = 0; i < newnumcols; i++)
                {
                    result[r][c] += matrix[r][i] * m[i][c];
                }
            }
        }
        return result;
    }

    /**
     * @brief Writes the object to the output stream.
     * @param os the output stream to write to
     * @param m the object to write out
     * @return the output stream
     */
    friend ostream &writeOut(ostream &os, const Matrix &m)
    {
        for(unsigned long r = 0; r < m.num_rows; r++)
        {
            for(unsigned long c = 0; c < m.num_cols; c++)
            {
                os << m.matrix[r][c] << ",";
            }
            os << "\n";
        }
        return os;
    }

    /**
     * @brief Reads in from the input stream.
     * @param is the input stream to read from
     * @param m the object to read into
     * @return
     */
    friend istream &readIn(istream &is, Matrix &m)
    {
        char delim;
        for(unsigned long r = 0; r < m.num_rows; r++)
        {
            for(unsigned long c = 0; c < m.num_cols; c++)
            {
                is >> m.matrix[r][c];
                is >> delim;
            }
        }
        return is;
    }

    /**
     * @brief Overloading the << operator for outputting to an output stream.
     * This can be used for writing to console or storing to file.
     * @param os the output stream.
     * @param m the matrix to output.
     * @return the output stream.
     */
    friend ostream &operator<<(ostream &os, const Matrix &m)
    {
        return writeOut(os, m);
    }

    /**
     * @brief Overloading the >> operator for inputting from an input stream.
     * This can be used for writing to console or storing to file.
     * @param is the input stream.
     * @param m the matrix to input to.
     * @return the input stream.
     */
    friend istream &operator>>(istream &is, Matrix &m)
    {
        return readIn(is, m);
    }

    /**
     * @brief Sets the value at the specified indices, including handling type conversion from char to the template
     * class.
     * @param x the x index.
     * @param y the y index.
     * @param value the value to set
     */
    void setValue(const unsigned long &x, const unsigned long &y, const char *value)
    {
        matrix[y][x] = static_cast<T>(*value);
    }

    /**
     * @brief Imports the matrix from a csv file.
     *
     * @throws runtime_error: if type detection for the filename fails.
     * @param filename the file to import.
     */
    virtual void import(const string &filename)
    {
        if(!importCsv(filename))
        {
            string s = "Type detection failed for " + filename + ". Check file_name is correct.";
            throw runtime_error(s);
        }
    }

    /**
     * @brief Imports the matrix from a csv file using the fast-csv-parser method.
     * @param filename the path to the file to import.
     */
#ifdef use_csv
    bool importCsv(const string &file_name)
    {
    if(file_name.find(".csv") != string::npos)
        {
            stringstream os;
            os  << "Importing " << file_name << " " << flush;
            writeInfo(os.str());
            // LineReader option
            io::LineReader in(file_name);
            // Keep track of whether we've printed to terminal or not.
            bool bPrint = false;
            // Initialies empty variable so that the setValue operator overloading works properly.
            unsigned int number_printed = 0;
            for(unsigned long i =0; i<num_rows; i++)
            {
                char* line = in.next_line();
                if(line == nullptr)
                {
                    if(!bPrint)
                    {
                        writeError("Input dimensions incorrect - read past end of file.");
                        bPrint = true;
                    }
                    break;
                }
                else
                {
                    char *dToken;
                    dToken = strtok(line,",");
                    for(unsigned long j = 0; j<num_cols; j++)
                    {
                        if(dToken == nullptr)
                        {
                            if(!bPrint)
                            {
                            writeError("Input dimensions incorrect - read past end of file.");
                                bPrint = true;
                            }
                            break;
                        }
                        else
                        {
                            // This function is overloaded to correctly determine the type of the template
                            setValue(j,i,dToken);
                            dToken = strtok(NULL,",");
                        }
                    }
                    // output the percentage complete
                    double dComplete = ((double)i/(double)num_rows)*20;
                    if( number_printed < dComplete)
                    {
                        stringstream os;
                        os  << "\rImporting " << file_name << " ";
                        number_printed = 0;
                        while(number_printed < dComplete)
                        {
                            os << ".";
                            number_printed ++;
                        }
                        os << flush;
                        writeInfo(os.str());
                    }

                }
            }
            writeInfo("done.\n");
            return true;
        }
        return false;
    }
#endif
#ifndef use_csv

    /**
     * @brief Imports the matrix from a csv file using the standard, slower method.
     * @deprecated this function should not be used any more as it is much slower.
     * @param filename the path to the file to import.
     * @return true if the csv can be imported.
     */
    bool importCsv(const string &filename)
    {
        if(filename.find(".csv") != string::npos)
        {
            stringstream os;
            os << "Importing" << filename << " " << flush;
            ifstream inputstream;
            inputstream.open(filename.c_str());
            unsigned long number_printed = 0;
            for(uint32_t j = 0; j < num_rows; j++)
            {
                string line;
                getline(inputstream, line);
                istringstream iss(line);
                for(uint32_t i = 0; i < num_cols; i++)
                {
                    char delim;
                    T val;
                    iss >> val >> delim;
                    matrix[j][i] = val;
                }
                double dComplete = ((double) j / (double) num_rows) * 5;
                if(number_printed < dComplete)
                {
                    os << "\rImporting " << filename << " " << flush;
                    while(number_printed < dComplete)
                    {
                        os << ".";
                        number_printed++;
                    }
                    os << flush;
                    writeInfo(os.str());

                }
            }
            stringstream os2;
            os2 << "\rImporting" << filename << "..." << "done." << "                          " << endl;
            inputstream.close();
            writeInfo(os2.str());
            return true;
        }
        return false;
    }

#endif // use_csv
};

#endif // MATRIX
