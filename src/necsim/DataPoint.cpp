//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Samuel Thompson
 * @date 30/08/2016
 * @file Datapoint.cpp
 *
 * @brief Contains the Datapoint class for storing objects during simulation run time.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * This class is only used during simulation runs and is not outputted to a database.
 * A Row of Datapoint objects is utilised by the main Tree objects.
 */
#include "Logging.h"
#include "DataPoint.h"

void DataPoint::setup(unsigned long x, unsigned long y, long xwrap_in, long ywrap_in, unsigned long reference_in,
                      unsigned long list_position_in, double min_max_in)
{
    xpos = x;
    ypos = y;
    xwrap = xwrap_in;
    ywrap = ywrap_in;
    next_lineage = 0;
    reference = reference_in;
    list_position = list_position_in;
    nwrap = 0;
    min_max = min_max_in;
}

void DataPoint::setup(unsigned long reference_in, unsigned long list_position_in, double min_max_in)
{
    setup(0, 0, 0, 0, reference_in, list_position_in, min_max_in);
}

void DataPoint::setup(DataPoint datin)
{
    xpos = datin.getXpos();
    ypos = datin.getYpos();
    xwrap = datin.getXwrap();
    ywrap = datin.getYwrap();
    next_lineage = datin.getNext();
//		last = datin.get_last(); // removed as of version 3.1
    reference = datin.getReference();
    list_position = datin.getListpos();
    nwrap = datin.getNwrap();
    min_max = datin.getMinmax();
}

void DataPoint::setReference(unsigned long z)
{
    reference = z;
}

void DataPoint::setNext(unsigned long x)
{
    next_lineage = x;
}

void DataPoint::setListPosition(unsigned long l)
{
    list_position = l;
}

void DataPoint::setNwrap(unsigned long n)
{
    nwrap = n;
}

void DataPoint::setMinmax(double d)
{
    min_max = d;
}

unsigned long DataPoint::getXpos()
{
    return xpos;
}

unsigned long DataPoint::getYpos()
{
    return ypos;
}

long DataPoint::getXwrap()
{
    return xwrap;
}

long DataPoint::getYwrap()
{
    return ywrap;
}

unsigned long DataPoint::getReference()
{
    return reference;
}

unsigned long DataPoint::getNext()
{
    return next_lineage;
}

unsigned long DataPoint::getListpos()
{
    return list_position;
}

unsigned long DataPoint::getNwrap()
{
    return nwrap;
}

double DataPoint::getMinmax()
{
    return min_max;
}

void DataPoint::decreaseNwrap()
{
    if(nwrap == 0)
    {
        throw out_of_range("ERROR_DATA_001: Trying to decrease  nwrap less than 0.");
    }
    else
    {
        nwrap--;
    }
}

void DataPoint::setEndpoint(long x, long y, long xwrapin, long ywrapin)
{
    xpos = x;
    ypos = y;
    xwrap = xwrapin;
    ywrap = ywrapin;
}

ostream &operator<<(ostream &os, const DataPoint &d)
{
    os << d.xpos << "," << d.ypos << "," << d.xwrap << "," << d.ywrap << "," << d.next_lineage << "," << d.reference
       << "," << d.list_position << "," << d.nwrap << ",";
    os << d.min_max << "\n";
    return os;
}

istream &operator>>(istream &is, DataPoint &d)
{
    //os << m.num_rows<<" , "<<m.num_cols<<" , "<<endl;
    char delim;
    //os << "datapoint" << endl;
    is >> d.xpos >> delim >> d.ypos >> delim >> d.xwrap >> delim >> d.ywrap >> delim >> d.next_lineage >> delim
       >> d.reference >> delim >> d.list_position >> delim >> d.nwrap >> delim;
    is >> d.min_max;
    return is;
}

#ifdef DEBUG
void DataPoint::logActive(const int &level)
{
    writeLog(50, "x, y, (x wrap, y wrap): " + to_string(xpos) + ", " + to_string(ypos) + ", (" +
                 to_string(xwrap) + ", " + to_string(ywrap) + ")");
    writeLog(50, "Lineage next: " + to_string(next_lineage));
    writeLog(50, "Reference: " + to_string(reference));
    writeLog(50, "List position: " + to_string(list_position));
    writeLog(50, "Number in wrapped lineages: " + to_string(nwrap));
    writeLog(50, "Minimum maximum: " + to_string(min_max));
}
#endif // DEBUG
