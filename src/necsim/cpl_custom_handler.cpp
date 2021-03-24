// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file CPLCustomHandlerNecsim.cpp
 * @brief Contains a custom CPL error handler.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#include <sstream>
#include "cpl_custom_handler.h"
#include "Logging.h"

namespace necsim
{
    void cplNecsimCustomErrorHandler(CPLErr eErrClass, int err_no, const char* msg)
    {
        std::stringstream error_msg;
        if(!loggerIsSetup())
        {
#ifndef use_rcoalescence
            std::cerr << "Logging object has not been set before CPL error thrown: " << err_no << ". " << msg << std::endl;
#endif // use_rcoalescence
        }
        else
        {
            if(eErrClass == CE_Fatal)
            {
                error_msg << "Critical gdal error: " << err_no << ". " << msg << std::endl;
                writeCritical(error_msg.str());
            }
            else if(eErrClass == CE_Failure)
            {
                error_msg << "Gdal error: " << err_no << ". " << msg << std::endl;
                writeError(error_msg.str());
            }
            else if(eErrClass == CE_Warning)
            {
                error_msg << "Gdal warning: " << err_no << ". " << msg << std::endl;
                writeWarning(error_msg.str());
            }
#ifdef DEBUG
            else
            {
                writeLog(10, error_msg.str());
            }
#endif // DEBUG
        }
    }
}