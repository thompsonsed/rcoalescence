// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file CPLCustomHandlerNecsim.h
 * @brief Contains a custom CPL error handler.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#ifndef SPECIATIONCOUNTER_cplCustomHandlerNecsim_H
#define SPECIATIONCOUNTER_cplCustomHandlerNecsim_H
#ifdef with_gdal

#include <cpl_error.h>

/**
 * @brief Contains a custom CPLErrorHandler for reporting errors from GDAL imports.
 * @param eErrClass the error class (CE_None, CE_Debug, CE_Warning, CE_Failure or CE_Fatal)
 * @param err_no the error number to report
 * @param msg the message to report
 */
void cplNecsimCustomErrorHandler(CPLErr eErrClass, int err_no, const char *msg);

#endif // with_gdal

#endif //SPECIATIONCOUNTER_cplCustomHandlerNecsim_H
