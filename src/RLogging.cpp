//This file is part of rcoalescence project which is released under BSD-3 license.
//Visit https://opensource.org/licenses/BSD-3-Clause) for full license details.

#include "RLogging.h"
using namespace std;
void writeInfo(string message)
{
	if(logging_mode)
	{
		Rcpp::Rcout << message << flush;
	}
}

void writeWarning(string message)
{
	writeInfo(message);
}

void writeError(string message)
{
	writeInfo(message);
}

void writeCritical(string message)
{
	writeInfo(message);
}


