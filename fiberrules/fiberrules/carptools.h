#ifndef CARPTOOLS_H_IS_INCLUDED
#define CARPTOOLS_H_IS_INCLUDED

#include <string>
#include <vector>

#include "types.h"

// Read carp pts file
void readPtsFile(const std::string filename, uint& numPts, Real*& data);

// Read carp elem file
void readElemFile(const std::string filename, std::vector<Element>& elements);

// Read carp data file
void readData(const std::string filename, const uint numData, Real* data);

#endif
