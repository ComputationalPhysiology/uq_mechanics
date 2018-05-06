#include "types.h"

#include <iostream>
#include <fstream>
#include <stdexcept>

//-----------------------------------------------------------------------------
void readPtsFile(const std::string filename, uint& numPts, 
		 Real*& data) 
{
  // Load the point file
  std::ifstream pts_file(filename.c_str());
  pts_file >> numPts;
  data = new double[3*numPts];
  
  // Read in the coordinate data
  if (pts_file)
  {
    // Iterate over the points
    for (uint ipt=0; ipt < numPts; ++ipt) 

      // Grab each coordinate
      for (uint idim=0; idim<3; ++idim)
	pts_file >> data[idim*numPts + ipt];
  }

  // Sanity check
  if (!pts_file) 
    throw std::runtime_error("error reading in pts file-- check the header?");
}
//-----------------------------------------------------------------------------
void readElemFile(const std::string filename, std::vector<Element>& elements)
{
  
  // Load the point filethis_element
  std::ifstream elem_file(filename.c_str());

  // Read number of elements and resize data
  uint numElem = 0;
  elem_file >> numElem;
  elements.resize(numElem);

  // Read in all elements
  std::string elemType;
  for (uint ielem=0; ielem < numElem; ++ielem) 
  {
    // Get current element
    Element& elem = elements[ielem];
    
    // Read type
    elem_file >> elemType;

    // Resize indices
    if (elemType == "Tt") 
      elem.indices.resize(4);
    else if (elemType == "Py")
      elem.indices.resize(5);
    else if (elemType == "Pr")
      elem.indices.resize(6);
    else if (elemType == "Oc")
      elem.indices.resize(6);
    else if (elemType == "Hx")
      elem.indices.resize(8);
    else
      throw std::runtime_error("Unknown element type: " + elemType);
    
    // Read vertex indices
    for (uint ii=0; ii < elem.indices.size(); ++ii)
      elem_file >> elem.indices[ii];

    // Read region info
    elem_file >> elem.region;
    
    // Assign the type
    elem.type = elemType;

  }

  // Sanity check
  if (!elem_file) 
    throw std::runtime_error("error reading in elem file-- check the header?");
}
//-----------------------------------------------------------------------------
void readData(const std::string filename, const uint numData , 
	      Real* data) 
{
  std::ifstream file(filename.c_str());
  for (uint ii=0; ii < numData; ii++)
    file >> data[ii];

  // Sanity check
  if (!file) 
    throw std::runtime_error("error reading in data file.");
}
//-----------------------------------------------------------------------------
