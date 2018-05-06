#ifndef TYPES_H_IS_INCLUDED
#define TYPES_H_IS_INCLUDED

#include <string>
#include <vector>

// Typdefs
typedef double Real;
typedef unsigned int uint;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

// The Quaternion data structure
struct Quaternion
{
  // Real value
  Real s;

  // Imaginary values
  Real i;
  Real j;
  Real k;

};

// Element data struct
struct Element 
{
  std::string type;
  std::vector<uint> indices;
  uint region;
};

// Data struct for easy communication 
class FiberRulesData 
{

public:
  
  // Constructor for default initialization
  FiberRulesData(uint numScalarData_): numScalarData(numScalarData_),
				       lv_scalar(numScalarData),
				       rv_scalar(numScalarData), 
				       epi_scalar(numScalarData),
				       apex_gradient(numScalarData*3),
				       lv_gradient(numScalarData*3),
				       rv_gradient(numScalarData*3),
				       epi_gradient(numScalarData*3),
				       fiber_sheet_tensor(numScalarData*9),
				       fiber_rotation_epi(50.0),
				       fiber_rotation_endo(40.0),
				       sheet_rotation_epi(65.0),
				       sheet_rotation_endo(25.0)
  {}
  
  // Empty destructor. 
  // Note: class does not owe its data
  ~FiberRulesData()
  { }

  // The number of data points
  uint numScalarData;
  
  // Scalar data
  std::vector<Real> lv_scalar;
  std::vector<Real> rv_scalar; 
  std::vector<Real> epi_scalar;

  // Gradient data
  std::vector<Real> apex_gradient;
  std::vector<Real> lv_gradient;
  std::vector<Real> rv_gradient;
  std::vector<Real> epi_gradient;

  // Fiber_sheet_tensor
  std::vector<Real> fiber_sheet_tensor;

  // Angle parameters
  Real fiber_rotation_epi;
  Real fiber_rotation_endo;
  Real sheet_rotation_epi;
  Real sheet_rotation_endo;

};


#endif
