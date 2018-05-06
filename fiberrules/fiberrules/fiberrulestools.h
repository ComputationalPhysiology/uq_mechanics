#ifndef FIBERRULESTOOLS_H_IS_INCLUDED
#define FIBERRULESTOOLS_H_IS_INCLUDED

#include "types.h"

// Dot product of two Quaternions 
Real quaternionDot(const Quaternion& a, const Quaternion& b);

// Multiplication of two Quaternions 
Quaternion quaternionMult(const Quaternion& a, const Quaternion& b);

// Interpolation of two Quaternions 
Quaternion bislerp(const Real a_coeff, const Quaternion& a, 
		   const Real b_coeff, const Quaternion& b);

// Quaternion to rotation tensor
void quaternionToRotation(const Quaternion& q, Real M[9]);

// Rotation tensor to Quaternion
Quaternion rotationToQuaternion(const Real M[9]);

// Create axissystem (as a Quaternion)
Quaternion axisSystem(const Real apex[3], const Real surface_inward[3]);

// Rotate an axis given as a Quaternion, with the two angles alpha and beta
// FIXME: What do these angle determines?
Quaternion rotateAxis(const Quaternion q, const Real alpha, const Real beta);

// Average scalar data value (no geometric average)
Real dataCentroid(const Element& element, const uint numData, const Real* phi_data); 

// Compute the coordinate of the center of an element
void computeCenterCoordinate(const Element& this_element,
			     const uint numPts, const Real* pts,
			     Real center_coor[3]);

#ifdef HAS_ARMADILLO
// Compute the gradient of a scalar field
void computeGradient(const Element& this_element,
                     const uint numPts, const Real* pts,
                     const uint numData, const Real* data,
                     Real gradient[3]);

// Compute gradients and center values from element, pts and scalar data
void computeGradientsAndCenterValues(const FiberRulesData& vertex_data,
				     const std::vector<Element>& elements,
				     const uint numPts, const Real* pts,
				     const uint numData, const Real* apexScalar,
				     FiberRulesData& center_data);
#endif

// Compute the fiber sheet tensor given the scalar solutions of three
// laplacians together with four gradients fields
void computeSingleFiberSheetSystem(const Real lv_center, const Real rv_center, 
				   const Real epi_center, const Real apex_gradient[3],
				   const Real lv_gradient[3], const Real rv_gradient[3],
				   const Real epi_gradient[3], const Real alpha_epi_arg, 
				   const Real alpha_endo_arg, const Real beta_epi_arg, 
                   const Real beta_endo_arg, Real R[9],
                   const Real noise, const Real noise2);

// Compute the fiber sheet tensor given the scalar solutions of three
// laplacians together with four gradients fields
void computeFiberSheetSystem(FiberRulesData& data,
                             const uint numNoise, const Real* noise, const uint numNoise2,const Real* noise2);

#endif
