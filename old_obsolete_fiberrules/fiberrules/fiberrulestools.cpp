#include "types.h"
#include "fiberrulestools.h"

#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <stdexcept>
#include <cassert>

#ifdef HAS_ARMADILLO
#include <armadillo>
#endif

//-----------------------------------------------------------------------------
void print_vec(const char* name, const Real v[3]) 
{
  printf("%s = [%g %g %g]\n", name, v[0], v[1], v[2]);
}
//-----------------------------------------------------------------------------
void print_quaternion(const char* name, const Quaternion q) 
{
  printf("%s = [%g %g %g %g]\n", name, q.s, q.i, q.j, q.k);
}
//-----------------------------------------------------------------------------
void print_matrix(const char* name, Real M[9]) 
{
  printf("%s = [\n", name);
  printf(" %g %g %g\n", M[0*3+0], M[0*3+1], M[0*3+2]);
  printf(" %g %g %g\n", M[1*3+0], M[1*3+1], M[1*3+2]);
  printf(" %g %g %g\n", M[2*3+0], M[2*3+1], M[2*3+2]);
  printf("]\n");
}
//-----------------------------------------------------------------------------
Real my_dot(const Real a[3], const Real b[3]) 
{
  Real sum = 0;
  for (int ii=0; ii < 3; ii++) 
    sum += a[ii] * b[ii];

  return sum;
}
//-----------------------------------------------------------------------------
Real my_norm(const Real v[3])
{
  return std::sqrt(my_dot(v, v));
}
//-----------------------------------------------------------------------------
void my_normalize(Real v[3])
{
  Real mag = my_norm(v);
  for (int ii=0; ii < 3; ii++)
    v[ii] /= mag;
}
//-----------------------------------------------------------------------------
void my_cross(const Real a[3], const Real b[3], Real c[3]) 
{
  c[0] = a[1]*b[2] - b[1]*a[2];
  c[1] = b[0]*a[2] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
} 
//-----------------------------------------------------------------------------
Real quaternionDot(const Quaternion& a, const Quaternion& b) 
{
  return a.s*b.s + a.i*b.i + a.j*b.j + a.k*b.k;
}
//-----------------------------------------------------------------------------
Quaternion quaternionMult(const Quaternion& a, const Quaternion& b) 
{
  Quaternion result;
  result.s = a.s*b.s - a.i*b.i - a.j*b.j - a.k*b.k;
  result.i = a.s*b.i + a.i*b.s + a.j*b.k - a.k*b.j;
  result.j = a.s*b.j - a.i*b.k + a.j*b.s + a.k*b.i;
  result.k = a.s*b.k + a.i*b.j - a.j*b.i + a.k*b.s;
  return result;
}
//-----------------------------------------------------------------------------
Quaternion bislerp(const Real a_coeff, const Quaternion& a, 
		   const Real b_coeff, const Quaternion& b) 
{

  const Real tol = 1e-12;
  Quaternion R[4];

  // HAKE: What are these 4 Quaternions representing

  // R0 = a
  // Original direction of a
  R[0] = a;

  // R1 = i*a
  // a
  R[1].s = -a.i;
  R[1].i =  a.s;
  R[1].j = -a.k;
  R[1].k =  a.j;

  // R2 = j*a
  R[2].s = -a.j;
  R[2].i =  a.k;
  R[2].j =  a.s;
  R[2].k = -a.i;

  // R3 = k*a
  R[3].s = -a.k;
  R[3].i = -a.j;
  R[3].j =  a.i;
  R[3].k =  a.s;

  // Find the quaternion from a candidates that is closest to b.
  Real max_dot_p = 0;
  int max_quat_index = -1;
  for (int ii=0; ii < 4; ii++) 
  {
    // Find some sort of direction of Quaternion
    Real dot_p = quaternionDot(R[ii], b);

    // Switch direction if negative
    if (dot_p < 0) 
    {
      dot_p = -dot_p;
      R[ii].s *= -1;
      R[ii].i *= -1;
      R[ii].j *= -1;
      R[ii].k *= -1;
    }

    // If dot product is larger we are closer 
    if (dot_p > max_dot_p) 
    {
      max_quat_index=ii;
      max_dot_p = dot_p;
    }
  }

  // The closest Quaternion
  const Quaternion& max_quat = R[max_quat_index];
  Quaternion result;

  // Check if Quaterions are very close (theta == sin_theta == 0)
  if (max_dot_p > 1 - tol) 
  {
    // If so just pick one.
    result = b;
  } 
  else 
  {
  
    // Compute the interpolation
    const Real theta = acos(max_dot_p);
    const Real sin_theta = std::sqrt(1-max_dot_p*max_dot_p);
    
    const Real adjusted_a_coeff = sin(theta*a_coeff)/sin_theta;
    const Real adjusted_b_coeff = sin(theta*b_coeff)/sin_theta;

    result.s = max_quat.s*adjusted_a_coeff + b.s*adjusted_b_coeff;
    result.i = max_quat.i*adjusted_a_coeff + b.i*adjusted_b_coeff;
    result.j = max_quat.j*adjusted_a_coeff + b.j*adjusted_b_coeff;
    result.k = max_quat.k*adjusted_a_coeff + b.k*adjusted_b_coeff;

  }

  // Return the interpolated result
  return result;
}
//-----------------------------------------------------------------------------
void quaternionToRotation(const Quaternion& q, Real M[9]) 
{

  // Help variables
  const Real w=q.s;
  const Real x=q.i;
  const Real y=q.j;
  const Real z=q.k;

  // Compute unit vectors
  M[0*3+0] = w*w+x*x-y*y-z*z; 
  M[0*3+1] = 2*x*y+2*w*z;
  M[0*3+2] = 2*x*z-2*w*y;

  M[1*3+0] = 2*x*y-2*w*z;
  M[1*3+1] = w*w-x*x+y*y-z*z;
  M[1*3+2] = 2*y*z+2*w*x;

  M[2*3+0] = 2*x*z+2*w*y;
  M[2*3+1] = 2*y*z-2*w*x;
  M[2*3+2] = w*w-x*x-y*y+z*z;
}
//-----------------------------------------------------------------------------
Quaternion rotationToQuaternion(const Real M[9]) 
{
  Quaternion q;
  const Real tol = 1e-12;
  

  // If real part is larger than 0
  const Real w2 = (1+M[0*3+0]+M[1*3+1]+M[2*3+2])/4;
  if (w2 > tol) 
  {
    q.s = std::sqrt(w2);
    q.i = (M[1*3+2]-M[2*3+1])/4/q.s;
    q.j = (M[2*3+0]-M[0*3+2])/4/q.s;
    q.k = (M[0*3+1]-M[1*3+0])/4/q.s;
  } 

  // If real part is zero
  else 
  {
    q.s = 0;

    // Check if first imaginary value (i) is zero
    const Real x2 = -(M[1*3+1]+M[2*3+2])/2;
    if (x2 > tol) 
    {
      q.i = std::sqrt(x2);
      q.j = M[0*3+1]/2/q.i;
      q.k = M[0*3+2]/2/q.i;
    }
    else 
    {
      q.i = 0;

      // Check if second imaginary value (j) is zero
      const Real y2 = (1-M[2*3+2])/2;
      if (y2 > tol) 
      {
        q.j = std::sqrt(y2);
        q.k = M[1*3+2]/2/q.j;
      } 
      else 
      {
        q.j = 0;
        q.k = 1;
      }
    }
  }

  return q;
}
//-----------------------------------------------------------------------------
Quaternion axisSystem(const Real apex[3], const Real surface_inward[3]) 
{

  // Compute the reference frame
  Real M[9];

  // Copy first component
  Real e0[3];
  for (int ii=0; ii < 3; ii++) 
    e0[ii] = apex[ii];

  // Normalize first component
  my_normalize(e0);
  for (int ii=0; ii < 3; ii++)
    M[ii*3+0] = e0[ii];

  // Create third component
  // Orthagonolizing the surface_inward vector to the apex unit vector 
  Real e2[3];
  for (int ii=0; ii < 3; ii++)
    e2[ii] = surface_inward[ii];
  my_normalize(e2);

  //Real dotp = my_dot(e0, surface_inward);
  Real dotp = my_dot(e0, e2);
  for (int ii=0; ii < 3; ii++)
    e2[ii] -= dotp*e0[ii];
  //e2[ii] = surface_inward[ii] - dotp*e0[ii];
  
  my_normalize(e2);
  for (int ii=0; ii < 3; ii++)
    M[ii*3+2] = e2[ii];
  
  // Create second component by cross product of first two components
  Real e1[3];
  my_cross(e2, e0, e1);
  my_normalize(e1);
  for (int ii=0; ii < 3; ii++)
    M[ii*3+1] = e1[ii];

  return rotationToQuaternion(M);
}
//-----------------------------------------------------------------------------
Quaternion rotateAxis(const Quaternion q, const Real alpha, const Real beta) 
{
  // Doing some cool magic to rotate the Quaternion
  const Quaternion a = {cos(-alpha/2), 0.0, 0.0, sin(-alpha/2)};
  const Quaternion b = {cos(-beta/2), sin(-beta/2), 0.0, 0.0};
  
  return quaternionMult(quaternionMult(b, a), q);
}
//-----------------------------------------------------------------------------
Real dataCentroid(const Element& element, const uint, 
		  const Real* phi_data)
{
  Real centroid = 0;
  

  // If pyramid
  if (element.type == "Py") 
  {
    // Average the values from all the vertices
    for (int ii=0; ii < 4; ii++)
      centroid += 3*phi_data[element.indices[ii]];

    centroid += 4*phi_data[element.indices[4]];
    centroid /= 16;

  } 
  else 
  {
    // Average the values from all the vertices
    for (uint ii=0; ii < element.indices.size(); ii++)
      centroid += phi_data[element.indices[ii]];

    centroid /= element.indices.size();
  }

  return centroid;
}
//-----------------------------------------------------------------------------
void computeCenterCoordinate(const Element& this_element,
			     const uint numPts, const Real* pts, 
			     Real center_coor[3])
{
  // Find the center of this element;
  for (int idim=0; idim < 3; idim++)
    center_coor[idim] = dataCentroid(this_element, numPts, &pts[idim*numPts/3]);
}
//-----------------------------------------------------------------------------
#ifdef HAS_ARMADILLO
void computeGradient(const Element& const_element, const uint numPts, const Real* pts,
		     const uint numData, const Real* data, Real gradient[3]) 
{

  // De constify the element. No real change is done...
  Element& this_element = const_cast<Element&>(const_element);

  if (numData != numPts/3)
    throw std::runtime_error("Num coordinates and data points have to be the same.");

  const uint numPts_loc = numPts/3;

  // Find the center of this element;
  Real center[3];
  computeCenterCoordinate(this_element, numPts, pts, center);
  
  // Find the center data value of this element
  const Real centroid = dataCentroid(this_element, numPts_loc, data);

  // Find the std deviation of the distance from the center.
  Real std_dev_distance = 0.0;
  for (uint ii=0; ii < this_element.indices.size(); ii++) 
  {
    Real v[3];
    for (int idim=0; idim < 3; idim++)
      v[idim] = pts[idim*numPts_loc + this_element.indices[ii]] - center[ii];
    
    std_dev_distance += my_dot(v,v);
  }
  
  std_dev_distance = std::sqrt(std_dev_distance);

  // Find the std dev of the data values.
  Real std_dev_values = 0;
  for (uint ii=0; ii < this_element.indices.size(); ii++)
    std_dev_values += (data[this_element.indices[ii]] - centroid) * \
      (data[this_element.indices[ii]] - centroid);
  
  std_dev_values = std::sqrt(std_dev_values);
  
  // DEBUG:
  if (std_dev_values < 1e-8)
    std_dev_values = 1.0e-8;

  /* Find c,x,y,z in least squares system
   * c + dot(pt(i)-center,[x,y,z]) = data(i) 
   */

  // Check that we have a positive structure
  double det = 1.0;
  uint tmp_ind = 0;
  if (this_element.type == "Tt")
  {
    // Tetrahedron
    
    // Fill coord vecs
    std::vector<arma::vec> coords(4);
    for (uint ii=0; ii < 4; ii++) 
    {
      coords[ii].set_size(3);
      for (uint idim=0; idim < 3; idim++)
	coords[ii](idim) = pts[idim*numPts_loc + this_element.indices[ii]];
    }

    // Calculate det([a-d, b-d, c-d]) == (a-d)*((b-d)X(c-d))
    det = arma::dot(coords[0]-coords[3],				\
		    arma::cross(coords[1]-coords[3], coords[2]-coords[3]));

    // If negative determinant
    if (det<0.0)
    {
      tmp_ind = this_element.indices[3];
      this_element.indices[3] = this_element.indices[2];
      this_element.indices[2] = tmp_ind;
    }
  }

  // Construct least square system
  arma::mat A(this_element.indices.size(), 4);
  arma::vec b(this_element.indices.size());

  for (uint ii=0; ii < this_element.indices.size(); ii++) 
  {
    b(ii) = data[this_element.indices[ii]] - centroid;
    A(ii, 0) = 1;
    for (int idim=0; idim < 3; idim++)
      A(ii, idim+1) = pts[idim*numPts_loc + this_element.indices[ii]] - center[idim];
  }
  
  //std::cout << "b:" << std::endl;
  //std::cout << b << std::endl;
  //
  //std::cout << "Constraint matrix:" << std::endl;
  //std::cout << A << std::endl;

  // Solve least squares system
  arma::Col<double> x = arma::solve(A, b);

  // Return the result.
  // Remember to normalize w/r/t to the mesh size and mesh units.
  for (int idim=0; idim < 3; idim++) 
    gradient[idim] = x(idim+1)*std_dev_distance/std_dev_values;

  // Change back index changes
  if (det<0.0)
  {
    this_element.indices[2] = this_element.indices[3];
    this_element.indices[3] = tmp_ind;
  }
}
//-----------------------------------------------------------------------------
void computeGradientsAndCenterValues(const FiberRulesData& vertex_data,
				     const std::vector<Element>& elements,
				     const uint numPts, const Real* pts,
				     const uint numData, const Real* apexScalar,
				     FiberRulesData& center_data)
{
  // Iterate over the elements
  for (uint ielem=0; ielem < elements.size(); ++ielem)
  {
    const Element& element = elements[ielem];

    // Compute the average value at the center of an element
    center_data.lv_scalar[ielem] = dataCentroid(element, numData, &vertex_data.lv_scalar[0]);
    center_data.rv_scalar[ielem] = dataCentroid(element, numData, &vertex_data.rv_scalar[0]);
    center_data.epi_scalar[ielem] = dataCentroid(element, numData, &vertex_data.epi_scalar[0]);
    
    // Compute gradients at the center of an element
    computeGradient(element, numPts, pts, numData, apexScalar, \
		    &center_data.apex_gradient[ielem*3]);
    computeGradient(element, numPts, pts, numData, &vertex_data.lv_scalar[0], \
		    &center_data.lv_gradient[ielem*3]);
    computeGradient(element, numPts, pts, numData, &vertex_data.rv_scalar[0], \
		    &center_data.rv_gradient[ielem*3]);
    computeGradient(element, numPts, pts, numData, &vertex_data.epi_scalar[0], \
		    &center_data.epi_gradient[ielem*3]);

    // Find the center of this element;
    Real center[3];
    computeCenterCoordinate(element, numPts, pts, center);
  }
}
#endif
//-----------------------------------------------------------------------------
void computeSingleFiberSheetSystem(const Real lv_scalar, 
				   const Real rv_scalar, 
				   const Real epi_scalar, 
				   const Real apex_gradient[3],
				   const Real lv_gradient[3], 
				   const Real rv_gradient[3],
				   const Real epi_gradient[3], 
				   const Real fiber_rotation_epi, 
				   const Real fiber_rotation_endo, 
				   const Real sheet_rotation_epi, 
      				   const Real sheet_rotation_endo, 
                   		   Real R[9],
                   		   const Real  noisei, const Real noisei2)
{
  const Real tol = 1e-7;
  const Real grad_tol = 1e-7;
  
  // Compute the quaternion for each surface
  Quaternion lv_axis;
  if (lv_scalar >= tol)
    lv_axis = axisSystem(apex_gradient, lv_gradient);

  Quaternion rv_axis;
  if (rv_scalar >= tol)
    rv_axis = axisSystem(apex_gradient, rv_gradient);

  Quaternion epi_axis;
  if (epi_scalar >= tol) 
    epi_axis = axisSystem(apex_gradient, epi_gradient);
  
  // Do the interpolation
  Quaternion q;

  // On the epicardium
  if (lv_scalar < tol && rv_scalar < tol)
  {
    q = rotateAxis(epi_axis, fiber_rotation_epi*M_PI/180. + noisei*M_PI/180., sheet_rotation_epi*M_PI/180 + noisei2*M_PI/180.);
  }
  else 
  {
    // We are anywhere in the mid-wall or at the one of the endo surfaces
    const Real septum_scalar = lv_scalar + rv_scalar;
    const Real fiber_rotation_septum = (fiber_rotation_endo * lv_scalar/septum_scalar + \
					(180-fiber_rotation_endo)* \
					rv_scalar/septum_scalar) * M_PI/180.;
    const Real sheet_rotation_septum = (sheet_rotation_endo * lv_scalar/septum_scalar -
					sheet_rotation_endo * rv_scalar/septum_scalar) * M_PI/180.;
    
    // We are in the lv or septum wall
    Quaternion lv_q;
    if (lv_scalar >= tol)
    {
      /** In the event the gradient is too small, use a sigmoid to
	  produce a fiber in the apex-to-base fiber direction
      */
      const Real papilary = 1/(1 + std::pow(my_norm(lv_gradient)/grad_tol, 2));
      lv_q = rotateAxis(lv_axis, 0*papilary + fiber_rotation_septum*(1-papilary) + noisei*M_PI/180.,
			0*papilary + sheet_rotation_septum *(1-papilary) + noisei2*M_PI/180.);
    }
    
    // We are in the rv or septum wall
    Quaternion rv_q;
    if (rv_scalar >= tol) 
    {
      const Real papilary = 1/(1 + std::pow(my_norm(rv_gradient)/grad_tol, 2));
      rv_q = rotateAxis(rv_axis, 0 * papilary - fiber_rotation_septum * (1-papilary) + noisei*M_PI/180.,
			0*papilary + sheet_rotation_septum * (1 - papilary) + noisei2*M_PI/180.);
    }

    // We are only in the septum wall
    Quaternion septum_q;

    // We are on the lv endo wall 
    if (lv_scalar >= tol && rv_scalar < tol) 
      septum_q = lv_q;
    
    // We are on the rv endo wall 
    else if (lv_scalar < tol && rv_scalar >= tol) 
      septum_q = rv_q;
    
    // We are in the septum mid wall
    else 
      septum_q = bislerp(lv_scalar/septum_scalar, lv_q,
			 rv_scalar/septum_scalar, rv_q);
    
    // If we are indeed in the septum mid wall (No effect of the epi wall)
    if (epi_scalar < tol)
    {
      q = septum_q;
    }
    
    // We are close to the epi wall
    else 
    {
      
      const Real alpha_wall = (fiber_rotation_epi*epi_scalar +
			       (180 - fiber_rotation_endo) * septum_scalar) * M_PI/180.;
      const Real beta_wall = (sheet_rotation_epi * epi_scalar - 
			      sheet_rotation_endo*septum_scalar) * M_PI/180.;
      
      const Real papilary = 1/(1 + std::pow(my_norm(epi_gradient)/grad_tol, 2));
      Quaternion epi_q = rotateAxis(epi_axis, 0*papilary + alpha_wall*(1 - papilary) + noisei*M_PI/180.,
				    0*papilary + beta_wall*(1-papilary) + noisei2*M_PI/180.);
      
      q = bislerp(epi_scalar, epi_q, septum_scalar, septum_q);

    }
  }
  
  // Convert the quaternion into a rotation matrix
  quaternionToRotation(q, R);
}
//-----------------------------------------------------------------------------
void computeFiberSheetSystem(FiberRulesData& data,
                             const uint numNoise, const Real* noise,  const uint numNoise2,const Real* noise2)
{
  // Check data
  if (data.numScalarData == 0)
    throw std::runtime_error("data.numScalarData cannot be 0.");
    
  printf("Fiber rotation| epi: %.f; endo: %.f\nSheet rotation| epi: %.f; endo: %.f\n", \
	 data.fiber_rotation_epi, data.fiber_rotation_endo, \
	 data.sheet_rotation_epi, data.sheet_rotation_endo);
  
  // Iterate over the data and generate the single fiber sheet tensors
  for (uint i=0; i < data.numScalarData; i++ )
    computeSingleFiberSheetSystem(data.lv_scalar[i], 
				  data.rv_scalar[i], 
				  data.epi_scalar[i], 
				  &data.apex_gradient[3*i],
				  &data.lv_gradient[3*i], 
				  &data.rv_gradient[3*i],
				  &data.epi_gradient[3*i], 
				  data.fiber_rotation_epi, 
				  data.fiber_rotation_endo, 
				  data.sheet_rotation_epi, 
				  data.sheet_rotation_endo, 
                  &data.fiber_sheet_tensor[9*i],
                  noise[i], noise2[2]);
}
//-----------------------------------------------------------------------------
