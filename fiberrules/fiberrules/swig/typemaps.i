/* -*- C -*- */
// Copyright (C) 2012 Johan Hake
//
// This file is part of FIBERRULES.
//
// FIBERRULES is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// FIBERRULES is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with FIBERRULES. If not, see <http://www.gnu.org/licenses/>.


//-----------------------------------------------------------------------------
// Argout typemap for the reading of ptr
//-----------------------------------------------------------------------------
%typemap(in, numinputs=0) (uint& numPts, Real*& data) (uint numPts_tmp, Real* data_tmp)
{
  $1 = &numPts_tmp;
  $2 = &data_tmp;
}

%typemap(argout) (uint& numPts, Real*& data)
{
  // Create 2 dimensional NumPy array
  npy_intp adims[2] = {*$1,3};
  PyArrayObject* array = reinterpret_cast<PyArrayObject*>(\
			  PyArray_SimpleNewFromData(2, adims, NPY_DOUBLE,
						    (char *)(*$2)));
  // Let array own data
#if NUMPY_VERSION_MINOR >= 7
  PyArray_ENABLEFLAGS(array, NPY_ARRAY_OWNDATA);
#else
  array->flags |= NPY_OWNDATA;
#endif

  // Append the output to $result
  %append_output(PyArray_Return(array));

}

//-----------------------------------------------------------------------------
// In typemap for pts and data arrays
//-----------------------------------------------------------------------------
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (const uint numPts, const Real* pts)
{
    $1 = PyArray_Check($input) ? 1 : 0;
}

%typemap(in) (const uint numPts, const Real* pts)
{

  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of 'double' expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");

  // Set output
  $1 = PyArray_SIZE(xa);
  $2 = (double *)PyArray_DATA(xa);
}

// Copy typemap
%typemap(in) (const uint numData, const Real* data) = (const uint numPts, const Real* pts);
%typemap(in) (const uint numData, const Real* apexScalar) = (const uint numPts, const Real* pts);
%typemap(in) (const uint numNoise, const Real* noise) = (const uint numPts, const Real* pts);
%typemap(in) (const uint numNoise2, const Real* noise2) = (const uint numPts, const Real* pts);


// Prevent setting numScalarData
%typemap(in) uint numScalarData
{
  SWIG_exception(SWIG_AttributeError, "$1_name attribute is not writable.");
}

//-----------------------------------------------------------------------------
// Macro for defining an in typemap for a [const] Real M[N] of Real
// The typemaps takes a NumPy array of that primitive
//
// TYPE       : The primitive type
// TYPE_UPPER : The SWIG specific name of the type used in the array type checks
//              values SWIG use: INT32 for integer, DOUBLE for double aso.
// ARG_NAME   : The name of the argument
// NUMPY_TYPE : The type of the NumPy array that will be returned
// TYPE_NAME  : The descriptor name of the pointer type
// DESCR      : The char descriptor of the NumPy type
// SIZE       : The size 
// CONST      : Const or not
//-----------------------------------------------------------------------------
%define FIXED_SIZED_ARRAY_IN_TYPEMAP(TYPE, TYPE_UPPER, ARG_NAME, \
				     NUMPY_TYPE, TYPE_NAME, DESCR, SIZE, CONST, )

// The typecheck
%typecheck(SWIG_TYPECHECK_ ## TYPE_UPPER ## _ARRAY) CONST TYPE ARG_NAME[SIZE]
{
    $1 = PyArray_Check($input) ? 1 : 0;
}

%typemap(in) CONST TYPE ARG_NAME[SIZE]
{

  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NUMPY_TYPE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of 'TYPE_NAME' expected."
           " Make sure the numpy array is contiguous, and uses dtype=DESCR.");

  // Check size of passed array
  if ( PyArray_SIZE(xa) != SIZE )
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of size SIZE.");
  
  $1 = (TYPE *)PyArray_DATA(xa);
}

%enddef

// Instantiate typemaps
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, M, NPY_DOUBLE, double, np.float_, 9, )
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, M, NPY_DOUBLE, double, np.float_, 9, const)
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, R, NPY_DOUBLE, double, np.float_, 9, )
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, apex, NPY_DOUBLE, double, np.float_, 3, const)
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, surface_inward, NPY_DOUBLE, double, np.float_, 3, const)
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, center_coor, NPY_DOUBLE, double, np.float_, 3, )
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, gradient, NPY_DOUBLE, double, np.float_, 3, )
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, lv_gradient, NPY_DOUBLE, double, np.float_, 3, const)
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, rv_gradient, NPY_DOUBLE, double, np.float_, 3, const)
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, epi_gradient, NPY_DOUBLE, double, np.float_, 3, const)
FIXED_SIZED_ARRAY_IN_TYPEMAP(Real, DOUBLE, apex_gradient, NPY_DOUBLE, double, np.float_, 3, const)
