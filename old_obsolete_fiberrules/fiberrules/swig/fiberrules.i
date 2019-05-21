%module(package="fiberrules.cpp", directors="1") cpp
%{
#include "fiberrules/fiberrules.h"

// NumPy includes
#define PY_ARRAY_UNIQUE_SYMBOL PyFIBERRULES
#include <numpy/arrayobject.h>

%}

%init%{
import_array();
%}

//%include "fiberrules/swig/shared_ptr_classes.i"
%include "fiberrules/swig/typemaps.i"

// Global exceptions
%include <exception.i>
%include "fiberrules/swig/exceptions.i"

// Do not expand default arguments in C++ by generating an extra
// function in the SWIG layer. This reduces code bloat.
 //%feature("compactdefaultargs");

// STL SWIG string and vector class
%include <std_string.i>
%include <std_vector.i>

// Include information about swig version
%include "fiberrules/swig/version.i"

// Enable automatic docstring generation
// FIXME: Consider generate from C++
%feature("autodoc", "1");

// Include the interface with pre and post modifications
%include "fiberrules/swig/pre.i"
%include "fiberrules/types.h"
%include "fiberrules/fiberrulestools.h"
%include "fiberrules/carptools.h"
%include "fiberrules/swig/post.i"
