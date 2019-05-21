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


%define FIBERRULESDATA_EXTENDER(NAME)

%ignore FiberRulesData::NAME;
%extend FiberRulesData
{
  
  PyObject* _ ## NAME ## _get()
  {
    // Get dimensions the array should have
    npy_intp adims[1] = {self->NAME.size()};
    
    // Create array
    return PyArray_SimpleNewFromData(1, &adims[0], NPY_DOUBLE,
  				     (char *)(&(self->NAME[0])));
  }
  
  // Extend Python interface
  %pythoncode%{
    def _ ## NAME ## _set(self, value):
        import numpy as np
        np_array = self._ ## NAME ## _get()
        if isinstance(value, np.ndarray):
            value = np.reshape(value, np_array.shape)
        np_array[:] = value

    NAME = property(_ ## NAME ## _get, _ ## NAME ## _set)
  %}
}
%enddef

%extend FiberRulesData
{
  %pythoncode%{
    def __setitem__(self, name, value):
        setattr(self, name, value)

    def __getitem__(self, name):
        return getattr(self, name)
  %}
  
}

FIBERRULESDATA_EXTENDER(lv_scalar)
FIBERRULESDATA_EXTENDER(rv_scalar)
FIBERRULESDATA_EXTENDER(epi_scalar)
FIBERRULESDATA_EXTENDER(apex_gradient)
FIBERRULESDATA_EXTENDER(lv_gradient)
FIBERRULESDATA_EXTENDER(rv_gradient)
FIBERRULESDATA_EXTENDER(epi_gradient)
FIBERRULESDATA_EXTENDER(fiber_sheet_tensor)

// Instantiate std::vector<Element> template
%template(ElementVector) std::vector<Element>;

// Extend Element
%extend Element 
{
  uint __getitem__(uint index)
  {
    if (index >= self->indices.size())
      throw std::runtime_error("Index out of range");
    
    return self->indices[index];
  }

  void __setitem__(uint index, double value)
  {
    if (index >= self->indices.size())
      throw std::runtime_error("Index out of range");
    
    self->indices[index] = value;
  }

  PyObject* set_indices(PyObject* indices)
  {
    
    // Check type
    if (!PyArray_Check(indices))
    {
      PyErr_SetString(PyExc_TypeError,"contiguous numpy array of 'size_t' expected. Make sure that the numpy array is contiguous, and uses dtype=np.uintp.");
      return NULL;
    }

    // Get PyArrayObject
    PyArrayObject *xa = reinterpret_cast<PyArrayObject*>(indices);

    // Check data type
    if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_UINTP))
    {
      PyErr_SetString(PyExc_TypeError, "contiguous numpy array of 'size_t' expected. Make sure that the numpy array is contiguous, and uses dtype=np.uintp.");
      return NULL;
    }

    // Get size
    const uint size = PyArray_SIZE(xa);

    //printf("Size: %d: [", size);

    // Resize vector if needed
    if (size != self->indices.size())
      self->indices.resize(size);

    // Get data
    const std::size_t* data = static_cast<std::size_t*>(PyArray_DATA(xa));
    
    // Set data
    for (uint i=0; i<size; i++)
    {
      //printf("%d, ", data[i]);
      self->indices[i] = data[i];
    }
    //printf("]\n");
    
    // Return
    Py_INCREF(Py_None);
    return Py_None;
  }
}

// Make readElemFile return an ElementVector
%pythoncode%{
  def readElemFile(filename):
      """
      Read a carp element file
      """
      elem_vector = ElementVector()
      _readElemFile(filename, elem_vector)
      return elem_vector
  %}
