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
// Include code to generate a __swigversion__ and a __fiberrulesversion__ 
// attributes, from defines during compile time, to the cpp module
//-----------------------------------------------------------------------------

%inline %{
int fiberrules_swigversion() { return SWIGVERSION; }
std::string fiberrules_version() { return FIBERRULES_VERSION; }
%}

%pythoncode %{
tmp = hex(fiberrules_swigversion())
__swigversion__ = "%d.%d.%d"%(tuple(map(int, [tmp[-5], tmp[-3], tmp[-2:]])))
__fiberrulesversion__ = fiberrules_version()
del tmp, fiberrules_swigversion, fiberrules_version
%}
