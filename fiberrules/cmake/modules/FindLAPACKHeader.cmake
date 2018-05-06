# - Try to find LAPACK header clapack.h
# Once done this will define
#
#  LAPACKHEADER_FOUND  - system has LAPACK
#  LAPACK_INCLUDE_DIRS - include directories for LAPACK
#

#=============================================================================
# Copyright (C) 2012 Johan Hake
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

message(STATUS "Checking for package 'LAPACK Header'")

# Check for header file
find_path(LAPACK_INCLUDE_DIRS clapack.h
 HINTS ${LAPACK_DIR}/include $ENV{LAPACK_DIR}/include 
 ${ATLAS_DIR}/include $ENV{ATLAS_DIR}/include 
 ${LAPACK_DIR}/atlas/include $ENV{LAPACK_DIR}/atlas/include
 DOC "Directory where the LAPACK header is located"
 )
mark_as_advanced(LAPACK_INCLUDE_DIRS)

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKHEADER
  "LAPACK C header could not be found. Be sure to set LAPACK_DIR or ATLAS_DIR."
  LAPACK_INCLUDE_DIRS)
