                                 PyGLPK Readme

   Copyright (c) 2007, Thomas W. Finley

Overview

   PyGLPK is a Python module which encapsulates the functionality of the GNU
   Linear Programming Kit (GLPK). The GLPK allows one to specify linear
   programs (LPs) and mixed integer programs (MIPs), and to solve them with
   either simplex, interior-point, or branch-and-bound algorithms. The goal
   of PyGLPK is to give one access to all documented functionality of GLPK.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
   for more details.

Availability

   To get the lastest version, see:
   http://www.cs.cornell.edu/~tomf/pyglpk/

Documentation

   The HTML documentation included with the release in the directory html
   contains information on building, testing, installation and documentation
   of all features of the module.

Building and Installing

   Building this module requires that the user have installed the GLPK 4.4 or
   later and GMP libraries. The module builds and appears to work on my
   simple test files in Python 2.3, 2.4, and 2.5. Earlier versions of Python
   will not work.

   Ideally, the following will work:

     * make
     * make test
     * make install

   See the HTML documentation for trouble shooting information.

Bugs and Commentary

   Please send information on issues of usage to Thomas Finley at
   tomf@cs.cornell.edu .

     ----------------------------------------------------------------------

   Thomas Finley, 2007
