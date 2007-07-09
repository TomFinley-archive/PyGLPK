/**************************************************************************
Copyright (C) 2007 Thomas Finley, tomf@cs.cornell.edu

This file is part of PyGLPK.

PyGLPK is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

PyGLPK is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PyGLPK; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**************************************************************************/

#ifndef _LP_H
#define _LP_H

#include <Python.h>
#include "glpk.h"

#define LPX_Check(op) PyObject_TypeCheck(op, &LPXType)

typedef struct {
  PyObject_HEAD

  // The LPX C object.
  LPX *lp;
  // The rows and columns (bar collections).
  PyObject *rows, *cols;
  // The objective.
  PyObject *obj;
  // The parameters and statistics.
  PyObject *params;
  // The last type of solver that was used.  -1 for none, 0 for
  // simplex or exact, 1 for interior point, and 2 for integer or
  // intopt.
  int last_solver:4;
  PyObject *weakreflist; // Weak reference list.
} LPXObject;

extern PyTypeObject LPXType;

/* Init the type and related types it contains. 0 on success. */
int LPX_InitType(PyObject *module);

#endif // _LP_H
