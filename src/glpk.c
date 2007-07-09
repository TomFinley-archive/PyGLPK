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

#include <Python.h>
#include "lp.h"

static PyMethodDef GLPKMethods[] = {
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initglpk(void) {
  PyObject *m;
  m = Py_InitModule("glpk", GLPKMethods);
  if (m==NULL) return;

  PyModule_AddObject(m, "version", Py_BuildValue
		     ("ii", GLPK_MAJOR_VERSION, GLPK_MINOR_VERSION));

  LPX_InitType(m);
}
