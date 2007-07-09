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

#ifndef _UTIL_H
#define _UTIL_H

#include <Python.h>
#include "glpk.h"

#if PY_MAJOR_VERSION == 2
#if PY_MINOR_VERSION < 5
typedef int Py_ssize_t;
typedef inquiry lenfunc;
typedef intargfunc ssizeargfunc;
typedef intobjargproc ssizeobjargproc;
#endif // PY_MINOR_VERSION < 5
#if PY_MINOR_VERSION < 4
#define Py_VISIT(o) do { if(o) { int v=visit((PyObject*)o,arg); \
			   if (v) return v; } } while (0)
#define Py_CLEAR(o) do { if(o) { PyObject *t=(PyObject*)o; (o)=NULL; \
			   Py_DECREF(t); } } while (0)
#define Py_RETURN_NONE return Py_INCREF(Py_None), Py_None
#endif // PY_MINOR_VERSION < 4
#endif // PY_MAJOR_VERSION == 2

#define GLPK_VERSION(MA, MI) ((MA)<GLPK_MAJOR_VERSION || (MA)==GLPK_MAJOR_VERSION && (MI)<=GLPK_MINOR_VERSION)


/* If "ob" is iterable. Returns 0 on failure with an appropriate
   exception set, 1 on success.
*/

int util_extract_if(PyObject *ob, PyObject *barcol,
		    int *len, int **ind, double **val);
int util_extract_iif(PyObject *ob, PyObject*lp,
		     int *len, int **ind1, int **ind2, double **val);

/* Generic type addition utility.  Returns 0 if the type was
   successfully added to the module, -1 if not. */
int util_add_type(PyObject *module, PyTypeObject *type);

/* On first call, stderr and stdout are redirected to /dev/null.  On
   second call, redirected back to old FD, and so on, and so forth.
   Sort of a toggle switch.  Returns 0 if this resulted in shutting
   stderr off, 1 if stderr is now "visible" again.
   
   In the event that the lp parameter LPX_K_MSGLEV has a value > 0,
   this function will do nothing. */
int util_toggle_output(LPX*lp);

#endif // # _UTIL_H
