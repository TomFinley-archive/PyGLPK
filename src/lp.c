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

#include "lp.h"
#include "structmember.h"
#include "barcol.h"
#include "obj.h"
#include "params.h"
#include "kkt.h"
#include "util.h"

#define LP (self->lp)

static int LPX_traverse(LPXObject *self, visitproc visit, void *arg) {
  Py_VISIT(self->rows);
  Py_VISIT(self->cols);
  Py_VISIT(self->obj);
  Py_VISIT(self->params);
  return 0;
}

static int LPX_clear(LPXObject *self) {
  if (self->weakreflist != NULL) {
    PyObject_ClearWeakRefs((PyObject*)self);
  }
  Py_CLEAR(self->rows);
  Py_CLEAR(self->cols);
  Py_CLEAR(self->obj);
  Py_CLEAR(self->params);
  return 0;
}

static void LPX_dealloc(LPXObject *self) {
  LPX_clear(self);
  if (LP) lpx_delete_prob(LP);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject * LPX_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  LPXObject *self;
  self = (LPXObject*) type->tp_alloc(type, 0);
  if (self != NULL) {
    self->lp = NULL;

    self->rows = NULL;
    self->cols = NULL;
    self->obj = NULL;
    self->params = NULL;
    self->weakreflist = NULL;

    self->last_solver = -1;
  }
  return (PyObject*)self;
}

static int LPX_init(LPXObject *self, PyObject *args, PyObject *kwds) {
  char *mps_n=NULL, *freemps_n=NULL, *cpxlp_n=NULL, *prob_n=NULL;
  PyObject *model_obj=NULL;
  static char *kwlist[] = {"gmp","mps","freemps","cpxlp","glp",NULL};
  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "|Ossss", kwlist,
       &model_obj, &mps_n, &freemps_n, &cpxlp_n, &prob_n)) {
    return -1;
  }
  int numargs = (mps_n?1:0)+(freemps_n?1:0)+(cpxlp_n?1:0)+(model_obj?1:0)+
    (prob_n?1:0);
  if (numargs>1) {
    PyErr_SetString(PyExc_TypeError, "cannot specify multiple data sources");
    return -1;
  }
  if (numargs==0) {
    // No arguments.  Create an empty problem.
    self->lp = lpx_create_prob();
  } else {
    // Some of these are pretty straightforward data reading routines.
    if (mps_n) {
      self->lp = lpx_read_mps(mps_n);
    } else if (freemps_n) {
#if GLPK_VERSION(4, 7)
      self->lp = lpx_read_freemps(freemps_n);
#else
      PyErr_SetString(PyExc_NotImplementedError,
		      "free MPS reading not present till GLPK 4.7"
#endif
    } else if (cpxlp_n) {
      self->lp = lpx_read_cpxlp(cpxlp_n);
    } else if (prob_n) {
#if GLPK_VERSION(4, 6) && !GLPK_VERSION(4, 15)
	self->lp = lpx_read_prob(prob_n);
#else
	PyErr_SetString
	  (PyExc_NotImplementedError,
	   "GNU LP reading not present till GLPK 4.6 or after 4.14");
	return -1;
#endif
    } else if (model_obj) {
      // This one can take a few possible values.
      char *model[] = {NULL,NULL,NULL};
      if (PyString_Check(model_obj)) {
	// Single string object.
	model[0] = PyString_AsString(model_obj);
	if (!model[0]) return -1;
      } else if (PyTuple_Check(model_obj)) {
	// Possibly module arguments.
	int i,size = PyTuple_Size(model_obj);
	if (size < -1) { return -1; }
	if (size >  3) { 
	  PyErr_SetString(PyExc_ValueError, "model tuple must have length<=3");
	  return -1; }
	for (i=0; i<size; ++i) {
	  PyObject *so = PyTuple_GET_ITEM(model_obj,i);
	  if (so==Py_None) continue;
	  model[i] = PyString_AsString(so);
	  if (model[i]==NULL) { return -1; }
	}
      } else {
	PyErr_SetString(PyExc_TypeError, "model arg must be string or tuple");
	return -1;
      }
      // Now, pass in that information.
      if (!model[0]) return -1;
      self->lp = lpx_read_model(model[0], model[1], model[2]);
    }
  }
  // Any of the methods above may have failed, so the LP would be null.
  if (LP == NULL) {
    PyErr_SetString(numargs?PyExc_RuntimeError:PyExc_MemoryError,
		    "could not create problem");
    return -1;
  }
  // Create those rows and cols and things.
  self->cols = (PyObject*)BarCol_New(self, 0);
  self->rows = (PyObject*)BarCol_New(self, 1);
  self->obj = (PyObject*)Obj_New(self);
  self->params = (PyObject*)Params_New(self);

  return 0;
}

static PyObject* LPX_Str(LPXObject *self) {
  // Returns a string representation of this object.
  return PyString_FromFormat
    ("<%s %d-by-%d at %p>", self->ob_type->tp_name,
     lpx_get_num_rows(LP), lpx_get_num_cols(LP), self);
}

PyObject *LPX_GetMatrix(LPXObject *self) {
  int row, numrows, listi, i, nnz, rownz;
  PyObject *retval;

  numrows = lpx_get_num_rows(LP);
  nnz = lpx_get_num_nz(LP);
  
  retval = PyList_New(nnz);
  if (nnz==0 || retval==NULL) return retval;

  // We don't really need this much memory, but, eh... 
  int *ind = (int*)calloc(nnz,sizeof(int));
  double *val = (double*)calloc(nnz,sizeof(double));

  listi = 0;
  for (row=1; row<=numrows; ++row) {
    rownz = lpx_get_mat_row(LP, row, ind-1, val-1);
    if (rownz==0) continue;
    for (i=0; i<rownz; ++i) {
      PyList_SET_ITEM(retval, listi++, Py_BuildValue
		      ("iid", row-1, ind[i]-1, val[i]));
    }
    // Continue to downscale these vectors, freeing memory in C even
    // as we use more memory in Python.
    nnz -= rownz;
    if (nnz) {
      ind = (int*)realloc(ind, nnz*sizeof(int));
      val = (double*)realloc(val, nnz*sizeof(double));
    }
  }
  free(ind);
  free(val);
  if (PyList_Sort(retval)) {
    Py_DECREF(retval);
    return NULL;
  }
  return retval;
}

int LPX_SetMatrix(LPXObject *self, PyObject *newvals) {
  int len1, len2, len, *ind1, *ind2;
  double*val;

  len1 = lpx_get_num_rows(LP);
  len2 = lpx_get_num_cols(LP);

  // Now, attempt to convert the input item.
  if (newvals == NULL || newvals == Py_None) {
    len = 0;
  } else if (!util_extract_iif(newvals, (PyObject*)self, &len,
			       &ind1, &ind2, &val)) {
    return -1;
  }

  // Input the stuff into the LP constraint matrix.
  lpx_load_matrix(LP, len, ind1-1, ind2-1, val-1);
  // Free the memory.
  if (len) {
    free(ind1);
    free(ind2);
    free(val);
  }
  return 0;
}

/****************** METHODS ***************/

/*static PyObject* LPX_OrderMatrix(LPXObject *self) {
  lpx_order_matrix(LP);
  Py_RETURN_NONE;
  }*/

static PyObject* LPX_Scale(LPXObject *self, PyObject*args) {
  PyObject *arg=NULL;
  int truth=1;
  PyArg_ParseTuple(args, "|O", &arg);
  if (arg!=NULL) {
    truth = PyObject_IsTrue(arg);
    if (truth==-1) return NULL;
  }
  (truth ? lpx_scale_prob : lpx_unscale_prob)(LP);
  Py_RETURN_NONE;
}

static PyObject* LPX_basis_std(LPXObject *self) {
  lpx_std_basis(LP); Py_RETURN_NONE; }

static PyObject* LPX_basis_adv(LPXObject *self) {
  lpx_adv_basis(LP); Py_RETURN_NONE; }

static PyObject* LPX_basis_cpx(LPXObject *self) {
#if GLPK_VERSION(4,10)
  lpx_cpx_basis(LP); Py_RETURN_NONE; 
#else
  PyErr_SetString(PyExc_NotImplementedError,
		  "Bixby's basis not present till GLPK 4.10");
  return NULL;
#endif
}

static PyObject* LPX_basis_read(LPXObject *self, PyObject *args) {
  char *bas_filename = NULL;
  if (!PyArg_ParseTuple(args, "s", &bas_filename)) {
    return NULL;
  }
  if (lpx_read_bas(LP, bas_filename)) {
    PyErr_SetString(PyExc_RuntimeError, "could not read basis file");
    return NULL;
  }
  Py_RETURN_NONE;
}

/**************** SOLVER METHODS **************/

static PyObject* solver_retval_to_message(int retval) {
  switch (retval) {
  case LPX_E_OK:	Py_RETURN_NONE;
  case LPX_E_FAULT:	return PyString_FromString("fault");
  case LPX_E_OBJLL:	return PyString_FromString("objll");
  case LPX_E_OBJUL:	return PyString_FromString("objul");
  case LPX_E_ITLIM:	return PyString_FromString("itlim");
  case LPX_E_TMLIM:	return PyString_FromString("tmlim");
  case LPX_E_SING:	return PyString_FromString("sing");

  case LPX_E_NOPFS:	return PyString_FromString("nopfs");
  case LPX_E_NODFS:	return PyString_FromString("nodfs");

  case LPX_E_NOFEAS:	return PyString_FromString("nofeas");
  case LPX_E_NOCONV:	return PyString_FromString("noconv");
  case LPX_E_INSTAB:	return PyString_FromString("instab");

  default:		return PyString_FromString("unknown?");
  }
}

static PyObject* LPX_solver_simplex(LPXObject *self) {
  int retval = lpx_simplex(LP);
  if (retval!=LPX_E_FAULT) self->last_solver = 0;
  return solver_retval_to_message(retval);
}

static PyObject* LPX_solver_exact(LPXObject *self) {
#if GLPK_VERSION(4, 13)  
  int retval = lpx_exact(LP);
  if (retval!=LPX_E_FAULT) self->last_solver = 0;
  return solver_retval_to_message(retval);
#else
  PyErr_SetString(PyExc_NotImplementedError,
		  "exact solver not present till GLPK 4.13");
  return NULL;
#endif
}

static PyObject* LPX_solver_interior(LPXObject *self) {
  util_toggle_output(LP);
  int retval = lpx_interior(LP);
  util_toggle_output(LP);

  if (retval!=LPX_E_FAULT) self->last_solver = 1;
  return solver_retval_to_message(retval);
}

static PyObject* LPX_solver_integer(LPXObject *self) {
  if (lpx_get_class(LP) != LPX_MIP) {
    PyErr_SetString(PyExc_RuntimeError,
		    "integer solver inapplicable to continuous problem");
    return NULL;
  }
  if (lpx_get_status(LP) != LPX_OPT) {
    PyErr_SetString(PyExc_RuntimeError, "integer solver requires "
		    "existing optimal basic solution");
    return NULL;
  }
  int retval = lpx_integer(LP);
  if (retval!=LPX_E_FAULT) self->last_solver = 2;
  return solver_retval_to_message(retval);
}

static PyObject* LPX_solver_intopt(LPXObject *self) {
#if GLPK_VERSION(4, 9)
  if (lpx_get_class(LP) != LPX_MIP) {
    PyErr_SetString(PyExc_RuntimeError,
		    "cannot apply integer solver to continuous problem");
    return NULL;
  }
  util_toggle_output(LP);
  int retval = lpx_intopt(LP);
  util_toggle_output(LP);
  if (retval!=LPX_E_FAULT) self->last_solver = 2;
  return solver_retval_to_message(retval);
#else
  PyErr_SetString(PyExc_NotImplementedError,
		  "intopt solver not present till GLPK 4.9");
  return NULL;
#endif
}

static KKTObject* LPX_kkt(LPXObject *self, PyObject *args) {
  // Cannot get for undefined primal or dual.
  if (lpx_get_prim_stat(LP)==LPX_P_UNDEF ||
      lpx_get_dual_stat(LP)==LPX_D_UNDEF) {
    PyErr_SetString(PyExc_RuntimeError, "cannot get KKT when primal or dual "
		    "basic solution undefined");
    return NULL;
  }
  // Check the Python arguments.
  int scaling = 0;
  PyObject *arg = NULL;
  if (!PyArg_ParseTuple(args, "|O", &arg)) return NULL;
  scaling = ((arg==NULL) ? 0 : PyObject_IsTrue(arg));
  if (scaling == -1) return NULL;
  // OK, all done with those checks.  Now get the KKT condition.
  KKTObject *kkt = KKT_New();
  if (!kkt) return NULL;
  lpx_check_kkt(LP, scaling, &(kkt->kkt));
  return kkt;
}

static KKTObject* LPX_kktint(LPXObject *self) {
#if GLPK_VERSION(4, 9)
  // Cannot get for an LP which is not a mixed integer program.
  if (lpx_get_class(LP) != LPX_MIP) {
    PyErr_SetString(PyExc_RuntimeError,
		    "cannot get int solution quality for non-MIP LP");
    return NULL;
  }
  // OK, all done with those checks.  Now get the KKT condition.
  KKTObject *kkt = KKT_New();
  if (!kkt) return NULL;
  lpx_check_int(LP, &(kkt->kkt));
  return kkt;
#else
  PyErr_SetString(PyExc_NotImplementedError,
		  "MIP feasibility checker not present till GLPK 4.9");
  return NULL;
#endif
}

static PyObject* LPX_write(LPXObject *self, PyObject *args, PyObject *keywds) {
  static char* kwlist[] = {"mps", "bas", "freemps", "cpxlp", "glp", "prob",
			   "sol", "sens_bnds", "ips", "mip", NULL};
  static int(*writers[])(LPX*,char*) = {
    lpx_write_mps, lpx_write_bas, 
#if GLPK_VERSION(4, 7)
    lpx_write_freemps, 
#else
    NULL,
#endif
    lpx_write_cpxlp,
#if GLPK_VERSION(4, 6) && !GLPK_VERSION(4, 15)
    lpx_write_prob,
#else
    NULL,
#endif
    lpx_print_prob, lpx_print_sol, lpx_print_sens_bnds, lpx_print_ips,
    lpx_print_mip};
  char* fnames[] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  int i;
  if (!PyArg_ParseTupleAndKeywords
      (args, keywds, "|ssssssssss", kwlist, fnames,fnames+1,fnames+2,fnames+3,
       fnames+4,fnames+5,fnames+6,fnames+7,fnames+8,fnames+9)) {
    return NULL;
  }
  for (i=0; i<10; ++i) {
    if (fnames[i]==NULL) continue;
    if (writers[i]==NULL) {
      PyErr_Format(PyExc_NotImplementedError,
		   "writer for '%s' format absent in this version of GLPK",
		   kwlist[i]);
      return NULL;
    }
    int retval = (writers[i])(LP, fnames[i]);
    if (retval==0) continue;
    PyErr_Format(PyExc_RuntimeError, "writer for '%s' failed to write to '%s'",
		 kwlist[i], fnames[i]);
    return NULL;
  }
  Py_RETURN_NONE;
}

/****************** GET-SET-ERS ***************/

static PyObject* LPX_getname(LPXObject *self, void *closure) {
  char *name = lpx_get_prob_name(LP);
  if (name==NULL) Py_RETURN_NONE;
  return PyString_FromString(name);
}
static int LPX_setname(LPXObject *self, PyObject *value, void *closure) {
  char *name;
  if (value==NULL || value==Py_None) {
    lpx_set_prob_name(LP, NULL);
    return 0;
  }
  name = PyString_AsString(value);
  if (name==NULL) return -1;
  if (PyString_Size(value) > 255) {
    PyErr_SetString(PyExc_ValueError, "name may be at most 255 chars");
    return -1;
  }
  lpx_set_prob_name(LP, name);
  return 0;
}

static PyObject* LPX_getobj(LPXObject *self, void *closure) {
  Py_INCREF(self->obj);
  return self->obj;
}

static PyObject* LPX_getnumnonzero(LPXObject *self, void *closure) {
  return PyInt_FromLong(lpx_get_num_nz(LP));
}

static PyObject* LPX_getmatrix(LPXObject *self, void *closure) {
  return LPX_GetMatrix(self);
}
static int LPX_setmatrix(LPXObject *self, PyObject *value, void *closure) {
  return LPX_SetMatrix(self, value);
}

static PyObject* status2string(int status) {
  switch (status) {
  case LPX_I_OPT:
  case LPX_T_OPT:
  case LPX_OPT:    return PyString_FromString("opt");
  case LPX_I_FEAS:
  case LPX_P_FEAS:
  case LPX_D_FEAS:
  case LPX_FEAS:   return PyString_FromString("feas");
  case LPX_P_INFEAS:
  case LPX_D_INFEAS:
  case LPX_INFEAS: return PyString_FromString("infeas");
  case LPX_I_NOFEAS:
  case LPX_P_NOFEAS:
  case LPX_D_NOFEAS:
  case LPX_NOFEAS: return PyString_FromString("nofeas");
  case LPX_UNBND:  return PyString_FromString("unbnd");
  case LPX_I_UNDEF:
  case LPX_T_UNDEF:
  case LPX_P_UNDEF:
  case LPX_D_UNDEF:
  case LPX_UNDEF:  return PyString_FromString("undef");
  default:         return PyString_FromString("unknown?");
  }
}
static PyObject* LPX_getstatus(LPXObject *self, void *closure) {
  int ls = self->last_solver, status;
  if      (ls <= 0) status=lpx_get_status(LP);
  else if (ls == 1) status=lpx_ipt_status(LP);
  else if (ls == 2) status=lpx_mip_status(LP);
  else {
    PyErr_SetString(PyExc_RuntimeError,
		    "bad internal state for last solver identifier");
    return NULL;
  }
  return status2string(status);
}
static PyObject* LPX_getspecstatus(LPXObject *self, int(*statfunc)(LPX*)) {
  return status2string(statfunc(LP));
}

static PyObject* LPX_getray(LPXObject *self, void *closure) {
  int ray = lpx_get_ray_info(LP), numrows;
  if (ray==0) Py_RETURN_NONE;
  numrows = lpx_get_num_rows(LP);
  ray--;
  if (ray < numrows) return PySequence_GetItem(self->rows, ray);
  return PySequence_GetItem(self->cols, ray - numrows);
}

static PyObject* LPX_getkind(LPXObject *self, void *closure) {
  PyObject *retval=NULL;
  switch (lpx_get_class(LP)) {
  case LPX_LP:  retval = (PyObject*)&PyFloat_Type; break;
  case LPX_MIP: retval = (PyObject*)&PyInt_Type;   break;
  default:
    PyErr_SetString(PyExc_RuntimeError, "unexpected problem kind encountered");
    return NULL;
  }
  Py_INCREF(retval);
  return retval;
}
static int LPX_setkind(LPXObject *self, PyObject *value, void *closure) {
  if (value==(PyObject*)&PyInt_Type) {
    lpx_set_class(LP, LPX_MIP);
    return 0;
  }
  if (value==(PyObject*)&PyFloat_Type) {
    lpx_set_class(LP, LPX_LP);
    return 0;
  }
  PyErr_SetString(PyExc_ValueError,"either the type float or int is required");
  return -1;
}

static PyObject* LPX_getnumint(LPXObject *self, void *closure) {
  return PyInt_FromLong(lpx_get_class(LP)==LPX_MIP?lpx_get_num_int(LP):0); }

static PyObject* LPX_getnumbin(LPXObject *self, void *closure) {
  return PyInt_FromLong(lpx_get_class(LP)==LPX_MIP?lpx_get_num_bin(LP):0); }

/****************** OBJECT DEFINITION *********/

int LPX_InitType(PyObject *module) {
  int retval;
  if ((retval=util_add_type(module, &LPXType))!=0) return retval;
  if ((retval=Obj_InitType(module))!=0) return retval;
  if ((retval=Params_InitType(module))!=0) return retval;
  if ((retval=BarCol_InitType(module))!=0) return retval;
  if ((retval=KKT_InitType(module))!=0) return retval;
  return 0;
}

static PyMemberDef LPX_members[] = {
  {"rows", T_OBJECT_EX, offsetof(LPXObject, rows), RO,
   "Row collection.  See the help on class BarCollection."},
  {"cols", T_OBJECT_EX, offsetof(LPXObject, cols), RO,
   "Column collection.  See the help on class BarCollection."},
  {"params", T_OBJECT_EX, offsetof(LPXObject, params), RO,
   "Control parameter collection.  See the help on class Params."},
  {NULL}
};

static PyGetSetDef LPX_getset[] = {
  {"name", (getter)LPX_getname, (setter)LPX_setname,
   "Problem name, or None if unset.", NULL},
  {"obj", (getter)LPX_getobj, (setter)NULL, 
   "Objective function object.", NULL},
  {"nnz", (getter)LPX_getnumnonzero, (setter)NULL,
   "Number of non-zero constraint coefficients.", NULL},
  {"matrix", (getter)LPX_getmatrix, (setter)LPX_setmatrix,
   "The constraint matrix as a list of three element (row index,\n"
   "column index, value) tuples across all non-zero elements of\n"
   "the constraint matrix.", NULL},
  // Solution status retrieval.
  {"status", (getter)LPX_getstatus, (setter)NULL,
   "The status of solution of the last solver.  This takes the\n"
   "form of a string with these possible values.\n\n"
   "opt    -- The solution is optimal.\n"
   "undef  -- The solution is undefined.\n"
   "feas   -- The solution is feasible, but not necessarily optimal.\n"
   "infeas -- The solution is infeasible.\n"
   "nofeas -- The problem has no feasible solution.\n"
   "unbnd  -- The problem has an unbounded solution.", NULL},
  {"status_s", (getter)LPX_getspecstatus, (setter)NULL,
   "The status of the simplex solver's solution.", (void*)lpx_get_status},
  {"status_i", (getter)LPX_getspecstatus, (setter)NULL,
   "The status of the interior point solver's solution.",
   (void*)lpx_ipt_status},
  {"status_m", (getter)LPX_getspecstatus, (setter)NULL,
   "The status of the MIP solver's solution.", (void*)lpx_mip_status},
  {"status_primal", (getter)LPX_getspecstatus, (setter)NULL,
   "The status of the primal solution of the simplex solver.\n"
   "Possible values are 'undef', 'feas', 'infeas', 'nofeas' in\n"
   "similar meaning to the .status attribute.",
   (void*)lpx_get_prim_stat},
  {"status_dual", (getter)LPX_getspecstatus, (setter)NULL,
   "The status of the dual solution of the simplex solver.\n"
   "Possible values are 'undef', 'feas', 'infeas', 'nofeas' in\n"
   "similar meaning to the .status attribute.",
   (void*)lpx_get_dual_stat},
  // Ray info.
  {"ray", (getter)LPX_getray, (setter)NULL,
   "A non-basic row or column the simplex solver has identified\n"
   "as causing primal unboundness, or None if no such variable\n"
   "has been identified.", NULL},
  // Setting for MIP.
  {"kind", (getter)LPX_getkind, (setter)LPX_setkind,
   "Either the type 'float' if this is a pure linear programming\n"
   "(LP) problem, or the type 'int' if this is a mixed integer\n"
   "programming (MIP) problem.", NULL},
  {"nint", (getter)LPX_getnumint, (setter)NULL,
   "The number of integer column variables.  Always 0 if this is\n"
   "not a mixed integer problem.", NULL},
  {"nbin", (getter)LPX_getnumbin, (setter)NULL,
   "The number of binary column variables, i.e., integer with 0\n"
   "to 1 bounds.  Always 0 if this is not a mixed integer problem.", NULL},
  {NULL}
};

static PyMethodDef LPX_methods[] = {
  /*{"order_matrix", (PyCFunction)LPX_OrderMatrix, METH_NOARGS,
   "order_matrix()\n\n"
   "Reorder internal rows and column entries in ascending index order."},*/
  {"scale", (PyCFunction)LPX_Scale, METH_VARARGS,
   "scale([toscale=True])\n\n"
   "Either scale or unscale depending on truth of the argument.\n"
   "Note that this only affects the internal state of the LP\n"
   "representation."},
  // Basis construction techniques for simplex solvers.
  {"std_basis", (PyCFunction)LPX_basis_std, METH_NOARGS,
   "std_basis()\n\n"
   "Construct the standard trivial inital basis for this LP."},
  {"adv_basis", (PyCFunction)LPX_basis_adv, METH_NOARGS,
   "adv_basis()\n\n"
   "Construct an advanced initial basis, triangular with as few\n"
   "variables as possible fixed."},
  {"cpx_basis", (PyCFunction)LPX_basis_cpx, METH_NOARGS,
   "cpx_basis()\n\n"
   "Construct an advanced Bixby basis.\n\n"
   "This basis construction method is described in:\n"
   "Robert E. Bixby. Implementing the Simplex Method: The Initial\n"
   "Basis.  ORSA Journal on Computing, Vol. 4, No. 3, 1992,\n"
   "pp. 267-84."},
  {"read_basis", (PyCFunction)LPX_basis_read, METH_VARARGS,
   "read_basis(filename)\n\n"
   "Reads an LP basis in the fixed MPS format from a given file."},
  // Solver routines.
  {"simplex", (PyCFunction)LPX_solver_simplex, METH_NOARGS,
   "simplex()\n\n"
   "Attempt to solve the problem using a simplex method.\n\n"
   "This returns None if the problem was successfully solved.\n"
   "Alternately, on failure it will return one of the following\n"
   "strings to indicate failure type.\n\n"
   "fault   -- There are no rows or columns, or the initial basis\n"
   "           is invalid, or the initial basis matrix is singular\n"
   "           or ill-conditioned.\n"
   "objll   -- The objective reached its lower limit.\n"
   "objul   -- The objective reached its upper limit.\n"
   "itlim   -- Iteration limited exceeded.\n"
   "tmlim   -- Time limit exceeded.\n"
   "sing    -- The basis matrix became singular or ill-conditioned.\n"
   "nopfs   -- No primal feasible solution. (Presolver only.)\n"
   "nodfs   -- No dual feasible solution. (Presolver only.)\n" },

  {"exact", (PyCFunction)LPX_solver_exact, METH_NOARGS,
   "exact()\n\n"
   "Attempt to solve the problem using an exact simplex method.\n\n"
   "This returns None if the problem was successfully solved.\n"
   "Alternately, on failure it will return one of the following\n"
   "strings to indicate failure type.\n\n"
   "fault   -- There are no rows or columns, or the initial basis\n"
   "           is invalid, or the initial basis matrix is singular\n"
   "           or ill-conditioned.\n"
   "itlim   -- Iteration limited exceeded.\n"
   "tmlim   -- Time limit exceeded." },

  {"interior", (PyCFunction)LPX_solver_interior, METH_NOARGS,
   "interior()\n\n"
   "Attempt to solve the problem using an interior-point method.\n\n"
   "This returns None if the problem was successfully solved.\n"
   "Alternately, on failure it will return one of the following\n"
   "strings to indicate failure type.\n\n"
   "fault   -- There are no rows or columns.\n"
   "nofeas  -- The problem has no feasible (primal/dual) solution.\n"
   "noconv  -- Very slow convergence or divergence.\n"
   "itlim   -- Iteration limited exceeded.\n"
   "instab  -- Numerical instability when solving Newtonian system." },

  {"integer", (PyCFunction)LPX_solver_integer, METH_NOARGS,
   "intopt()\n\n"
   "MIP solver based on branch-and-bound.\n\n"
   "This method requires a mixed-integer problem where an optimal\n"
   "solution to an LP relaxation (either through simplex() or\n"
   "exact()) has already been found.  Alternately, try intopt().\n\n"
   "This returns None if the problem was successfully solved.\n"
   "Alternately, on failure it will return one of the following\n"
   "strings to indicate failure type.\n\n"
   "fault   -- There are no rows or columns, or it is not a MIP\n"
   "           problem, or integer variables have non-int bounds.\n"
   "nopfs   -- No primal feasible solution.\n"
   "nodfs   -- Relaxation has no dual feasible solution.\n"
   "itlim   -- Iteration limited exceeded.\n"
   "tmlim   -- Time limit exceeded.\n"
   "sing    -- Error occurred solving an LP relaxation subproblem." },

  {"intopt", (PyCFunction)LPX_solver_intopt, METH_NOARGS,
   "intopt()\n\n"
   "More advanced MIP branch-and-bound solver than integer(). This\n"
   "variant does not require an existing LP relaxation.\n\n"
   "This returns None if the problem was successfully solved.\n"
   "Alternately, on failure it will return one of the following\n"
   "strings to indicate failure type.\n\n"
   "fault   -- There are no rows or columns, or it is not a MIP\n"
   "           problem, or integer variables have non-int bounds.\n"
   "nopfs   -- No primal feasible solution.\n"
   "nodfs   -- Relaxation has no dual feasible solution.\n"
   "itlim   -- Iteration limited exceeded.\n"
   "tmlim   -- Time limit exceeded.\n"
   "sing    -- Error occurred solving an LP relaxation subproblem." },

  {"kkt", (PyCFunction)LPX_kkt, METH_VARARGS,
   "kkt([scaled=False])\n\n"
   "Return an object encapsulating the results of a check on the\n"
   "Karush-Kuhn-Tucker optimality conditions for a basic (simplex)\n"
   "solution.  If the argument 'scaled' is true, return results \n"
   "of checking the internal scaled instance of the LP instead."},
  {"kktint", (PyCFunction)LPX_kktint, METH_NOARGS,
   "kktint()\n\n"
   "Similar to kkt(), except analyzes solution quality of an\n"
   "mixed-integer solution.  Note that only the primal components\n"
   "of the KKT object will have meaningful values."},

  // Data writing
  {"write", (PyCFunction)LPX_write, METH_VARARGS | METH_KEYWORDS,
   "write(format=filename)\n\n"
   "Output data about the linear program into a file with a given\n"
   "format.  What data is written, and how it is written, depends\n"
   "on which of the format keywords are used.  Note that one may\n"
   "specify multiple format and filename pairs to write multiple\n"
   "types and formats of data in one call to this function.\n\n"
   "mps       -- For problem data in the fixed MPS format.\n"
   "bas       -- The current LP basis in fixed MPS format.\n"
   "freemps   -- Problem data in the free MPS format.\n"
   "cpxlp     -- Problem data in the CPLEX LP format.\n"
   "glp       -- Problem data in the GNU LP format.\n"
   "prob      -- Problem data in a plain text format.\n"
   "sol       -- Basic solution in printable format.\n"
   "sens_bnds -- Bounds sensitivity information.\n"
   "ips       -- Interior-point solution in printable format.\n"
   "mip       -- MIP solution in printable format."},
  
  {NULL}
};

PyTypeObject LPXType = {
  PyObject_HEAD_INIT(NULL)
  0,					/* ob_size */
  "glpk.LPX",				/* tp_name */
  sizeof(LPXObject),			/* tp_basicsize*/
  0,					/* tp_itemsize*/
  (destructor)LPX_dealloc,		/* tp_dealloc*/
  0,					/* tp_print*/
  0,					/* tp_getattr*/
  0,					/* tp_setattr*/
  0,					/* tp_compare*/
  (reprfunc)LPX_Str,			/* tp_repr*/
  0,					/* tp_as_number*/
  0,					/* tp_as_sequence*/
  0,					/* tp_as_mapping*/
  0,					/* tp_hash */
  0,					/* tp_call*/
  (reprfunc)LPX_Str,			/* tp_str*/
  0,					/* tp_getattro*/
  0,					/* tp_setattro*/
  0,					/* tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,/* tp_flags*/
  "LPX()                 --> Empty linear program.\n"
  "LPX(gmp=filename)     --> Linear program with data read from a\n"
  "    GNU MathProg file containing model and data.\n"
  "LPX(mps=filename)     --> Linear program with data read from a\n"
  "    datafile in fixed MPS format.\n"
  "LPX(freemps=filename) --> Linear program with data read from a\n"
  "    datafile in free MPS format.\n"
  "LPX(cpxlp=filename)   --> Linear program with data read from a\n"
  "    datafile in fixed CPLEX LP format.\n"
  "LPX(glp=filename)     --> Linear program with data read from a\n"
  "    datafile in GNU LP format.\n"
  "LPX(gmp=(model_filename,[data_filename,[output_filename]])-->\n"
  "    Linear program from GNU MathProg input files.  The first\n"
  "    element is a path to the model second, the second to the\n"
  "    data section.  If the second element is omitted or is None\n"
  "    then the model file is presumed to also hold the data.\n"
  "    The third elment holds the output data file to write\n"
  "    display statements to.  If omitted or None, the output\n"
  "    is instead put through to standard output.\n"
  "\n"
  "This represents a linear program object.  It holds data and\n"
  "offers methods relevant to the whole of the linear program.\n"
  "There are many members in this class, but the most important\n"
  "are:\n"
  "  obj     Represents the objective function.\n"
  "  rows    A collection over which one can access rows.\n"
  "  cols    Same, but for columns.\n"
  "  params  Holds control parameters and statistics.",
	/* tp_doc */
  (traverseproc)LPX_traverse,		/* tp_traverse */
  (inquiry)LPX_clear,			/* tp_clear */
  0,					/* tp_richcompare */
  offsetof(LPXObject, weakreflist),	/* tp_weaklistoffset */
  0,					/* tp_iter */
  0,					/* tp_iternext */
  LPX_methods,				/* tp_methods */
  LPX_members,				/* tp_members */
  LPX_getset,				/* tp_getset */
  0,					/* tp_base */
  0,					/* tp_dict */
  0,					/* tp_descr_get */
  0,					/* tp_descr_set */
  0,					/* tp_dictoffset */
  (initproc)LPX_init,			/* tp_init */
  0,					/* tp_alloc */
  LPX_new,				/* tp_new */
};
