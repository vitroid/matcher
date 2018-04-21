//http://acooke.org/cute/ExampleCod0.html
//https://github.com/numpy/numpy/blob/master/numpy/core/src/dummymodule.c
//https://qiita.com/junkoda/items/17df11d7a20dc9d50e7d

#include "Python.h"
#include <numpy/arrayobject.h>
#include <math.h>
#include "matcher2.h"

static PyObject *matcher2(PyObject *self, PyObject* args);

static PyMethodDef module_methods[] = {
  {"matcher2",
   matcher2,
   METH_VARARGS,
   "Match atomic environments."},
  {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "cmatcher2",
  NULL,
  -1,
  module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};


//my initializer
PyMODINIT_FUNC PyInit_cmatcher(void) {
  PyObject *m;
  import_array();
  m = PyModule_Create(&moduledef);
  if (!m)
    return NULL;
  return m;
}



//taken from C_arraytest.c in Scipy.

/* ==== Check that PyArrayObject is a double (Float) type and a matrix ==============
    return 1 if an error and raise exception */ 
int  not_doublematrix(PyArrayObject *mat)  {
  if (mat->descr->type_num != NPY_DOUBLE || mat->nd != 2)  {
    PyErr_SetString(PyExc_ValueError,
                    "In not_doublematrix: array must be of type Float and 2 dimensional (n x m).");
    return 1;  }
  return 0;
}
//taken from C_arraytest.c in Scipy.

/* ==== Check that PyArrayObject is a double (Float) type and a vector ==============
    return 1 if an error and raise exception */ 
int  not_doublevector(PyArrayObject *mat)  {
  if (mat->descr->type_num != NPY_DOUBLE || mat->nd != 1)  {
    PyErr_SetString(PyExc_ValueError,
                    "In not_doublevector: array must be of type Float and 1 dimensional (n).");
    return 1;  }
  return 0;
}


static PyObject *matcher2(PyObject *self, PyObject* args)
{
  // arguments: Oatoms, cell, radius, rprox, every
  PyArrayObject *pos, *cell, *unitatoms, *unitcell;
  int adjdens;
  
  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "O!O!O!O!i",
			&PyArray_Type, &pos,
			&PyArray_Type, &cell,
			&PyArray_Type, &unitatoms,
			&PyArray_Type, &unitcell,
			&adjdens)) return NULL;
  if (NULL == pos)       return NULL;
  if (NULL == cell)      return NULL;
  if (NULL == unitatoms) return NULL;
  if (NULL == unitcell)  return NULL;
  if (not_doublematrix(pos))       return NULL;
  if (not_doublematrix(unitatoms)) return NULL;
  if (not_doublevector(cell))      return NULL;
  if (not_doublevector(unitcell))  return NULL;
  /* Get the dimensions of the input */
  int n=pos->dimensions[0];
  int nu=unitatoms->dimensions[0];

  double* a = (double*)pos->data;
  double* c = (double*)cell->data;
  double* au = (double*)unitatoms->data;
  double* cu = (double*)unitcell->data;
  int nostore=0;
  matchtype* matches = matcher2_core(n, a, c, nu, au, cu, adjdens, nostore);
  int nmatches = match_len(matches);
  //Py_DECREF(pos);
  //Py_DECREF(cell);
  //return value is a tuple of tuples.
  PyObject* result = PyTuple_New(nmatches);
  while ( matches != NULL ){
    matchtype* s = matches;
    //in reverse order
    nmatches -= 1;
    PyObject* list = PyTuple_New(s->natom);
    for(int i=0;i<s->natom;i++){
      PyTuple_SetItem(list,i,Py_BuildValue("i",s->atoms[i]));
    }
    PyObject* mat = PyTuple_New(9);
    for(int i=0;i<9;i++){
      PyTuple_SetItem(mat,i,Py_BuildValue("f",s->mat[i]));
    }
    PyTuple_SetItem(result,
		    nmatches,
		    Py_BuildValue("(fiiOO)",
				  s->rmsd,
				  s->atom_gro,
				  s->atom_unitcell,
				  mat,
				  list));
    matches = s->next;
    free(s);
  }
  return result;
}
