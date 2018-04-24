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
PyMODINIT_FUNC PyInit_cmatcher2(void) {
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
  PyArrayObject *gatoms, *gcell, *uatoms, *ucell;
  int adjdens;
  int nostore;
  
  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "O!O!O!O!ii",
			&PyArray_Type, &gatoms,
			&PyArray_Type, &gcell,
			&PyArray_Type, &uatoms,
			&PyArray_Type, &ucell,
			&adjdens,
                        &nostore)) return NULL;
  if (NULL == gatoms) return NULL;
  if (NULL == gcell)  return NULL;
  if (NULL == uatoms) return NULL;
  if (NULL == ucell)  return NULL;
  if (not_doublematrix(gatoms)) return NULL;
  if (not_doublematrix(uatoms)) return NULL;
  if (not_doublevector(gcell))  return NULL;
  if (not_doublevector(ucell))  return NULL;
  /* Get the dimensions of the input */
  match* matches = matcher2_core(gatoms->dimensions[0],
                                 (double*)gatoms->data,
                                 (double*)gcell->data,
                                 uatoms->dimensions[0],
                                 (double*)uatoms->data,
                                 (double*)ucell->data,
                                 adjdens,
                                 nostore);
  int nmatches = match_len(matches);
  //Py_DECREF(gatoms);
  //Py_DECREF(cell);
  //return value is a tuple of tuples.
  PyObject* result = PyTuple_New(nmatches);
  while ( matches != NULL ){
    match* s = matches;
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
				  s->gcenter,
				  s->ucenter,
				  mat,
				  list));
    matches = s->next;
    free(s);
  }
  return result;
}
