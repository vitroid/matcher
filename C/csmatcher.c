//http://acooke.org/cute/ExampleCod0.html
//https://github.com/numpy/numpy/blob/master/numpy/core/src/dummymodule.c
//https://qiita.com/junkoda/items/17df11d7a20dc9d50e7d




#include "Python.h"
#include <numpy/arrayobject.h>
#include <math.h>
#include "smatcher.h"

static PyObject *smatcher(PyObject *self, PyObject* args);

static PyMethodDef module_methods[] = {
  {"smatcher", smatcher, METH_VARARGS,
   "Match atomic environments."},
  {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "csmatcher",
        NULL,
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL
};


//my initializer
PyMODINIT_FUNC PyInit_csmatcher(void) {
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


static PyObject *smatcher(PyObject *self, PyObject* args) {
  // arguments: Oatoms, cell, radius, rmsdmax, every
  PyArrayObject *pos, *cell;
  float radius, rmsdmax;
  int every;

  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "O!O!ffi", &PyArray_Type,
                        &pos, &PyArray_Type, &cell, &radius, &rmsdmax, &every)) return NULL;
  if (NULL == pos) return NULL;
  if (not_doublematrix(pos)) return NULL;
  if (not_doublevector(cell)) return NULL;
  /* Get the dimensions of the input */
  int n=pos->dimensions[0];

  double* a = (double*)pos->data;
  double* c = (double*)cell->data;
  smatchtype* smatch = smatcher_core(n, a, c, radius, rmsdmax, every);
  int nsmatch = smatch_len(smatch);
  //Py_DECREF(pos);
  //Py_DECREF(cell);
  //return value is a tuple of tuples.
  PyObject* result = PyTuple_New(nsmatch);
  while ( smatch != NULL ){
    smatchtype* s = smatch;
    //in reverse order
    nsmatch -= 1;
    PyTuple_SetItem(result,
                    nsmatch,
                    Py_BuildValue("(iif(fff)f)",
                                  s->i,
                                  s->j,
                                  s->radius,
                                  s->d[0], s->d[1], s->d[2],
                                  s->rmsd)
                    );
    smatch = s->next;
    free(s);
  }
  return result;
}
