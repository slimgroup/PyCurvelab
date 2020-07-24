%module fdct3_wrapper

%{
#include <fdct3_wrapper.hpp>
%}


extern PyObject* fdct3_wrapper(PyObject* input,int nbs,int nbz,int ac);

extern PyObject* ifdct3_wrapper(PyObject* in,int nbs,int nbz,int ac,int n1,int n2,int n3);

extern PyObject* fdct3_param_wrapper(int n1,int n2,int n3,int nbs,int nbz,int ac);

extern void fdct_init();


