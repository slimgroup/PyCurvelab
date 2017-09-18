%module fdct2_wrapper

%{
#include <fdct2_wrapper.hpp>
%}


extern PyObject* fdct2_wrapper(PyObject* input,int nbs,int nba,int ac);

extern PyObject* ifdct2_wrapper(PyObject* in,int nbs,int nba,int ac,int n1,int n2);

extern PyObject* fdct2_param_wrapper(int n1,int n2,int nbs,int nba,int ac);

extern void fdct_init();

