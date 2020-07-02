#include <fdct_wrapper.hpp>

#include <fdct3d.hpp>

/*
This software was written by Darren Thomson and Gilles Hennenfent.
Copyright owned by The University of British Columbia, 2006.  This
software is distributed and certain rights granted under a License
Agreement with The University of British Columbia. Please contact the
authors or the UBC University Industry Liaison Office at 604-822-8580
for further information.
*/

PyObject* fdct3_wrapper(PyObject* input,int nbs,int nbz,int ac);

PyObject* ifdct3_wrapper(PyObject* in,int nbs,int nbz,int ac,int n1,int n2,int n3);

PyObject* fdct3_param_wrapper(int n1,int n2,int n3,int nbs,int nbz,int ac);






