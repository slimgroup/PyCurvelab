#include <fdct_wrapper.hpp>

#include <fdct_wrapping.hpp>
using namespace fdct_wrapping_ns;

/*
This software was written by Darren Thomson and Gilles Hennenfent.
Copyright owned by The University of British Columbia, 2006.  This
software is distributed and certain rights granted under a License
Agreement with The University of British Columbia. Please contact the
authors or the UBC University Industry Liaison Office at 604-822-8580
for further information.
*/

PyObject* fdct2_wrapper(PyObject* input,int nbs,int nba,int ac);

PyObject* ifdct2_wrapper(PyObject* in,int nbs,int nba,int ac,int n1,int n2);

PyObject* fdct2_param_wrapper(int n1,int n2,int nbs,int nba,int ac);






