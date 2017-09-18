#include <fdct2_wrapper.hpp>

/*
This software was written by Darren Thomson and Gilles Hennenfent.
Copyright owned by The University of British Columbia, 2006.  This
software is distributed and certain rights granted under a License
Agreement with The University of British Columbia. Please contact the
authors or the UBC University Industry Liaison Office at 604-822-8580
for further information.
*/

PyObject* fdct2_wrapper(PyObject* input,int nbs,int nba,int ac)
{
	CpxNumMat f;
	vector<vector<CpxNumMat> > g;

	// cast the input to a numpy array
	PyArrayObject *in;
	in = (PyArrayObject*) input;

	// make a CpxNumMat with the data from in
	f._m = in->dimensions[1];
	f._n = in->dimensions[0];
	f._data = (cpx*) PyArray_DATA(in);

	// do the forward FDCT
	fdct_wrapping(f.m(),f.n(),nbs,nba,ac,f,g);

	// zero the size of the CpxNumMat and dereference its pointer
	f._m = 0;
	f._n = 0;
	f._data = NULL;

	// set up the list output
	PyListObject* out;
	PyArrayObject* output;
	out = (PyListObject*) PyList_New(0);
	npy_intp dims[2];
	int i,j;

	for(i=0;i<g.size();i++)
	{
		// append a list for this scale
		PyList_Append((PyObject*) out,PyList_New(0));

		for(j=0;j<g[i].size();j++)
		{
			// set the dimensions for this scale & angle
			dims[0] = g[i][j].n();
			dims[1] = g[i][j].m();

			// make an array for this scale & angle's data
			output = (PyArrayObject*) PyArray_SimpleNewFromData(2, dims, PyArray_CDOUBLE,g[i][j]._data);
			Py_INCREF((PyObject*) output);
		
			// append this angle's data to the list for this scale			
			PyList_Append(PyList_GetItem((PyObject*) out,i),(PyObject*) output);
									
			// zero the CpxNumMat for this scale & angle, hand ownership to numpy
			g[i][j]._data = NULL;
			g[i][j]._m = 0;
			g[i][j]._n = 0;
			output->flags = output->flags | NPY_OWNDATA;
		}
	}

	return (PyObject*) out;	
}


PyObject* ifdct2_wrapper(PyObject* in,int nbs,int nba,int ac,int n1,int n2)
{
	// Create objects for ifdct
	CpxNumMat f;
	vector<vector<CpxNumMat> > g;

	// Prepare to loop over wedges
	int i,j,angles;
	PyArrayObject* wedge;
	g.resize(nbs);
	
	for(i=0;i<PyList_Size(in);i++)
	{
		// Set the number of angles at this scale
		angles = PyList_Size(PyList_GetItem(in,i));
		g[i].resize(angles);

 		for(j=0;j<angles;j++)
		{
			// Get this wedge's data, put it in g 
			wedge = (PyArrayObject*) PyList_GetItem(PyList_GetItem(in,i),j);
			g[i][j]._data = (cpx*) wedge->data;			
			g[i][j]._m = wedge->dimensions[1];
			g[i][j]._n = wedge->dimensions[0];
		}
	}

	// Do ifdct on g, get f
	ifdct_wrapping(n1,n2,nbs,nba,ac,g,f);

	// Clean up g without deallocating
	for(i=0;i<PyList_Size(in);i++)
	{
 		for(j=0;j<PyList_Size(PyList_GetItem(in,i));j++)
		{
			g[i][j]._data = NULL;
			g[i][j]._m = 0;
			g[i][j]._n = 0;
		}		
	}

	// Prepare to create output array	
	npy_intp dims[2];
	dims[0] = n2;
	dims[1] = n1;	

	// Make a numpy array from f
	PyArrayObject* output;
	output = (PyArrayObject*)  PyArray_SimpleNewFromData(2, dims, PyArray_CDOUBLE,f._data);
	output->flags = output->flags | NPY_OWNDATA;
	Py_INCREF((PyObject*) output);

	// Clean f without deallocating
	f._data = NULL;
	f._m = 0;
	f._n = 0;

	return PyArray_Transpose(output,NULL);
}



PyObject* fdct2_param_wrapper(int n1,int n2,int nbs,int nba,int ac)
{
	PyObject *out;
	PyObject *psx,*psy,*pfx,*pfy,*pnx,*pny;
	
	vector< vector<double> > sx,sy,fx,fy;
	vector< vector<int> > nx,ny;
	
	fdct_wrapping_param(n1,n2,nbs,nba,ac,sx,sy,fx,fy,nx,ny);
	
	int i,j;

	psx = PyList_New(sx.size());
	for(i=0;i<sx.size();i++)
	{
		PyList_SetItem(psx,i,PyList_New(sx[i].size()));
		for(j=0;j<sx[i].size();j++)
		{
			PyList_SetItem(PyList_GetItem(psx,i),j,PyFloat_FromDouble(sx[i][j]));
			Py_INCREF(PyList_GetItem(PyList_GetItem(psx,i),j));
		}
	}
	psy = PyList_New(sy.size());
	for(i=0;i<sy.size();i++)
	{
		PyList_SetItem(psy,i,PyList_New(sy[i].size()));
		for(j=0;j<sy[i].size();j++)
		{
			PyList_SetItem(PyList_GetItem(psy,i),j,PyFloat_FromDouble(sy[i][j]));
			Py_INCREF(PyList_GetItem(PyList_GetItem(psy,i),j));
		}
	}
	pfx = PyList_New(fx.size());
	for(i=0;i<fx.size();i++)
	{
		PyList_SetItem(pfx,i,PyList_New(fx[i].size()));
		for(j=0;j<fx[i].size();j++)
		{
			PyList_SetItem(PyList_GetItem(pfx,i),j,PyFloat_FromDouble(fx[i][j]));
			Py_INCREF(PyList_GetItem(PyList_GetItem(pfx,i),j));
		}
	}
	pfy = PyList_New(fy.size());
	for(i=0;i<fy.size();i++)
	{
		PyList_SetItem(pfy,i,PyList_New(fy[i].size()));
		for(j=0;j<fy[i].size();j++)
		{
			PyList_SetItem(PyList_GetItem(pfy,i),j,PyFloat_FromDouble(fy[i][j]));
			Py_INCREF(PyList_GetItem(PyList_GetItem(pfy,i),j));
		}
	}
	pnx = PyList_New(nx.size());
	for(i=0;i<nx.size();i++)
	{
		PyList_SetItem(pnx,i,PyList_New(nx[i].size()));
		for(j=0;j<nx[i].size();j++)
		{
			PyList_SetItem(PyList_GetItem(pnx,i),j,PyInt_FromLong(nx[i][j]));
			Py_INCREF(PyList_GetItem(PyList_GetItem(pnx,i),j));
		}
	}
	pny = PyList_New(ny.size());
	for(i=0;i<ny.size();i++)
	{
		PyList_SetItem(pny,i,PyList_New(ny[i].size()));
		for(j=0;j<ny[i].size();j++)
		{
			PyList_SetItem(PyList_GetItem(pny,i),j,PyInt_FromLong(ny[i][j]));
			Py_INCREF(PyList_GetItem(PyList_GetItem(pny,i),j));
		}
	}
		

	out = PyTuple_New(6);
	PyTuple_SetItem(out,1,psx);
	PyTuple_SetItem(out,0,psy);
	PyTuple_SetItem(out,3,pfx);
	PyTuple_SetItem(out,2,pfy);
	PyTuple_SetItem(out,5,pnx);
	PyTuple_SetItem(out,4,pny);
	
	return out;
}


void fdct_init()
{
	// Load numpy API at the start
	import_array();
}


