#include <fdct3_wrapper.hpp>

/*
This software was written by Darren Thomson and Gilles Hennenfent.
Copyright owned by The University of British Columbia, 2006.  This
software is distributed and certain rights granted under a License
Agreement with The University of British Columbia. Please contact the
authors or the UBC University Industry Liaison Office at 604-822-8580
for further information.
*/

PyObject* fdct3_wrapper(PyObject* input,int nbs,int nbz,int ac)
{
	CpxNumTns f;
	vector<vector<CpxNumTns> > g;

	// cast the input to a numpy array
	PyArrayObject* in;
	in = (PyArrayObject*) input;

	// make a CpxNumTns with the data from in
	f._m = in->dimensions[2];
	f._n = in->dimensions[1];
	f._p = in->dimensions[0];
	f._data = (cpx*) PyArray_DATA(in);

	// do the forward FDCT
	fdct3d_forward(f.m(),f.n(),f.p(),nbs,nbz,ac,f,g);

	// zero the size of the CpxNumTns and dereference its pointer
	f._m = 0;
	f._n = 0;
	f._p = 0;
	f._data = NULL;

	// set up the list output
	PyListObject* out;
	PyArrayObject* output;
	out = (PyListObject*) PyList_New(0);
	npy_intp dims[3];
	int i,j;

	for(i=0;i<g.size();i++)
	{
		// Append a list for this scale.  The inner list
		// starts with a refcount of 1, and then appending it
		// to the outer list immediately increments the
		// refcount to 2.
		PyList_Append((PyObject*) out,PyList_New(0));

		for(j=0;j<g[i].size();j++)
		{
			// set the dimensions for this scale & angle
			dims[0] = g[i][j].p();
			dims[1] = g[i][j].n();
			dims[2] = g[i][j].m();

			// make an array for this scale & angle's
			// data.  Starts with initial refcount of 1
                        output = (PyArrayObject*) PyArray_SimpleNewFromData(3, dims, PyArray_CDOUBLE,g[i][j]._data);
			// Append this angle's data to the list for
                        // this scale.  This increases the reference
                        // count of 'output' to 2
                        PyList_Append(PyList_GetItem((PyObject*) out,i),(PyObject*) output);

                        // Decrement the reference count of 'output'
                        // back to 1.  Now the above list holds the
                        // only reference count
                        Py_DECREF((PyObject*) output);
			
			// zero the CpxNumMat for this scale & angle, hand ownership to numpy
			g[i][j]._data = NULL;
			g[i][j]._m = 0;
			g[i][j]._n = 0;
			g[i][j]._p = 0;
			output->flags = output->flags | NPY_OWNDATA;
		}

                // Decrement the refcount of the inner list back to 1.
                // Now the outer list holds the only refcount.
                Py_DECREF((PyObject*) PyList_GetItem((PyObject*) out,i));

	}

        // In the end, the only remaining refcount we hold is to the
        // outer list.
	return (PyObject*) out;	
}


PyObject* ifdct3_wrapper(PyObject* in,int nbs,int nbz,int ac,int n1,int n2,int n3)
{
	// Create objects for ifdct
	CpxNumTns f;
	vector<vector<CpxNumTns> > g;

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
			g[i][j]._m = wedge->dimensions[2];
			g[i][j]._n = wedge->dimensions[1];
			g[i][j]._p = wedge->dimensions[0];
		}
	}
	
	// Do ifdct on g, get f
	fdct3d_inverse(n1,n2,n3,nbs,nbz,ac,g,f);

	// Clean up g without deallocating
	for(i=0;i<PyList_Size(in);i++)
	{
 		for(j=0;j<PyList_Size(PyList_GetItem(in,i));j++)
		{
			g[i][j]._data = NULL;
			g[i][j]._m = 0;
			g[i][j]._n = 0;
			g[i][j]._p = 0;
		}		
	}

	// Prepare to create output array	
	npy_intp dims[3];
	dims[0] = n3;
	dims[1] = n2;
	dims[2] = n1;	

	// Make a numpy array from f
	PyArrayObject* output;
	output = (PyArrayObject*)  PyArray_SimpleNewFromData(3, dims, PyArray_CDOUBLE,f._data);
	output->flags = output->flags | NPY_OWNDATA;
	Py_INCREF((PyObject*) output);

	// Clean f without deallocating
	f._data = NULL;
	f._m = 0;
	f._n = 0;
	f._p = 0;

	return PyArray_Transpose(output,NULL);
}



PyObject* fdct3_param_wrapper(int n1,int n2,int n3,int nbs,int nbz,int ac)
{
	PyObject *out;
	PyObject *pfx,*pfy,*pfz,*pnx,*pny,*pnz;
	
	vector< vector<double> > fx,fy,fz;
	vector< vector<int> > nx,ny,nz;
	
	fdct3d_param(n1,n2,n3,nbs,nbz,ac,fx,fy,fz,nx,ny,nz);
	
	int i,j;

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
	pfz = PyList_New(fz.size());
	for(i=0;i<fz.size();i++)
	{
		PyList_SetItem(pfz,i,PyList_New(fz[i].size()));
		for(j=0;j<fz[i].size();j++)
		{
			PyList_SetItem(PyList_GetItem(pfz,i),j,PyFloat_FromDouble(fz[i][j]));
			Py_INCREF(PyList_GetItem(PyList_GetItem(pfz,i),j));
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
	pnz = PyList_New(nz.size());
	for(i=0;i<nz.size();i++)
	{
		PyList_SetItem(pnz,i,PyList_New(nz[i].size()));
		for(j=0;j<nz[i].size();j++)
		{
			PyList_SetItem(PyList_GetItem(pnz,i),j,PyInt_FromLong(nz[i][j]));
			Py_INCREF(PyList_GetItem(PyList_GetItem(pnz,i),j));
		}
	}
		

	out = PyTuple_New(6);
	PyTuple_SetItem(out,2,pfx);
	PyTuple_SetItem(out,1,pfy);
	PyTuple_SetItem(out,0,pfz);
	PyTuple_SetItem(out,5,pnx);
	PyTuple_SetItem(out,4,pny);
	PyTuple_SetItem(out,3,pnz);
	
	return out;
}


void fdct_init()
{
	// Load numpy API at the start
	import_array();
}


