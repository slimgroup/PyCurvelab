#!/usr/bin/env python

# This software was written by Darren Thomson and Gilles Hennenfent.
# Copyright owned by The University of British Columbia, 2006.  This
# software is distributed and certain rights granted under a License
# Agreement with The University of British Columbia. Please contact
# the authors or the UBC University Industry Liaison Office at
# 604-822-8580 for further information.

import sys
import os
import string

FFTW = os.environ['FFTW']
FDCT = os.environ['FDCT']

for arg in sys.argv:
	a = string.split(arg,'=')
	if len(a)==2:
		if a[0] in ['FFTW','FDCT']:
			exec '%s = "%s"' %tuple(a)

fftw_inc = os.path.join(FFTW,'include')
fftw_lib = os.path.join(FFTW,'lib')
fdct2 = os.path.join(FDCT,'fdct_wrapping_cpp','src')
fdct3 = os.path.join(FDCT,'fdct3d','src')
pycl_inc = 'src'

try:
	import numpy
	npy_inc = os.path.join(os.path.split(numpy.__file__)[0],'core','include','numpy')
except:
	print "ERROR: numpy installation is necessary to install CurveLab"


from distutils.core import setup,Extension


setup(name='pyct',
		version='1.0',
		description='Python Wrappers for CurveLab-2.0',
		author='Darren Thomson and Gilles Hennenfent',
		author_email='dthomson,ghennenfent@eos.ubc.ca',
		url='http://slim.eos.ubc.ca',
		ext_package='pyct',
		ext_modules=[Extension('_fdct2_wrapper',[os.path.join('src','fdct2_wrapper.cpp'),os.path.join('src','fdct2_wrapper.i')],
						include_dirs=[fdct2,fftw_inc,npy_inc,pycl_inc],
						library_dirs=[fdct2,fftw_lib],libraries=['fdct_wrapping','fftw']),
					Extension('_fdct3_wrapper',[os.path.join('src','fdct3_wrapper.cpp'),os.path.join('src','fdct3_wrapper.i')],
						include_dirs=[fdct3,fftw_inc,npy_inc,pycl_inc],
						library_dirs=[fdct3,fftw_lib],libraries=['fdct3d','fftw'])],
		package_dir = {'pyct':'src'},
		packages=['pyct']
	 )			
