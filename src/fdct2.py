# This software was written by Darren Thomson and Gilles Hennenfent.
# Copyright owned by The University of British Columbia, 2006.  This
# software is distributed and certain rights granted under a License
# Agreement with The University of British Columbia. Please contact
# the authors or the UBC University Industry Liaison Office at
# 604-822-8580 for further information.

import sys as __sys
import fdct2_wrapper as _fdct2__f2

import numpy as _fdct2__n
import fdct 

_fdct2__f2.fdct_init()

class fdct2(fdct.fdct):
	"""
	2D FDCT using wrapping

	from fdct_wrapping_cpp in CurveLab 2.0 

	Object to define forward and inverse transform
	operations, as well as operations on data in 
	the transformed domain
	"""
	def __init__(self,n,nbs,nba,ac,norm=False,vec=True,cpx=False):
		""" 
		n - a tuple with the shape of the 2D input array
		nbs - number of scales
		nba - number of angles at 2nd coarsest scale
		ac - use curvelets at coarsest scale (true/false)
		norm - normalize columns of the transform (true/FALSE)
		vec - return curvelet coefficients in a vector (TRUE/false)
				-> vector type is a specially defined "CLarray"
				-> CLarray has contiguous data, can be referenced as a 
					long array or by scales & angles 
				-> see doc for CLarray and CLsarray for more info
				-> false will return curvelet coefficients in a 
					list(list(array)) with non-contiguous data
                cpx - assume complex inputs
                                -> does not convert complex coefficients to
                                        real ones after the forward transform
		"""
		fdct.fdct.__init__(self,n,nbs,nba,ac,norm,vec,cpx)
		self.s,self.f,self.sizes = self.param()
	
	def param(self):
		"""
		Get transform parameters
		
		returns (s,f,n) where:
		s - 
		f - 
		n - 
		"""
		parm = __f2.fdct2_param_wrapper(self.n[0],self.n[1],self.nbs,self.nba,self.ac)
		
		s = []
		f = []
		n = []
		
		for i in range(len(parm[0])):
			s.append([])
			f.append([])
			n.append([])
			
			for j in range(len(parm[0][i])):
				s[i].append((parm[0][i][j],parm[1][i][j]))
				f[i].append((parm[2][i][j],parm[3][i][j]))
				n[i].append((parm[4][i][j],parm[5][i][j]))
		
		return s,f,n
	
	def fwd(self,a):
		"""
		apply forward transform to a,
		return coefficients
		"""
		a = __n.complex_(a)	
		
		if a.flags.carray:
			a = a.transpose().copy()

		c = __f2.fdct2_wrapper(a,self.nbs,self.nba,self.ac)
		
		if self.norm:
			self.normalize(c)
		
		if not self.cpx: 
			self.c2r(c)

		if self.vec:
			return self.vect(c)
		else:
			return c				

	def inv(self,c):
		"""
		apply inverse transform to coefficients c,
		return data
		"""
		if self.vec:
			c = self.struct(c)
		
		if self.norm:
			c = self.normalize(c,inplace=False)
		
		if not self.cpx:
			self.r2c(c)

		a =  __f2.ifdct2_wrapper(c,self.nbs,self.nba,self.ac,self.n[0],self.n[1])

		if not self.cpx: 
			a = a.real
			self.c2r(c)

		return a
		

