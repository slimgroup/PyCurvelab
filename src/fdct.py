# This software was written by Darren Thomson and Gilles Hennenfent.
# Copyright owned by The University of British Columbia, 2006.  This
# software is distributed and certain rights granted under a License
# Agreement with The University of British Columbia. Please contact
# the authors or the UBC University Industry Liaison Office at
# 604-822-8580 for further information.

import numpy as _fdct__n
import math as _fdct__math
from scipy import fftpack as _fdct__fft

from CLarray import *

class fdct(object):
	"""
	FDCT Object 
	
	Should not be called directly
	
	Call subclasses fdct2 and fdct3
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
		self.n = n
		self.nbs = nbs
		self.nba = nba
		self.ac = ac
		self.norm = norm
		self.vec = vec
		self.cpx = cpx
		self.E = None

	def domain(self):
		try:
			return self.domain_
		except:
			self.domain_ = __n.prod(self.n)
			return self.domain_

	def range(self):
		try:
			return self.range_
		except:
			self.range_ = 0
			for sc in self.sizes:
				for an in sc:
					self.range_ += __n.prod(an)			
			return self.range_

	def vect(self,c):
		"""
		create a curvelet vector from a 
		curvelet structure c, return
		the vector
		"""
		vlist = []
		for i in c:
			for j in i:
				vlist.append(j.ravel())
		return CLarray(__n.concatenate(vlist),fdct_op=self)

	def struct(self,v):
		"""
		create a curvelet structure from a 
		curvelet vector v, return
		the structure
		"""
		ind = 0
		c = []
		for i in self.sizes:
			c.append([])
			for j in i:
				if not self.cpx:
					c[len(c)-1].append(__n.array(v[ind:(ind+__n.prod(j))].reshape(j)))
				else:	
					c[len(c)-1].append(__n.array(__n.complex_(v[ind:(ind+__n.prod(j))].reshape(j))))
				ind += __n.prod(j)
		return c 


	def c2r(self,c):
		"""
		convert coefficients in c from complex
		to real curvelet coefficients
		
		this is an in-place operation
		"""
		c[0][0] = c[0][0].real

		if self.ac:
			end = len(c)
		else:
			end = len(c) - 1
			c[self.nbs-1][0] = c[self.nbs-1][0].real
	
		for i in range(1,end):
			for j in range(len(c[i])/2):
				c[i][j] = __n.sqrt(2) * c[i][j].real
				c[i][j+len(c[i])/2] = __n.sqrt(2) * c[i][j+len(c[i])/2].imag
		c[i][0] = c[i][0].real		

	def r2c(self,c):		
		"""
		convert coefficients in c from real
		to complex curvelet coefficients
		
		this is an in-place operation
		"""
		c[0][0] = __n.complex_(c[0][0])
		
		if self.ac:
			end = len(c)
		else:
			end = len(c) - 1
			c[self.nbs-1][0] = __n.complex_(c[self.nbs-1][0])

		for i in range(1,len(c)):
			for j in range(len(c[i])/2):
				c[i][j] = __n.complex_(c[i][j])
				c[i][j].imag = -c[i][j+len(c[i])/2]
				c[i][j] /= __n.sqrt(2)
				c[i][j+len(c[i])/2] = c[i][j].conj()
		
	def getindex(self,loc):
		"""
		returns the index in the curvelet 
		vector for a location given as a tuple
		or list of length 4, i.e. 
		(scale,angle,row,column)
		"""
		index = 0
		for i in range(loc[0]):
			for j in range(len(self.sizes[i])):
				index += self.sizes[i][j][0] * self.sizes[i][j][1]
		
		for i in range(loc[1]):
			index += self.sizes[loc[0]][i][0] * self.sizes[loc[0]][i][1]
		
		index += self.sizes[loc[0]][loc[1]][1] * loc[2] + loc[3]				
		return index

	def index(self,*loc):	
		"""
		returns the index in the curvelet 
		vector for a location
		can take 0-4 arguments
		for 1-3 arguments, returns the starting
		and ending points for that scale or
		scale/angle or scale/angle/row
		"""
		loc = list(loc)
		loclen = len(loc)
		while len(loc)<4:
			loc.append(0)
		
		index = self.getindex(loc)
		
		if loclen==2 and loc[1]==len(self.sizes[loc[0]])-1:
			if loc[0]==len(self.sizes)-1:
				outdex = self.range()
			else:
				loc[1] = 0
				loc[0] += 1
				outdex = self.getindex(loc)
			return (index,outdex)
		if loclen==1 and loc[0]==len(self.sizes)-1:
			return (index,self.range())
		elif 0<loclen<4:
			loc[loclen-1] += 1
			outdex = self.getindex(loc)
			return (index,outdex)
		else:
			return index

	def loc(self,index):
		"""
		returns the (scale,angle,row,column) 
		for a given index in a curvelet vector
		
		CURRENTLY NOT IMPLEMENTED
		"""
		pass
		
	def normstruct(self):
		"""
		returns a structure where E[sc][an] is the energy
		of curvelets at that scale and angle
		"""
		X = __n.ones(self.n)
		x = __fft.fftshift(__fft.ifft2(X)) * __n.sqrt(__n.prod(X.shape));

		myvec = self.vec
		
		sn = self.norm
		self.norm = False
		sv = self.vec
		self.vec = False
		sc = self.cpx
		self.cpx = True
		
		c = self.fwd(x)
		
		self.norm = sn
		self.vec = sv
		self.cpx = sc

		E = []
		for i in range(len(self.sizes)):
			E.append([])
			for j in range(len(self.sizes[i])):
				val = __n.sqrt(sum(sum(c[i][j]*c[i][j].conj())) / __n.prod(self.sizes[i][j]))
				E[len(E)-1].append(val)
		
		self.vec = myvec
		return E
				
	def normalize(self,c,inplace=True):
		"""
		normalizes the curvelet coefficients in c
			-> c must be complex!!!
		
		this is an in-place operation, unless forced
		otherwise, in which case it returns a new 
		structure
		"""
		if self.E is None:
			self.E = self.normstruct()

		if not inplace:
			cout = []

		for i in range(len(self.sizes)):
			if not inplace:
				cout.append([])
			for j in range(len(self.sizes[i])):
				if inplace:
					c[i][j] /= self.E[i][j]		
				else:
					cout[i].append(c[i][j].copy())
					cout[i][j] /= self.E[i][j]
				
		if not inplace:
			return cout
				
	def normvec(self):
		"""
		return a curvelet vector with the energy of each coefficient
		"""
		if self.E is None:
			self.E = self.normstruct()
				
		nv = __n.zeros(self.range())
		indi = 0
		indf = 0
		for i in range(len(self.sizes)):
			for j in range(len(self.sizes[i])):
				indf += __n.prod(self.sizes[i][j])
				nv[indi:indf] = self.E[i][j]
				indi = indf
		
		return CLarray(nv,fdct_op=self)
		
		
