# This software was written by Darren Thomson and Gilles Hennenfent.
# Copyright owned by The University of British Columbia, 2006.  This
# software is distributed and certain rights granted under a License
# Agreement with The University of British Columbia. Please contact
# the authors or the UBC University Industry Liaison Office at
# 604-822-8580 for further information.

import pyct as ct
import numpy as np
import scipy as sc

norm = np.linalg.norm


def test(dim=2,clen=10):

	for i in xrange(clen):
		print "-----------------------------------"
		if dim==2:
			sz = np.arange(256,513)
		elif dim==3:
			sz = np.arange(64,128)

		np.random.shuffle(sz)
		iscplx = [True,False]
		np.random.shuffle(iscplx)
		if iscplx[0]:
			print "Complex input"
			f = np.array(sc.randn(*sz[:dim])+sc.randn(*sz[:dim])*1j)
		else:
			print "Real input"
			f = np.array(sc.randn(*sz[:dim]))
			
		isac = [True,False]
		np.random.shuffle(isac)
		if isac[0]:
			print "All curvelets"
		else:
			print "Wavelets at finest scale"


		print f.shape
		
		if dim==2:
			A = ct.fdct2(f.shape,6,32,isac[0],cpx=iscplx[0])
		elif dim==3:
			A = ct.fdct3(f.shape,4,8,isac[0],cpx=iscplx[0])

		x = A.fwd(f)
		
		if np.allclose(norm(f.flatten(),ord=2),norm(x,ord=2) ):
			print 'Energy check ok!'
		else:
			print 'Problem w energy test'

		fr = A.inv(x)
		if np.allclose(f.flatten(),fr.flatten() ):
			print 'Inverse check ok!'
		else:
			print 'Problem w inverse test'

		print "||f|| = ",norm(f.flatten(),ord=2),f.dtype
		print "||x|| = ",norm(x,ord=2),x.dtype
		print "||fr|| = ",norm(fr.flatten(),ord=2),fr.dtype


def normtest(dim=2,clen=10):
	for i in xrange(clen):
		print "-----------------------------------"
		if dim==2:
			sz = np.arange(256,513)
		else:
			sz = np.arange(64,128)

		np.random.shuffle(sz)
		print sz[:dim]

		iscplx = [True,False]
		np.random.shuffle(iscplx)

		isac = [True,False]
		np.random.shuffle(isac)
		if isac[0]:
			print "All curvelets"
		else:
			print "Wavelets at finest scale"

		if dim==2:		
			A = ct.fdct2(sz[:2],6,32,isac[0],norm=True,cpx=iscplx[0])
			pos = np.arange(A.range())
			np.random.shuffle(pos)
			if iscplx[0]:
				print "Complex input"
				x = np.zeros(A.range(),dtype='complex')
				v = np.random.rand()
				x[pos[0]] = v +np.sqrt(1-v**2)*1j
			else:
				print "Real input"
				x = np.zeros(A.range())
				x[pos[0]] = 1.
		elif dim==3:
			A = ct.fdct3(sz[:3],4,8,isac[0],norm=True,cpx=iscplx[0])
			pos = np.arange(A.range())
			np.random.shuffle(pos)
			if iscplx[0]:
				print "Complex input"
				x = np.zeros(A.range(),dtype='complex')
				v = np.random.rand()
				x[pos[0]] = v +np.sqrt(1-v**2)*1j
			else:
				print "Real input"
				x = np.zeros(A.range())
				x[pos[0]] = 1.

		f = A.inv(x)

		if np.allclose(norm(f.flatten(),ord=2),norm(x,ord=2) ):
			print 'Norm check ok!'
		else:
			print 'Problem w norm test'

		print "||f|| = ",norm(f.flatten(),ord=2),f.dtype
		print "||x|| = ",norm(x,ord=2),x.dtype
