"""
Fast Discrete Curvelet Transform

Uses libraries from CurveLab 2.0 

developed by:
Emmanuel Candes, Laurent Demanet, and Lexing Ying
California Institute of Technology

this wrapper developed by:
Darren Thomson and Gilles Hennenfent
Seismic Laboratory for Imaging and Modeling
University of British Columbia
"""

# This software was written by Darren Thomson and Gilles Hennenfent.
# Copyright owned by The University of British Columbia, 2006.  This
# software is distributed and certain rights granted under a License
# Agreement with The University of British Columbia. Please contact
# the authors or the UBC University Industry Liaison Office at
# 604-822-8580 for further information.

try:
	from fdct2 import *
except:
	pass

try:
	from fdct3 import *
except:
	pass
		
from test import test,normtest
					

