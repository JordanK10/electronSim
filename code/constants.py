import numpy as np
import itertools
import matplotlib.pyplot as plt
import math
import matplotlib.animation as animation
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import HTML
import matplotlib.gridspec as gridspec
from random import *
import sys
import time
import os

#Space dimension Parameters
LEN = 20
WID = 100
DEP = 5
A = LEN*WID
tau = 3
DT = .1
f = 1

RANK = 2
if RANK == 2:
    TOT = LEN*WID
    DIM = [LEN,WID]
    E= np.array([f,0])
else:
    TOT = LEN*WID*DEP
    DIM = [LEN,WID,DEP]
    E=np.array([f,0,0])
#
# #Particle Parameters
# PROBSCAT=True
VEL = 2
SIZE = int(100)
Me = 1 #MeV/c (mass of electron)
Ce= 1
Se = .5 #

Cp = 1
Cb = 10


# #Scatterer Parameters
BULKSCAT = True
SSCAT = True
DETSCAT = False
if DETSCAT:
    SSCAT = False
P=1
SCATTERDENS = 1 #How many scatterers per unit area
SRAD=1
DSRAND = True


THERM = False
THERMGRAD = False
THERMRES = 3
TAUtempDEP = True
FIXEDEDGES = False

SMIL = True

T0=1
TGRAD=50

LGD = True


#FOR PERIODIC SCATTERERS
#xct controls number of scatterers in  row, and yct controls spacing between rows
#Scatdepth controls the sharpness of the scatterers
    #Scatterers are only symmetric if the scattering counts are symmetric
    #Implementing a duty factor?
#Diag controls the diagonality of the wave
SCATPAT = "periodic"

#Sim Parameters
FPS = 30
FRAMES =500
#
# #Block Parameters
# LEG=1
# ANG=np.pi/2
# RAD = 1
#
# V="V"
# S="S"
# #y = cint(y)%DIV+((y%(2*DIV))/DIV)*(DIV-2*(y%DIV))
