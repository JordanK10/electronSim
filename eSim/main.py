#########
##INCLUDES
#########

import numpy as np
import itertools
import matplotlib.pyplot as plt
import math
import matplotlib.animation as animation
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import HTML
import matplotlib.gridspec as gridspec
from random import *
import sys
import time
import os
import json
from bokeh.plotting import figure, curdoc
from bokeh.driving import linear
from bokeh.layouts import gridplot,row,column
from bokeh.models import ColumnDataSource,Slider,CustomJS, Range1d
from bokeh.io import output_file, show
from bokeh.embed import json_item
from bokeh.resources import CDN
from threading import Thread

from flask import Flask,render_template
from jinja2 import Template
import random

#Space dimension Parameters
LEN = 20
WID = 10
DEP = 5
A = LEN*WID
tau = [2]
DT = .01
f = 3

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
SIZE = int(250)
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
# RAD = 1matp
#
# V="V"
# S="S"
# #y = cint(y)%DIV+((y%(2*DIV))/DIV)*(DIV-2*(y%DIV))

#Rounds x to nearest VAL
def roundToValue(x,VAL,rval):
    if rval == "u":
        return (VAL * float(int((float(x)/VAL*100))/100))+1
    else:
        return (VAL * float(int((float(x)/VAL*100))/100))

#Normalizes a numpy array - intended for converting vectors into unit vectors
def unitVector(vector):
    return (vector/np.linalg.norm(vector))

#Generates a random vector. C determines positive/neg. Vel determines if velocity
def randV(vel,c,scat):
    #Generates a vector of random components.
    if vel:
        vect = np.array([uniform(c,1) for dim in DIM])
        #Converts to unit vector
        return unitVector(vect)*VEL
    else:
        vect = np.array([uniform(c*dim,dim) for dim in DIM])

        return vect

#Generates a sinusoidal function. X is an array of x values, y and array of y values
# a is an index determining which dimension of the graph to call offsets from
#b indicates the slope of the sin
def cint(x):
    try:
        return int(x)
    except:
        return x.astype(int)

#Sin function
def sin(x,a,y):
    if RANK == 2:
        return SCATDC[a]+SCATDEPTH[a]*np.cos((x+(y*b))*SCATCT[a])
    if RANK == 3:
        if HBONE:
            y = cint(y)%DIV+DIV*((cint(y)%(2*DIV))/DIV)%2*(DIV-2*(y%DIV))
        return SCATDC[2]-SCATDC[a]+SCATDEPTH[a]*np.cos((x+(y*DIAG))*SCATCT[a])

#Generates a sloped linear function with same variable properties as the previous function
def linear(x,a,y):
    if RANK == 2:
        return SCATDC[a]+SCATDEPTH[a]*x
    if RANK ++ 3:
        return SCATDC[2]-SCATDC[a]+SCATDEPTH[a]*x

def scatterProb(p):
    global tau
    return uniform(0,1) < dt/tau[0]

def convertPosToKey(pos):

    key = []

    for val in pos:
        key = np.append(key,np.absolute(roundToValue(val,np.float(1./THERMRES),0)))

    return key

def buildGrid():
    return np.round(np.divide(np.indices((LEN*THERMRES,WID*THERMRES,DEP*THERMRES)).T.reshape(-1,3).astype(float),
            THERMRES),2)

def set_temp(body):
    if body.pos[0] <= LEN/1.3:
        body.T = T0+TGRAD
    else:
        body.T = T0

#Field class
class Field(object):

    #Particle Class
    class Prtcl(object):

        #Wraps protons around a given dimension
        def wrap(self,i,bulk):

            val = np.absolute(self.pos[i])
            if bulk:
                bulk.I +=(self.pos[i]/val)/dt

            self.pos[i] = np.absolute(val-DIM[i])

        #specularly reflects a particle off of a boundary
        def specularReflect(self,i,bulk,vert):

            #Only scatters off the top of the space IF surface scattering is on
            # if  vert and (((SSCAT and self.pos[i]<DIM[i]/2 and uniform(0,1) < (1-P)) or \
            #     (DETSCAT and self.pos[i]<DIM[i]/2 and bulk.scatgrid.isScatter(self)))):
            #     self.scatter(bulk)
            #     self.color="blue"


            self.vel[i] = -self.vel[i]

            if self.pos[i] > 105:
                self.pos[i] = 95
            self.pos[i] += 2*(int(self.pos[i])-self.pos[i])


        def update2DPos(self,E,bulk):

            self.prevPos=np.append(self.prevPos,[self.pos],axis=0)

            self.pos += self.vel*dt + self.accel
            self.color="red"

            if self.pos[0] > DIM[0] or self.pos[0] <0:
                self.wrap(0,bulk)
                # self.specularReflect(0,bulk,True)
            if self.pos[1] > DIM[1] or self.pos[1] <0:
                self.specularReflect(1,bulk,True)

            #BulkScatter
            if scatterProb(self.tau_init) and BULKSCAT:
                self.vel = randV(1,-1,False)
                if THERM:
                    bulk.thermField.exchangeHeat(self)
                self.color="blue"

            self.vel += E*dt


        def __init__(self,p,v,m,c,sf):

            self.pos =      p
            self.prevPos = [p]
            self.vel =      v
            self.mass =     m
            self.chrg =     c
            self.scatf =    sf
            self.color =    "red"
            if THERMGRAD:
                set_temp(self)
            else:
                self.T =    50.
            self.E     =    self.T*Cp
            self.dE    =    0
            self.dT    =    0
            self.tau_init = tau[0]
            self.Etau0 =     tau[0]*self.E
            self.C     =    Cp
            self.accel =    .5*E*dt**2
            if RANK == 3:
                if THERM:
                    self.update = self.update3DPosTherm
                else:
                    self.update = self.update3DPos
            else:
                self.update = self.update2DPos


    #Generates a dictionary of particles. If rand is tru, particles have random initial position. Otherwise given in a specified manner.
    def genParticles(self,sf,rand,size):
        if rand:
            return {(i+len(self.particles)):self.Prtcl(randV(0,0,False),randV(1,-1,False),Me,Ce,sf) for i in range(size)}
        else:
            return {(i+len(self.particles)):self.Prtcl(np.array([float(i%LEN),float((i/WID))/2.,1]), randV(1,-1,False),Me,Ce,sf)
                for i in range(size)}

    def plot2DParticles(self,ax,fig):
        ax.set_xlim(0,LEN)
        ax.set_ylim(0,WID)
        pos = [np.array([]) for i in range(RANK)]
        pos = [np.append(p.pos,p.color) for p in self.particles.values()]
        pos = np.transpose(pos)

        ax.scatter(np.array(pos[0]).astype(float), np.array(pos[1]).astype(float), alpha=.5,
           color=pos[2])
        plt.draw()

    def avgPos(self):

        pos = np.array([0,0,0])

        for key,val in self.particles.items():
            pos = np.add(pos,val.pos)

        return np.divide(pos,len(self.particles))

    def avgVel(self,val):
        if val == "x":
            return( np.average([elec.vel[0] for elec in self.particles.values()]))
        if val == "y":
            return np.average([elec.vel[1] for elec in self.particles.values()])
        if val == "z":
            return np.average([elec.vel[2] for elec in self.particles.values()])

    def setE(self,E):
        self.E = E

    def update(self):
        self.it +=1
        for elec in self.particles.values():
            elec.update(self.E,self)

        return

    def __init__(self,E,tau_init,dep,d_t,p,s,size):

        global   dt;           dt  = d_t
        global   P;            P   = p
        if DSRAND:
            s   = int(s*LEN*WID)
        # global   DEP;          DEP = dep
        # global   DIM;          DIM = [LEN,WID,DEP]
        # global   TOT;          TOT = LEN*WID*DEP


        if RANK == 3:
            self.plot = self.plot3DParticles
            self.thermPlot = self.plot3DTherm

        else:
            self.plot = self.plot2DParticles

        self.E   = E
        self.I   = 0
        self.tau_init = tau_init
        self.it  = 0

        if DETSCAT:
            self.scatgrid = self.ScatterGrid(s)

        if THERM:
            self.thermField = self.ThermalGrid()

        self.particles = {}
        self.particles = self.genParticles(None,True,size)
        self.potential = None

        #PRESENTATION
        self.avgPos=[]
        self.count = True

###################
### BEGIN CODE
################


t = [0]
Vx = [0];Vy = [0];Vz = [0]
Ep = [0];Eb = [0];Ea = [0]

v = [0]

i = 0


dE = [0]
grad = [0]

app = Flask(__name__)


def genConstants(P,tau,thickness):
    T=LEN/VEL
    if tau:
        P=0
    elif thickness:
        tau = 1
    dt = .01*tau

    return dt,tau,P

def calcVars(metal):
    J = metal.avgVel("x")*SIZE/(TOT/DEP)
    sigma = J/(E[0]+1e-15)
    return 1/sigma

def metalCycle(metal):
    metal.I=0
    c = metal.update()
    return c

dt = DT
metal = Field(E,tau[0],DEP,dt,0,0,SIZE)


#Defining the parameters of our plots
p1 = figure(width=500,plot_height=250, title=None, x_range=[0,LEN],y_range=[0,WID],name="line")
p2 = figure(width=500,plot_height=150, x_axis_label='Time (t)',y_axis_label='Velocity(dx/dt)',name='dot')

#Defining data structures for our plot
pos = [np.array([]) for i in range(RANK)]
pos = [np.append(p.pos,p.color) for p in metal.particles.values()]
pos = np.transpose(pos)

#Assembling out datastructures into a data dictionary for Bokeh
data = {"x":np.array(pos[0]).astype(float),\
        "y":np.array(pos[1]).astype(float),\
        "c":pos[2]}
source=ColumnDataSource(data=data)

#Datastreams for our plots
r1 = p1.circle(x="x",y="y",color="c",alpha=.5,source=source)
r2 = p2.line([],[],color="red")
r3 = p2.line([],[],color="blue")
ds1=r1.data_source
ds2=r2.data_source
ds3=r3.data_source

#Inputs for our plot in the form of sliders
sliderX = Slider(start=0, end=10, value=1, step=.1, title="X Electric Field Strength")
sliderY = Slider(start=-10, end=10, value=1, step=.1, title="Y Electric Field Strength")
slider2 = Slider(start=0.001, end=3, value=1, step=.05, title="Scattering Time")


@app.route('/')
def root():
    return render_template('index.html',resources=CDN.render())

@app.route('/initSim')
def initSim():
    #Specifying the layout we are passing
    p=gridplot([[p1],[p2]])
    layout = column(p,sliderX,sliderY,slider2)
    curdoc().add_root(layout)
    return json.dumps(json_item(layout,"esim"))
    # p = figure()
    # p.circle(1, 1)
    # item_text = json.dumps(json_item(p, "esim"))

def update():
    # ax2.text(0, -.5,"\nResistivity:"+str(np.round(rho,4)),horizontalalignment='left',verticalalignment='top')

    #Simulation content here
    t.append(t[-1]+1)
    ds2.data['x'].append(t[-1])
    ds3.data['x'].append(t[-1])
    ds2.data['y'].append(metal.avgVel("x"))
    ds3.data['y'].append(metal.avgVel("y"))

    rho =calcVars(metal)

    pos = [np.array([]) for i in range(RANK)]
    pos = [np.append(p.pos,p.color) for p in metal.particles.values()]
    pos = np.transpose(pos)

    #Updating our datastreams with new data
    ds1.data['x']=np.array(pos[0]).astype(float)
    ds1.data['y']=np.array(pos[1]).astype(float)
    ds1.data['c']=np.array(pos[2])

    #Triggering new data
    ds1.trigger('data', ds1.data, ds1.data)
    ds2.trigger('data', ds2.data, ds2.data)
    ds3.trigger('data', ds3.data, ds3.data)
    metalCycle(metal)
def blocking_task():
    while True:
        time.sleep(.1)
        curdoc().add_next_tick_callback(update())


def eX_handler(attr, old, new):
    global E
    E[0]=sliderX.value
def eY_handler(attr, old, new):
    global E
    E[1]=sliderY.value
def tau_handler(attr, old, new):
    global tau
    tau[0]=slider2.value

app.config['ENV'] = 'development'
app.config['DEBUG'] = True
app.config['TESTING'] = True

if __name__=='__main__':
    print("running")
    app.run(debug=True)

#Updating slider
# sliderX.on_change('value',eX_handler)
# sliderY.on_change('value',eY_handler)
# slider2.on_change('value',tau_handler)

#Exporting our app to the HTML doc or whatever
# curdoc().add_root(layout)
thread = Thread(target=blocking_task)
thread.start()
