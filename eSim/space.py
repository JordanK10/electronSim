from constants import *

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
