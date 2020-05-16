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

def scatterProb(tau):
    return uniform(0,1) < dt/tau

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
            if scatterProb(self.tau) and BULKSCAT:
                self.vel = randV(1,-1,False)
                if THERM:
                    bulk.thermField.exchangeHeat(self)
                self.color="blue"

            self.vel += E*dt

        def update3DPos(self,E,bulk):

            self.pos += self.vel*dt + .5*E*dt**2
            self.color="red"

            for i in range(RANK):
                while self.pos[0] >= DIM[0] or self.pos[0] <=0:
                    self.wrap(0,bulk)
                while self.pos[1] >= DIM[1] or self.pos[1] <=0:
                    self.wrap(1,bulk)
                if self.pos[2] >= DIM[2] or self.pos[2] <=0:
                    self.specularReflect(2,bulk,True)

            self.bulkScatter(E,bulk)
            self.vel += E*dt

        def update3DPosTherm(self,E,bulk):

            self.prevPos = self.pos[0]
            self.pos += self.vel*dt + .5*E*dt**2
            self.color="red"

            for i in range(RANK):
                while self.pos[0] >= DIM[0] or self.pos[0] <=0:
                    self.specularReflect(0,bulk,False)
                while self.pos[1] >= DIM[1] or self.pos[1] <=0:
                    self.specularReflect(1,bulk,False)
                if self.pos[2] >= DIM[2] or self.pos[2] <=0:
                    self.specularReflect(2,bulk,True)

            self.bulkScatter(E,bulk)
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
            self.tau =      tau
            self.Etau0 =     tau*self.E
            self.C     =    Cp
            self.accel =    .5*E*dt**2
            if RANK == 3:
                if THERM:
                    self.update = self.update3DPosTherm
                else:
                    self.update = self.update3DPos
            else:
                self.update = self.update2DPos

    #Scatterer grid contains the scatterers. Particles land on the
    class ScatterGrid(object):
        #Scatterer Class
        class Scatterer(object):

            def __init__(self,pos,radius):

                self.radius = radius
                self.pos    = pos

        #Determines if a given particle is incident on a scatterer
        def isScatter(self,prtcl):

            string = "REPORT--"
            #Rounds the location to the nearest unit of the grid (.5)
            for roundVal in [["u","u"],["u","d"],["d","u"],["d","d"]]:

                x = roundToValue(prtcl.pos[0],SRAD,roundVal[0])
                y = roundToValue(prtcl.pos[1],SRAD,roundVal[1])

                #Checks dictionary of scatterers at the rounded gridpoint
                if y in self.scatterers.get(x,{}):

                    #Iterates through the associated array to check if particle is within a scatterer's radius
                    for scat in self.scatterers[x][y]:
                        if np.linalg.norm(scat.pos-prtcl.pos[0:2]) <= SRAD+.05:
                            return True

            return False

        def addScatterer(self,loc,rad):

            self.appendScatterer(self.Scatterer(loc,rad))

        #Appends a scatterer to the dictionary depending on rounded grid point
        def appendScatterer(self,scat):

            #Rounds the gridpoint
            for roundVal in [["u","u"],["u","d"],["d","u"],["d","d"]]:
                x = roundToValue(scat.pos[0],SRAD,roundVal[0])
                y = roundToValue(scat.pos[1],SRAD,roundVal[1])

                #Fills in radius around the gridpoint
                for xx in range(-1,1,1):
                    for yy in range(-1,1,1):
                        xxx = x + xx*SRAD
                        yyy = y + yy*SRAD

                            #Appends the scatterer. Mod operation for wrapping
                        try:
                            if not scat in self.scatterers[float(np.absolute((xxx)))][float(np.absolute(yyy))]:
                                try:
                                    self.scatterers[float(np.absolute((xxx)))][float(np.absolute(yyy))].append(scat)
                                except:
                                    continue
                        except:
                            continue

        #Initializes the scatter grid
        def __init__(self,num):

            self.scatterers = {i*SRAD:{j*SRAD:[] for j in range(LEN/SRAD+1)
                                       } for i in range(WID/SRAD+1)}
            self.scatlist = []

            #Fills grid with designated number of scatterers at random locations

            # Divides grid up into nxn cells where n is the celing of the square root of the number of scatterers to be input
            if not DSRAND:

                grdLen = int(np.ceil(np.sqrt(num)))

                try:
                    grdIndex = 1/float(grdLen-1)
                except:
                    grdIndex = 0

            for j in range(num):

                # Inserts to random location
                if DSRAND:
                    x = uniform(0,LEN); y = uniform(0,WID)
                # Inserts to location on grid
                else:
                    x = (((j/(grdLen)))*grdIndex)*LEN
                    y = (((j%(grdLen)))*grdIndex)*WID

                scat = self.Scatterer([x,y],SRAD)
                self.appendScatterer(scat)
                self.scatlist.append(scat)

    #Thermal field contains different grid buckets on the lattice with unique energy
    class ThermalGrid(object):


        #Calculates the average temperature at every cross section given by X
        def TofX(self):

            #Generates array of cross sectional slices set by THERMRES
            T = [i/THERMRES for i in range(LEN*THERMRES)]

            grid = buildGrid()

            #Iterates through every member of the thermgrid dictionary, summing temperatures
            for point in grid:
                T[int(point[0]*THERMRES)] += self.field[repr(point)].T

            #Returns the avg by dividing by cross sectional area
            return(np.divide(np.array(T),WID*DEP*THERMRES**2))

        def avgGrad(self):

            E = [0]

            for key,value in self.field.items():
                E.append(value.E)

            np.sort(E)

            x = 0
            size = len(E)
            for i in range(size/2):
                x += np.absolute(E[i]-E[size-i-1])/6

            return x


        #A bucket of uniform energy
        class ThermalPoint(object):

            def updateThermValues(self,item,T):
                item.dT = T-item.T
                item.dE = item.dT*item.C
                item.T = T
                item.E = T*item.C

            #Averages the energy when an electron interacts with it
            def interact(self, electron):

                T = np.absolute((self.E+electron.E)/(Cb+Cp))

                if FIXEDEDGES:
                    if not (self.pos[0]==0.): #or self.pos[0]==LEN-1./float(THERMRES)):
                        self.updateThermValues(self,T)
                else:
                        self.updateThermValues(self,T)

                self.updateThermValues(electron,T)

                if TAUtempDEP:
                    electron.tau = float(electron.Etau0/float(electron.E*tau))


            #Initializes energy
            def __init__(self, loc):

                self.pos = loc
                if THERMGRAD:
                    set_temp(self)
                else:
                    self.T = 0
                self.E = self.T*Cb
                self.dE = 0
                self.dT = 0
                self.C = Cb

        def exchangeHeat(self,electron):

            coord = convertPosToKey(electron.pos)

            try:
                self.field[repr(coord)].interact(electron)
            except:
                return

        #Generates the grid of energy buckets
        def genThermalGrid(self):

            j = 0

            BlockArr = {}

            grid = buildGrid()

            for point in grid:
                BlockArr[repr(point)] = self.ThermalPoint(point)

            return BlockArr

        def __init__(self):

            self.size = TOT*(THERMRES**3)
            self.field = self.genThermalGrid()

    def interpolate(self,values,rng):

        r = []

        for ar in values:

            p = np.poly1d(np.polyfit(np.linspace(0,rng-1,rng),
                          [ar[int(val)] for val in np.linspace(-rng,-1,rng)],1))
            for val in np.linspace(0,rng-1,rng):
                print(val,p(val),val-rng)
                ar[int(val-rng)] = p(val)
            print("\n")

    def avgdE(self):

        e = 0

        for ID,thermPt in self.thermField.field.items():
            for bound,dE in LOC_E.items():
                if thermPt.pos[0] <= int(bound):
                    LOC_E[bound] += thermPt.dE
                else:
                    LOC_E[bound] -= thermPt.dE
            e += thermPt.dE
        for ID,particle in self.particles.items():
            for bound,dE in LOC_E.items():
                if particle.pos[0] <= int(bound):
                    LOC_E[bound] += particle.dE
                else:
                    LOC_E[bound] -= particle.dE
            e += particle.dE

        return [LOC_E["2"],LOC_E["6"],LOC_E["10"],e]

    def heatFlowX(self):

        heat1 = 0; heat2 = 0;heat3 = 0; heat = 0;

        #Tracking particles across different cross sectional areas

        for key,electron in self.particles.items():
            vel = electron.vel[0]
            if int(electron.pos[0]) == 2:
                heat1 += electron.E*vel/np.absolute(vel)
            if int(electron.pos[0]) == 6:
                heat2 += electron.E*vel/np.absolute(vel)
            if int(electron.pos[0]) == 9:
                heat3 += electron.E*vel/np.absolute(vel)
            heat += 0


        return heat1, heat2, heat3, heat

    def totalDE(self):

        var = 0

        for ID,part in self.particles.items():
            var += part.dE

        return var

    #Generates a dictionary of particles. If rand is tru, particles have random initial position. Otherwise given in a specified manner.
    def genParticles(self,sf,rand,size):
        if rand:
            return {(i+len(self.particles)):self.Prtcl(randV(0,0,False),randV(1,-1,False),Me,Ce,sf) for i in range(size)}
        else:
            return {(i+len(self.particles)):self.Prtcl(np.array([float(i%LEN),float((i/WID))/2.,1]), randV(1,-1,False),Me,Ce,sf)
                for i in range(size)}

    def plot3DParticles(self,ax):

        pos = [np.array([]) for i in range(RANK)]
        for p in self.particles.values():
            for i in range(RANK):
                pos[i] = np.append(pos[i],p.pos[i])

        ax.scatter(pos[0], pos[1], pos[2], zdir='z', alpha=.5, depthshade=True,
                   color=[p.color for p in self.particles.values()])
        # if SURFSCAT:
        #     ax.plot_surface(self.scatterFuncPlot[0],self.scatterFuncPlot[1],
        #                     self.scatterFuncPlot[2], alpha=.1,color="green")

    def plot3DTherm(self,ax):

        pos = [np.array([]) for i in range(RANK)]
        clr = np.array([])
        for key,therm in (self.thermField.field).items():
            for i in range(RANK):
                pos[i] = np.append(pos[i], therm.pos[i])
            clr = np.append(clr,therm.T)

        os.system('cls' if os.name == 'nt' else 'clear')
        d=self.thermField.field
        # print(np.array([d[i].T for i in d]))
        # print('\n\n\n\n')
        # print(np.array([d[i].E for i in d]))
        ax.clear()
        sct = ax.scatter(pos[0], pos[1], pos[2], c=clr, s=100, vmin=T0, vmax=T0+TGRAD, alpha=.25,cmap=plt.cm.get_cmap('RdBu'),
         zdir='z', depthshade=True)

        global LGD
        if not LGD:
            plt.colorbar(sct)
            LGD = True

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
        #PLOTTING FOR THE PRESENTATION
        # self.count = True
        # if self.it > 75 and self.it < 150:
        #     if self.it % 20 == 0:
        #         self.particles.update(self.genParticles(None,True,1))
        #         self.count=False
        # if self.it >= 150 and len(self.particles )< 500:
        #     if self.it % 3 == 0:
        #         self.particles.update(self.genParticles(None,True,1))
        #         self.count = False
        # if self.it == 650:
        #     global E
        #     E[1]=1
        #     global SMIL
        #     SMIL = False
        #
        # r = np.transpose([p.pos for p in self.particles.items()])
        # x = r[0]
        # y = r[1]
        #
        # self.avgPos.append([np.average(x),np.average(y)])
        # return plt.Circle((np.average(x),np.average(y)), 1.5, color='b',zorder=1)


    def __init__(self,E,tau,dep,d_t,p,s,size):

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
        self.tau = tau
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

    # def updateForces(self):
    #
    #     for IDb0,block0 in self.space.items():
    #
    #         for IDp0,ptcl0 in block0.particles.items():
    #             ptcl0.fmag = 0
    #             ptcl0.fdir = ptcl0.fdir - ptcl0.fdir
    #
    #             for IDb1,block1 in self.space.items():
    #
    #                 if not IDb1==IDb0:
    #                     rdiff = np.subtract(block1.potential.pos,ptcl0.pos)
    #                     Veff = block1.potential.mag/np.linalg.norm(rdiff)
    #                     ptcl0.fmag += Veff
    #                     ptcl0.fdir += np.divide(block1.potential.mag,rdiff)
    #
    #             for IDp1,ptcl1 in block0.particles.items():
    #
    #                 if not IDp0==IDp1:
    #                     ptcl0.fmag += Ce/np.linalg.norm(np.subtract(ptcl1.pos,ptcl0.pos))
    #                     ptcl0.fdir += np.divide(block1.potential.mag,
    #                                             np.subtract(ptcl1.pos,ptcl0.pos))
    #
    #             ptcl0.computeForce()

    # def genScatterFunc(self):
    #         if SCATPAT =="periodic":
    #             return np.array([sin,sin])
    #         if SCATPAT =="diagperiodic":
    #             return np.array([sin,linear])
    #
    # def genScatterFuncPlot(self):
    #     if RANK == 3:
    #         X=np.arange(0,LEN,.25)
    #         Y=np.arange(0,WID,.25)
    #
    #         X,Y = np.meshgrid(X,Y)
    #         rng = np.add(self.scatterFunc[0](X,0,Y),self.scatterFunc[1](Y,1,0))
    #         return [X,Y,rng]
    #
    #     if RANK == 2:
    #         x = (np.lignspace(0,LEN*WID,LEN*WID*RES))%LEN
    #         rng = self.scatterFunc[0](x,0)
    #         return [x,rng]
    #
    # def setScatterFuncPlot(self):
        # self.scatterFuncPlot = self.genScatterFuncPlot()



"""
perform FFT on the time series plot of average velocity

Reproduce:
    Produce a conductive constant
    Linear relationship between characteristic time and conductivity
    Fuchs Model

"""
# array([0.75, 1.  , 0.  ])
