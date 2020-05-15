from space import *
from constants import *

t = [0]
Vx = [0]
Vy = [0]
Vz = [0]
v = [0]

i = 0

Ep = [0]
Eb = [0]
Ea = [0]

dE = [0]
grad = [0]

def showScatGrid(metal):

    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    scatz= metal.scatgrid.scatlist
    scatx = [scatz[i].pos[0] for i in range(len(scatz))]
    scaty = [scatz[i].pos[1] for i in range(len(scatz))]
    ax.plot(scatx,scaty, 'bo')
    ax.plot(LEN,WID, 'ro')
    ax.plot(0,0, 'ro')

    for i in range(len(scatx)):
        circle1 = plt.Circle((scatx[i], scaty[i]), 1, color='r',alpha=.25)
        ax.add_artist(circle1)


def genConstants(P,tau,thickness):
    T=LEN/VEL
    if tau:
        P=0
    elif thickness:
        tau = 1
    dt = .01*tau

    return dt,tau,P

def calcJ(metal):
    J = metal.avgVel("x")*SIZE/(TOT/DEP)
    return(J)
    # return metal.I/(A)

def calcSigma(J):
    return J/(E[0]+1e-15)

def calcRho(sigma):
    return 1/sigma

def calcVars(metal):
    J = calcJ(metal)
    sigma = calcSigma(J)
    return 1/sigma

def metalCycle(metal):
    metal.I=0
    c = metal.update()
    return c

def createAnimation():

    global tau
    global DT
    dt = DT

    metal = Field(E,tau,DEP,dt,0,0,SIZE)
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    plt.xlabel('Average Velocity (dx/dt)')
    plt.ylabel('Time (t)')

    global RANK
    if RANK==3:
        ax = fig.add_subplot(211, projection='3d')

    i = 0

    def plotMisc(metal,rho):

        plt.subplot(212)
        plt.xlabel('Average Velocity (dx/dt)')
        plt.ylabel('Time (t)')
        ax2.text(0, 6,"\nResistivity:"+str(rho),horizontalalignment='left',verticalalignment='top')
        ax.plot([-5,DIM[0]+5],[0-1,DIM[1]+1], color="white")
        ax2.plot(t,Vx, color="red")
        ax2.plot(t,Vy, color="blue")
        ax2.plot(t,Vz, color="green")
        t.append(t[-1]+1)
        Vx.append(metal.avgVel("x"))
        Vy.append(metal.avgVel("y"))
        if RANK==3:
            Vz.append(metal.avgVel("z"))

    def update(n):
        if(n>1):
            ax.cla()
            ax2.cla()
            global E
            # rho =calcVars(metal)
            # plotMisc(metal,rho)
            metal.plot(ax,fig)
            metalCycle(metal)
            plt.autoscale(False)
            # ax.set_zlim3d(-1, DEP+1)
            return fig
    # Creates and displays the animation
    anim = animation.FuncAnimation(fig, update,frames=int(FRAMES), interval=1)
    plt.show()
    # anim.save('3DTestDATAPRINT02.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])

#Determines the resistivity dependence on tau
def conductivityTest():

    #Declares data arrays
    rho = []
    x = []

    #Iterates over values of tau
    st = time.time()


    for val in np.linspace(10,10,1):
        #Generates dt, tau, and P values (only dt is relevant)
        dt,tau,P = genConstants(None,float(val),None)
        global DT
        dt = DT

        #Creates a metal nanocrystal based on constants in constants.py
        metal = Field(E,float(val),DEP,dt,P,0,SIZE)

        #Generates a set of rho values to average over
        rhoVals = []

        #Generates a number of iterations (frames/(frames/iteration)
        frames = int(FRAMES/dt*.35)

        start_time = time.time()
        #Cycles through frames number of iterations
        for j in range(frames):
            #Iterates the simulation
            metalCycle(metal)
            sys.stdout.flush()
            sys.stdout.write("\r{0}".format("Iter:"+str(float(float(j)/float(frames)))+"  :  Time:"+str(time.time()-start_time)+"   :   Resistivity:"+str(calcRho(calcSigma(calcJ(metal))))))
        print('\n\n')

        rho.append(calcVars(metal))
        x.append(1/float(val))
    print(time.time()-st)

    print(x)
    print(rho)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel("dt/tau")
    plt.ylabel("Resistivity")
    ax.plot(x,rho,'bo')
    plt.show()

def thicknessTest():
    #Declares data arrays
    rho = []
    x = []
    print(SIZE)

    for val in np.linspace(.6,1.75,5):

        #Generates dt, tau, and P values (only dt is relevant)
        global tau
        global DT
        dt = DT
        metal = Field(E,tau,float(val),dt,P,0,SIZE)

        #Generates a set of rho values to average over
        rhoVals = []

        #Generates a number of iterations (frames/(frames/iteration)
        frames = int(FRAMES/dt**.35)
        start_time = time.time()

        for j in range(frames):
            metalCycle(metal)
            sys.stdout.flush()
            sys.stdout.flush()
            sys.stdout.write("\r{0}".format("Iter:"+str(float(float(j)/float(frames)))+"  :  Time:"+str(time.time()-start_time)+"   :   Current:"+str((calcRho(calcSigma(calcJ(metal)))))))
        print('\n\n')

        rho.append(calcVars(metal))
        x.append(float(val))

    print(rho)
    print(x)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel("t")
    plt.ylabel("Resistivity")
    ax.plot(x,rho,'bo')
    plt.show()

def fermiAnim():
    #Declares data arrays

    #Generates dt, tau, and P values (only dt is relevant)
    global tau
    global DT
    global DEP
    dt = DT
    metal = Field(E,tau,DEP,dt,0,0,SIZE)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel("Bin")
    plt.ylabel("Count")
    bcktCt = 200
    distrib = np.linspace(1,bcktCt,bcktCt)

    def update(n):
        print(n)
        if(n>1):
            count = np.array([0 for i in range(bcktCt)])
            metalCycle(metal)
            for key,particle in metal.particles.iteritems():
                vel = np.around(np.linalg.norm(particle.vel),decimals=0)
                count[int(vel)-1]+=1
            ax.cla()
            plt.ylim(ymax=1)
            ax.plot(distrib,unitVector(count),'r')
            return fig

    anim = animation.FuncAnimation(fig, update,frames=FRAMES, interval=1)
    plt.show()

def fermiThickAnim():
    #Declares data arrays
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel("t")
    plt.ylabel("Resistivity")


    upLim=20
    DEP = np.linspace(0.5,5,upLim)
    bcktCt = 320
    distrib = np.linspace(1,bcktCt,bcktCt)


    def update(n):

        #Generates dt, tau, and P values (only dt is relevant)
        global tau
        global DT
        dt = DT
        metal = Field(E,tau,DEP[n],dt,P)

        #Generates buckets
        count = np.array([0 for i in range(bcktCt)])

        #Generates a number of iterations (frames/(frames/iteration)
        frames = int(FRAMES/dt*.8)
        start_time = time.time()

        for j in range(int(frames)):
            metalCycle(metal)
            # sys.stdout.flush()
            # sys.stdout.write("\r{0}".format("Iter:"+str(float(float(j)/float(frames)))))
        print('\n\n')

        for key,particle in metal.particles.iteritems():
            vel = np.around(np.linalg.norm(particle.vel),decimals=0)
            if vel < bcktCt:
                count[int(vel)-1]+=1
        ax.cla()
        z = np.polyfit(distrib, unitVector(count), 3)
        f = np.poly1d(z)
        y=f(distrib)

        plt.ylim(ymax=.6)
        plt.annotate('Thickness: '+str(DEP[n]), xy=(0,.6), xytext=(0, .6))
        plt.annotate('AvgVel: '+str(metal.avgVel("x")), xy=(150,.6), xytext=(150, .6))

        ax.plot(distrib,y,'b')

    anim = animation.FuncAnimation(fig, update,frames=upLim, interval=1)
    # plt.show()
    anim.save('FermiDistrib_t1.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])

#Varies the scatter parameter/density,  plots as function of resistivity
def surfscatTest():

    #Declares data arrays
    rho = []
    x = []

    #P is probablistic, d for determinate.
    def test(p,d):

        #Generates dt, tau, and P values (only dt is relevant)
        global tau
        global DT
        dt = DT
        #Creates a metal nanocrystal based on constants in constants.py
        metal = Field(E,tau,DEP,dt,p,d,SIZE)
        # showScatGrid(metal)
        # # plt.show()

        #Generates a number of iterations (frames/(frames/iteration)
        frames = int(FRAMES/dt*.6)

        start_time = time.time()
        print("\n")
        rhoVals = []
        for j in range(frames+100):

            metalCycle(metal)
            # sys.stdout.flush()
            # sys.stdout.write("\r{0}".format("Iter:"+str(float(float(j)/float(frames)))+"  :  Time:"+str(time.time()-start_time)+"   :   Current:"+str(calcRho(calcSigma(calcJ(metal))))))
            if j > frames:
                rhoVals.append(calcRho(calcSigma(calcJ(metal))))
        print(np.average(rhoVals))
        rho.append(np.average(rhoVals))
        if SSCAT:
            x.append(float(1-val))
        if DETSCAT:
            x.append(float(val))

    #If a probablistic or deterministic, a different variable is iterated over
    if SSCAT:
        for val in np.linspace(.00,1,5):
            test(float(val),0)
    elif DETSCAT:
        if DSRAND:
            for val in np.linspace(1,10,10):
                test(0,float(val))
        else:
            for val in np.array([i**2 for i in range(0,16)]):
                test(0,int(val))

    print(rho)
    print(x)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel("1-P")
    plt.ylabel("Resistivity")
    ax.plot(x,rho,'bo')
    plt.show()

def timingTest():
    #Declares data arrays
    Runtime = []
    Density = [1,5,10,25,50,75,100]

    #Iterates over values of tau
    st = time.time()

    for val in Density:
        #Generates dt, tau, and P values (only dt is relevant)
        global DENS;
        DENS = val
        global SIZE
        SIZE = val*TOT

        global DT
        dt = DT

        #Creates a metal nanocrystal based on constants in constants.py
        metal = Field(E,10,DEP,dt,P,0,SIZE)

        #Generates a set of rho values to average over
        rhoVals = []

        #Generates a number of iterations (frames/(frames/iteration)
        frames = int(FRAMES/dt*.35)

        start_time = time.time()
        #Cycles through frames number of iterations
        for j in range(frames):
            #Iterates the simulation
            metalCycle(metal)
            sys.stdout.flush()
            sys.stdout.write("\r{0}".format("Iter:"+str(float(float(j)/float(frames)))+"  :  Time:"+str(time.time()-start_time)+"   :   Size"+str(SIZE)))
        print('\n\n')

        run_time = time.time()-start_time

        Runtime.append(float(run_time))
    print(time.time()-st)

    print(Runtime)
    print(Density)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel("Density")
    plt.ylabel("Runtime")
    ax.plot(Density,Runtime,'bo')
    plt.show()

def precTest():
    #Declares data arrays
    Variance = []
    Density = [1,5,10,25,50,75,100]
    Runtime = []
    #Iterates over values of tau
    st = time.time()

    for val in Density:
        #Generates dt, tau, and P values (only dt is relevant)
        global DENS;
        DENS = val
        global SIZE
        SIZE = val*TOT

        global DT
        dt = DT

        #Creates a metal nanocrystal based on constants in constants.py
        metal = Field(E,10,DEP,dt,P,0,SIZE)

        #Generates a number of iterations (frames/(frames/iteration)
        frames = int(FRAMES/dt*.35)

        start_time = time.time()
        #Cycles through frames number of iterations
        for j in range(frames):
            #Iterates the simulation
            metalCycle(metal)
            sys.stdout.flush()
            sys.stdout.write("\r{0}".format("Iter:"+str(float(float(j)/float(frames)))+"  :  Time:"+str(time.time()-start_time)+"   :   Size"+str(SIZE)))
        print('\n\n')

        run_time = time.time()-start_time

        Runtime.append(float(run_time))
        Variance.append(calcVars(metal)-(1/(float(val)*10)))
        print(val,Runtime[-1],Variance[-1])
    print(time.time()-st)

    print(Runtime)
    print(Variance)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel("Density")
    plt.ylabel("Runtime")
    ax.plot(Runtime,Variance,'bo')
    plt.show()

def thermAnim():

    metal = Field(E,tau,DEP,DT,1,SCATTERDENS,SIZE)
    fig = plt.figure(figsize=(8, 6))
    grid = plt.GridSpec(4,4,wspace=0.4,hspace=0.3)
    ax = fig.add_subplot(212,projection='3d')
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(221)

    i = 0

    def calcTotE(metal):

        ep = 0
        eb = 0
        for key,value in metal.particles.items():
            ep += value.E
        for key,value in metal.thermField.field.items():
            eb += value.E
        ea =(ep + eb)
        ep = ep
        eb = eb
        return ep,eb,ea

    def plotMisc(metal):

        ep,eb,ea = calcTotE(metal)
        ax.cla()
        ax2.cla()
        ax3.cla()
        ax2.text(0, 6,"\Total Energy:"+str(ea),horizontalalignment='left',verticalalignment='top')
        ax2.plot(t,Ep, color="red")
        ax2.plot(t,Eb, color="blue")
        ax2.plot(t,Ea, color="green")
        print(dE)
        # ax3.plot(t,dE, color="red")
        ax3.plot(t,grad,color="blue")
        t.append(t[-1]+1)
        Ep.append(ep)
        Eb.append(eb)
        Ea.append(ea)
        dE.append(metal.avgdE())
        grad.append(metal.thermField.avgGrad())

        if RANK==3:
            Vz.append(metal.avgVel("z"))

    def update(n):
        if(n>1):
            global E
            # plotMisc(metal)
            # metal.plot(ax)
            metal.thermPlot(ax)

            metalCycle(metal)
            plt.autoscale(False)
            ax.set_zlim3d(-1, DEP+1)
            return fig
    # Creates and displays the animation
    anim = animation.FuncAnimation(fig, update,frames=FRAMES*2, interval=1)
    plt.show()
    # anim.save('3DThermTest1.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])

def thermPosTest():

    metal = Field(E,tau,DEP,DT,1,SCATTERDENS,SIZE)
    fig = plt.figure(figsize=(8, 6))
    grid = plt.GridSpec(4,4,wspace=0.4,hspace=0.3)
    ax = fig.add_subplot(212,projection='3d')
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(221)

    def calcTotE(metal):

        ep = 0
        eb = 0
        for key,value in metal.particles.items():
            ep += value.E
        for key,value in metal.thermField.field.items():
            eb += value.E
        ea =(ep + eb)
        ep = ep
        eb = eb
        return ep,eb,ea

    def plotMisc(metal):

        ep,eb,ea = calcTotE(metal)
        ax2.cla()
        ax3.cla()
        ax2.text(0, 6,"\Total Energy:"+str(ea),horizontalalignment='left',verticalalignment='top')
        ax2.plot(t,Ep, color="red")
        ax2.plot(t,Eb, color="blue")
        ax2.plot(t,Ea, color="green")
        # ax3.plot(t,dE, color="red")
        ax3.plot(t,grad,color="blue")
        t.append(t[-1]+1)
        Ep.append(ep)
        Eb.append(eb)
        Ea.append(ea)
        dE.append(metal.avgdE())
        grad.append(metal.thermField.avgGrad())

        if RANK==3:
            Vz.append(metal.avgVel("z"))

    def update(n):
        if(n>1):
            global E
            # plotMisc(metal)
            pos = metal.avgPos()
            ax.scatter(pos[0], pos[1], pos[2], c=[t[-1]],vmin=0, vmax=1000 ,s=100, alpha=.5,cmap=plt.cm.get_cmap('RdBu'), zdir='z', depthshade=True)

            metalCycle(metal)
            plt.autoscale(False)
            return fig
    # Creates and displays the animation
    anim = animation.FuncAnimation(fig, update,frames=FRAMES*2, interval=1)
    plt.show()
    # anim.save('3DThermTest1.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])

def thermGradTest():
    metal = Field(E,tau,DEP,DT,1,SCATTERDENS,SIZE)
    fig = plt.figure(figsize=(16,9))
    grid = plt.GridSpec(2,2,wspace=0.4,hspace=0.3)
    ax3 = fig.add_subplot(111)
    plt.ylabel('Temperature')
    plt.xlabel('x')
    t0 = [0,150,500,1500,5000]
    colors = ["red","green","blue","orange","purple"]
    labels = ["x=3","x=6","x=9"]
    x = np.linspace(0,LEN-THERMRES,LEN*THERMRES)

    electron_DE = []
    rcParams['legend.loc'] = 'best'

    def plotMisc(metal,i):
        plt.ylim(ymin=260,ymax = 335)
        electron_DE.append(metal.totalDE())


        for key,val in LOC_E.items():
            LOC_E[key]=0
        if t[-1]==t0[i]:
            T = metal.thermField.TofX()
            ax3.plot(x,T,color=colors[i],label="t="+str(t0[i])+"s")
            i+=1
            legend()
            print(T)
        t.append(t[-1]+1)


        sys.stdout.flush()
        sys.stdout.write("\r{0}".format("Time:"+str(t[-1])))
        return i

    def update(n):
        if(n>1):
            global E; global i
            i = plotMisc(metal,i)
            metalCycle(metal)
            plt.autoscale(False)
            return fig
    # Creates and displays the animation
    anim = animation.FuncAnimation(fig, update,frames=FRAMES*2, interval=1)
    plt.show()
    # anim.save('3DThermTest1.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])

def heatFlowTest():

    metal = Field(E,tau,DEP,DT,1,SCATTERDENS,SIZE)
    fig = plt.figure(figsize=(16,9))
    grid = plt.GridSpec(2,2,wspace=0.4,hspace=0.3)
    ax = fig.add_subplot(212,projection='3d')
    ax2 = fig.add_subplot(211)
    plt.ylabel('Temperature')
    plt.xlabel('x')
    t0 = [0,150,500,1500,5000]
    colors = ["red","green","blue","orange","purple"]
    labels = ["x=3","x=6","x=9"]
    x = np.linspace(0,LEN-THERMRES,LEN*THERMRES)

    dEdt =[]
    r1 = []; r2 = [];r3 = [];r4 = []
    electron_DE = []

    rcParams['legend.loc'] = 'best'

    def plotMisc(metal,i):
        plt.ylim(ymin=260,ymax = 335)

        dEdt.append(metal.heatFlowX())

        # r1.append(dEdt[-1][0])
        # r2.append(dEdt[-1][1])
        # r3.append(dEdt[-1][2])
        # electron_DE.append(metal.totalDE())

        #AVGING OVER FIVE RUNS

        interpRange = 100
        for i in range(interpRange):
            metalCycle(metal)
            dEdt.append(metal.heatFlowX())
            r1.append(dEdt[-1][0])
            r2.append(dEdt[-1][1])
            r3.append(dEdt[-1][2])
            t.append(t[-1]+1)

        metal.interpolate([r1,r2,r3],interpRange)

        ax2.cla()
        ax2.plot(t[:-1],r1,color="red",label="x=2")
        ax2.plot(t[:-1],r2,color="blue",label="x=6")
        ax2.plot(t[:-1],r3,color="green",label="x=10")
        ax2.legend()

        sys.stdout.flush()
        sys.stdout.write("\r{0}".format("Time:"+str(t[-1])))
        return i

    def update(n):
        if(n>1):
            global E; global i
            i = plotMisc(metal,i)
            metal.thermPlot(ax)
            plt.autoscale(False)
            return fig
    # Creates and displays the animation
    anim = animation.FuncAnimation(fig, update,frames=FRAMES*2, interval=1)
    plt.show()
    # anim.save('3DThermTest1.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])

def randomWalk():


    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(313)
    plt.xlabel('Average Velocity (dx/dt)')
    plt.ylabel('Time (t)')

    metal = Field(E,tau,DEP,DT,0,0,SIZE)
    metalCycle(metal)
    circle1= metalCycle(metal)
    circle1= metalCycle(metal)

    metal.plot(ax,fig)
    ax.add_artist(circle1)
    ax2.plot()
    ax2.cla()

    vel = [[metal.avgPos[-1][0]],[metal.avgPos[-1][1]]]
    t = [0]

    def plotMisc(metal,rho,vel,circle1):
        global t
        ax.plot([-5,DIM[0]+5],[0-1,DIM[1]+1], color="white",alpha=.01)

        ax2.plot(t,vel,color="r")
        ax.add_artist(circle1)
        ax.plot()

    def update(n):
        if(n>1):
            ax.cla()
            circle1 = metalCycle(metal)
            t.append(t[-1]+1)

            if metal.count:
                vel[0].append(metal.avgPos[-1][0])
                vel[1].append(metal.avgPos[-1][1])
            else:
                t.pop(-1)
            # metal.plot(ax,fig)

            rho =calcVars(metal)
            # plotMisc(metal,rho,vel,circle1)
            return fig
    # Creates and displays the animation
    anim = animation.FuncAnimation(fig, update,frames=int(FRAMES), interval=1)
    anim.save('smiley1.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])
    print(vel[0])
    print(vel[1])


if __name__ == "__main__":

    # timingTest()
    # conductivityTest()
    # thicknessTest()
    # fermiAnim()
    # fermiThickAnim()
    # surfscatTest()
    # directBlastScatterTest()
    createAnimation()
    # thermAnim()
    # thermPosTest()
    # thermGradTest()
    # heatFlowTest()
    # precTest()
    # randomWalk()
    global DEP
    print(DEP)
