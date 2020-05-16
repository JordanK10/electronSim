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


global tau
global DT
dt = DT
metal = Field(E,tau,DEP,dt,0,0,SIZE)

global RANK

p1 = figure(width=1000,plot_height=500, title=None, x_range=[0,LEN],y_range=[0,WID])
p2 = figure(width=500,plot_height=150, x_axis_label='Time (t)',y_axis_label='Velocity(dx/dt)')

pos = [np.array([]) for i in range(RANK)]
pos = [np.append(p.pos,p.color) for p in metal.particles.values()]
pos = np.transpose(pos)
r1 = p1.circle(np.array(pos[0]).astype(float), np.array(pos[1]).astype(float), alpha=.5,color=pos[2])
r2 = p2.line([],[],color="red")
r3 = p2.line([],[],color="blue")

ds1=r1.data_source
ds2=r2.data_source
ds3=r3.data_source

p=gridplot([[p1,p2]])

@linear()
def update(step):
    # ax2.text(0, -.5,"\nResistivity:"+str(np.round(rho,4)),horizontalalignment='left',verticalalignment='top')
    ds2.data['x'].append(t[-1]+1)
    ds3.data['x'].append(t[-1]+1)
    ds2.data['y'].append(metal.avgVel("x"))
    ds3.data['y'].append(metal.avgVel("y"))

    rho =calcVars(metal)

    pos = [np.array([]) for i in range(RANK)]
    pos = [np.append(p.pos,p.color) for p in metal.particles.values()]
    pos = np.transpose(pos)

    ds1.data['x']=np.array(pos[0]).astype(float)
    ds1.data['y']=np.array(pos[1]).astype(float)

    ds1.trigger('data', ds1.data, ds1.data)
    ds2.trigger('data', ds2.data, ds2.data)
    ds1.trigger('data', ds1.data, ds1.data)
    metalCycle(metal)

curdoc().add_root(p)
curdoc().add_periodic_callback(update, 50)
