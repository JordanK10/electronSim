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
metal = Field(E,tau[0],DEP,dt,0,0,SIZE)

global RANK

p1 = figure(width=500,plot_height=250, title=None, x_range=[0,LEN],y_range=[0,WID])
p2 = figure(width=500,plot_height=150, x_axis_label='Time (t)',y_axis_label='Velocity(dx/dt)')

pos = [np.array([]) for i in range(RANK)]
pos = [np.append(p.pos,p.color) for p in metal.particles.values()]
pos = np.transpose(pos)
data = {"x":np.array(pos[0]).astype(float),\
        "y":np.array(pos[1]).astype(float),\
        "c":pos[2]}
source=ColumnDataSource(data=data)
r1 = p1.circle(x="x",y="y",color="c",alpha=.5,source=source)
r2 = p2.line([],[],color="red")
r3 = p2.line([],[],color="blue")

ds1=r1.data_source
ds2=r2.data_source
ds3=r3.data_source

sliderX = Slider(start=0, end=10, value=1, step=.1, title="X Electric Field Strength")
sliderY = Slider(start=-10, end=10, value=1, step=.1, title="Y Electric Field Strength")
slider2 = Slider(start=0.001, end=3, value=1, step=.05, title="Scattering Time")

p=gridplot([[p1],[p2]])
layout = column(p,sliderX,sliderY,slider2)
@linear()
def update(step):
    # ax2.text(0, -.5,"\nResistivity:"+str(np.round(rho,4)),horizontalalignment='left',verticalalignment='top')
    t.append(t[-1]+1)
    ds2.data['x'].append(t[-1])
    ds3.data['x'].append(t[-1])
    ds2.data['y'].append(metal.avgVel("x"))
    ds3.data['y'].append(metal.avgVel("y"))

    rho =calcVars(metal)

    pos = [np.array([]) for i in range(RANK)]
    pos = [np.append(p.pos,p.color) for p in metal.particles.values()]
    pos = np.transpose(pos)
    ds1.data['x']=np.array(pos[0]).astype(float)
    ds1.data['y']=np.array(pos[1]).astype(float)
    ds1.data['c']=np.array(pos[2])

    ds1.trigger('data', ds1.data, ds1.data)
    ds2.trigger('data', ds2.data, ds2.data)
    ds3.trigger('data', ds3.data, ds3.data)
    metalCycle(metal)

def eX_handler(attr, old, new):
    global E
    E[0]=sliderX.value
def eY_handler(attr, old, new):
    global E
    E[1]=sliderY.value
def tau_handler(attr, old, new):
    global tau
    tau[0]=slider2.value

sliderX.on_change('value',eX_handler)
sliderY.on_change('value',eY_handler)
slider2.on_change('value',tau_handler)


curdoc().add_root(layout)
curdoc().add_periodic_callback(update, 100)
