import numpy as np
import matplotlib.pyplot as plt

grid = np.linspace(1.8181,20,11)
print(grid)


LEN = 20
WID = 20
num = 121

grdLen = int(np.ceil(np.sqrt(num)))
grdIndex = 1/float(grdLen-1)
fig, ax = plt.subplots()


for j in range(num):

    x = (((j/(grdLen)))*grdIndex)*LEN
    y = (((j%(grdLen)))*grdIndex)*WID
    print(x,y)
    circle1 = plt.Circle((x,y), 1, color='r',alpha=.2)
    ax.add_artist(circle1)

plt.xlim(xmax=20)
plt.ylim(ymax=20)

plt.show()
