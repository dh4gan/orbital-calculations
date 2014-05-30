from orbitalelements import orbitalElements
from vector import Vector3D
import matplotlib.pyplot as plt
from numpy import amin, amax

G = 1.0
totalmass = 1.0
npoints = 100

position = Vector3D(1.0,0.0,0.0)
velocity = Vector3D(0.0,1.0,0.0)

orbit = orbitalElements(0.0,0.0,0.0,0.0,0.0,0.0)
orbit.calcOrbitFromVector(position, velocity, G, totalmass)

xorb,yorb,zorb = orbit.calcOrbitTrack(G,totalmass, npoints)

print "Orbit is \n", orbit
xmin = amin(xorb)
xmin = xmin - 0.1*abs(xmin)
xmax = amax(xorb)
xmax = xmax + 0.1*abs(xmax)

ymin = amin(yorb)
ymin = ymin - 0.1*abs(ymin)
ymax = amax(yorb)
ymax = ymax + 0.1*abs(ymax)



fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('x')
ax.set_ylabel('y')

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

#Plot orbital path
orbitline = ax.plot(xorb,yorb)
#Plot a vector showing the velocity magnitude and direction
velocityarrow = ax.quiver(position.x, position.y, velocity.x,velocity.y, color='black')
# Plot planet location
planet = ax.scatter(position.x, position.y, s=50, color='red')
star = ax.scatter(0.0,0.0, s=100, color='yellow')
plt.show()