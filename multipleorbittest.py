from orbitalelements import orbitalElements
from vector import Vector3D
import matplotlib.pyplot as plt
from numpy import amin, amax, sin,cos, pi

G = 1.0
totalmass = 1.0
npoints = 1000

positionangle = 315.0*pi/180.0
velangle = 180.0*pi/180.0

magpos = 1.0
velmin = 0.5
velmax = 1.5

nvel = 10
dvel = (velmax- velmin)/float(nvel)

for i in range(nvel):

    magvel = velmin + i*dvel
    position = Vector3D(magpos*cos(positionangle), magpos*sin(positionangle),0.0)
    velocity = Vector3D(magvel*cos(velangle), magvel*sin(velangle), 0.0)
    
    orbit = orbitalElements(0.0,0.0,0.0,0.0,0.0,0.0, position,velocity, G, totalmass)
    orbit.calcOrbitFromVector()
    
    xorb,yorb,zorb = orbit.calcOrbitTrack(npoints)
    
    print "Orbit is \n", orbit
    
    orbit.calcVectorFromOrbit()
    
    print "Position: ",position
    print "Velocity: ", velocity
    print "--"
    print "Derived position:", orbit.position
    print "Derived Velocity:", orbit.velocity
    print "Derived Distance:", orbit.position.mag()
    print "Derived Speed:", orbit.velocity.mag()
    
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