from orbitalelements import orbitalElements
from vector import Vector3D
import matplotlib.pyplot as plt
from numpy import amin, amax, sin,cos, pi, random

G = 1.0
totalmass = 1.0
npoints = 1000
ntries = 200

magpos = 1.0
magvel = 1.0

fail = 0

for i in range(ntries):

    positionangle = random.random()*2.0*pi
    velangle = -2.0*pi + random.random()*4.0*pi

    position = Vector3D(magpos*cos(positionangle), magpos*sin(positionangle),0.0)
    velocity = Vector3D(magvel*cos(velangle), magvel*sin(velangle), 0.0)
    
    orbit = orbitalElements(0.0,0.0,0.0,0.0,0.0,0.0, position,velocity, G, totalmass)
    orbit.calcOrbitFromVector()
    
    xorb,yorb,zorb = orbit.calcOrbitTrack(npoints)
    
    orbit2 = orbit.clone()
    orbit2.argper = 2.0*pi - orbit2.argper
    
    xorb2,yorb2,zorb2 = orbit2.calcOrbitTrack(npoints)
    
    print "Orbit is \n", orbit
    
    orbit.calcVectorFromOrbit()
    
    print "Position: ",position
    print "Velocity: ", velocity
    print "--"
    print "Derived position:", orbit.position
    print "Derived Velocity:", orbit.velocity
    print "Derived Distance:", orbit.position.mag()
    print "Derived Speed:", orbit.velocity.mag()
    
    if position.subtract(orbit.position).mag() >1.0e-4:
        print "Position FAIL"
        fail +=1
    if velocity.subtract(orbit.velocity).mag() >1.0e-4:
        print "Velocity FAIL"
    print "--"
    
    xmin = amin(xorb)
    xmin = xmin - 0.1*abs(xmin)
    xmax = amax(xorb)
    xmax = xmax + 0.1*abs(xmax)
    
    ymin = amin(yorb)
    ymin = ymin - 0.1*abs(ymin)
    ymax = amax(yorb)
    ymax = ymax + 0.1*abs(ymax)
    
    if(xmin > position.x): xmin = position.x - 1.0
    if(xmax < position.x): xmax = position.x + 1.0
    
    if(ymin > position.y): ymin = position.y - 1.0
    if(ymax < position.y): ymax = position.y + 1.0
    
    #fig1 = plt.figure()
    #ax = fig1.add_subplot(111)
    #ax.set_xlabel('x')
    #ax.set_ylabel('y')
  
    #ax.set_xlim(xmin,xmax)
    #ax.set_ylim(ymin,ymax)
    
    #Plot orbital path
    #orbitline = ax.plot(xorb,yorb)
    #orbitline2 = ax.plot(xorb2,yorb2, ls='dashed')
    #Plot a vector showing the velocity magnitude and direction
    #velocityarrow = ax.quiver(position.x, position.y, velocity.x,velocity.y, color='black')
    
    # Plot planet location
    #planettry = ax.scatter(orbit.position.x,orbit.position.y, s=50, edgecolor='red', facecolor = 'none')
    #planet = ax.scatter(position.x, position.y, s=50, color='green')
    
    #star = ax.scatter(0.0,0.0, s=100, color='yellow')
    #plt.show()
    
#plt.show()

print "NUMBER OF FAILS ", fail
    