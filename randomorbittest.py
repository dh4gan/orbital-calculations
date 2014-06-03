from orbitalelements import orbitalElements
from vector import Vector3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import amin, amax, sin,cos, pi, random

G = 1.0
totalmass = 1.0
npoints = 1000
ntries = 500

magpos = 1.0
magvel = 0.7

fail = 0

for i in range(ntries):

    # Randomly select spherical coordinate angles for position and velocity vectors
    postheta = random.random()*pi
    positionangle = random.random()*2.0*pi
    
    veltheta = random.random()*pi
    velangle = -2.0*pi + random.random()*4.0*pi
    
    # Use this code to fix the vector instead
    #postheta = pi/2.0
    #positionangle = pi/2.0
    #veltheta = pi/2.0
    #velangle = 0.0
    
    # Set positions and velocities
    position = Vector3D(magpos*sin(postheta)*cos(positionangle), magpos*sin(postheta)*sin(positionangle),magpos*cos(postheta))
    velocity = Vector3D(magvel*sin(veltheta)*cos(velangle), magvel*sin(veltheta)*sin(velangle), magvel*cos(veltheta))
    
    # Calculate the orbital elements
    orbit = orbitalElements(0.0,0.0,0.0,0.0,0.0,0.0, position,velocity, G, totalmass)
    orbit.calcOrbitFromVector()
    
    # Plot the orbital path
    xorb,yorb,zorb = orbit.calcOrbitTrack(npoints)
    
    print "--"
    print "Orbit is \n", orbit
    
    # Rederive the position and velocity vectors as a check
    orbit.calcVectorFromOrbit()
    
    print "--"
    print "Position: ",position
    print "Velocity: ", velocity
    print "--"
    print "Derived position:", orbit.position
    print "Derived Velocity:", orbit.velocity
    
    if position.subtract(orbit.position).mag() >1.0e-4:
        print "Position FAIL"
        fail +=1
    if velocity.subtract(orbit.velocity).mag() >1.0e-4:
        print "Velocity FAIL"
    print "--"
    
    # Plot to a graph (define x and y limits first)
    # Begin commenting out of code here if you want to just do lots of tests
    
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
    
    fig1 = plt.figure()
    ax = fig1.add_subplot(111, projection ='3d')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
  
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    
    #Plot orbital path
    orbitline = ax.plot(xorb,yorb,zorb, color= 'green')

    # Plot planet location, and where code thinks planet is
    planettry = ax.scatter(orbit.position.x,orbit.position.y,orbit.position.z, s=50, edgecolor='red', facecolor = 'none')
    planet = ax.scatter(position.x, position.y,position.z, s=50, color='green')
    
    star = ax.scatter(0.0,0.0, s=100, color='yellow')

    plt.show()
    # End commenting out of code here if you just want to do lots of tests

print "NUMBER OF FAILS ", fail
    