from orbitalelements import orbitalElements
from vector import Vector3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import amin, amax, sin,cos, pi, random

G = 1.0
totalmass = 1.0
npoints = 1000
ntries = 200

magpos = 1.0
magvel = 0.7

fail = 0



for i in range(ntries):

    postheta = random.random()*pi
    #postheta = pi/2.0
    positionangle = random.random()*2.0*pi
    #positionangle = pi/2.0
    
    veltheta = random.random()*pi
    #veltheta = pi/2.0
    velangle = -2.0*pi + random.random()*4.0*pi
    #velangle = 0.0
    
    position = Vector3D(magpos*sin(postheta)*cos(positionangle), magpos*sin(postheta)*sin(positionangle),magpos*cos(postheta))
    velocity = Vector3D(magvel*sin(veltheta)*cos(velangle), magvel*sin(veltheta)*sin(velangle), magvel*cos(veltheta))
    
    orbit = orbitalElements(0.0,0.0,0.0,0.0,0.0,0.0, position,velocity, G, totalmass)
    orbit.calcOrbitFromVector()
    
    xorb,yorb,zorb = orbit.calcOrbitTrack(npoints)
    
    orbit2 = orbit.clone()
    orbit3 = orbit.clone()
    orbit4 = orbit.clone()
  
    orbit2.argper = 2.0*pi - orbit2.argper
    
    orbit3.longascend - 2.0*pi - orbit3.longascend
    
    orbit4.argper = 2.0*pi - orbit4.argper
    orbit4.longascend - 2.0*pi - orbit4.longascend
    
    xorb2,yorb2,zorb2 = orbit2.calcOrbitTrack(npoints)
    xorb3,yorb3,zorb3 = orbit3.calcOrbitTrack(npoints)
    xorb4,yorb4,zorb4 = orbit4.calcOrbitTrack(npoints)
    
    print "Orbit is \n", orbit
    
    orbit.calcVectorFromOrbit()
    
    print "Position: ",position
    print "Velocity: ", velocity
    print "--"
    print "Derived position:", orbit.position
    print "Derived Velocity:", orbit.velocity
    print "Derived Distance:", orbit.position.mag()
    print "Derived Speed:", orbit.velocity.mag()
    print "Position Error::", position.subtract(orbit.position).mag()
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
    
    fig1 = plt.figure()
    ax = fig1.add_subplot(111, projection ='3d')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
  
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    
    #Plot orbital path
    orbitline = ax.plot(xorb,yorb,zorb, color= 'black')
    orbitline2 = ax.plot(xorb2,yorb2,zorb2, color = 'red' ,ls='-', label = 'Flip argper')
    orbitline3 = ax.plot(xorb3,yorb3,zorb3, color = 'green', ls='--', label = 'Flip longascend')
    #orbitline4 = ax.plot(xorb4,yorb4,zorb4, color = 'blue', ls=':', label = 'Flip both')
    #orbitplane= ax.quiver(0.0,0.0,cos(orbit.longascend), sin(orbit.longascend),minlength = 5.0, color = 'blue')
    
    #orbitplane= ax.quiver(0.0,0.0,cos(orbit.argper), sin(orbit.argper),minlength = 5.0, color = 'red')
    #Plot a vector showing the velocity magnitude and direction
    #velocityarrow = ax.quiver(position.x, position.y, velocity.x,velocity.y, color='black')
    
    # Plot planet location
    planettry = ax.scatter(orbit.position.x,orbit.position.y,orbit.position.z, s=50, edgecolor='red', facecolor = 'none')
    planet = ax.scatter(position.x, position.y,position.z, s=50, color='green')
    
    star = ax.scatter(0.0,0.0, s=100, color='yellow')
    ax.legend(loc='best')
    #plt.show()
    
#plt.show()

print "NUMBER OF FAILS ", fail
    