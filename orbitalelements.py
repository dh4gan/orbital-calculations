# Written 30/5/14 by dh4gan
# Conversion of orbital state vectors into orbital elements
# And vice versa

import vector

pi = 3.141592654
twopi = 2.0*pi

class orbitalElements(object):
    """Set of orbital elements"""
    def __init__(self,a, e, i, longascend, argper, trueanom):
        self.a = a
        self.e = e 
        self.i = i 
        self.longascend = longascend 
        self.argper = argper 
        self.trueanom = trueanom
        
    def __str__(self):
        s= 'a= %f \ne= %f \ni= %f \nlongascend= %f \nargper= %f \n trueanom= %f \n ' % (self.a, self.e, self.i, self.longascend, self.argper, self.trueanom)
        return s

    def calcOrbitFromVector(self,position, velocity, G, totalmass):
        """Takes input state vectors and calculates orbits"""
        
        # Calculate orbital angular momentum
        
        angmomvec = position.cross(velocity)
        
        angmom = angmomvec.mag()
        # Calculate Eccentricity Vector
        
        gravparam = G * totalmass;
        magpos = position.mag();
        magvel = velocity.mag();
        vdotr = velocity.dotProduct(position);

        if (magpos == 0.0):
            eccentricityVector = Vector3D(0.0,0.0,0.0);
        else:
            eccentricityVector = position.scaleVector(magvel*magvel).subtract(vdotr)
            eccentricityVector = eccentricityVector.scaleVector(gravparam). subtract(position.unitVector())

        eccentricity = eccentricityVector.mag();
        
        # Calculate Semi-latus rectum
        
        #
        
    
    def calcVectorFromOrbit(self):
        """Returns position and velocity vectors from orbital calculations"""
        
        # Calculate distance from orbital focus
        
        # Calculate position in the orbital plane
        
        # Do rotations 