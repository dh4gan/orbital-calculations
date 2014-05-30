# Written 30/5/14 by dh4gan
# Conversion of orbital state vectors into orbital elements
# And vice versa

from vector import Vector3D
import numpy as np

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
        self.angmom = 0
        
    def __str__(self):
        s= 'a= %f \ne= %f \ni= %f \nlongascend= %f \nargper= %f \ntrueanom= %f \n ' % (self.a, self.e, self.i, self.longascend, self.argper, self.trueanom)
        return s

    def calcOrbitFromVector(self,position, velocity, G, totalmass):
        """Takes input state vectors and calculates orbits"""
        
        # Calculate orbital angular momentum

        angmomvec = position.cross(velocity)
        self.angmom = angmomvec.mag()
        
        # Calculate Eccentricity Vector
        
        gravparam = G * totalmass 
        magpos = position.mag() 
        magvel = velocity.mag() 
        vdotr = velocity.dot(position) 
        
        if (magpos == 0.0):
            eccentricityVector = Vector3D(0.0,0.0,0.0) 
        else:
            eccentricityVector = position.scalarmult(magvel*magvel).subtract(velocity.scalarmult(vdotr))
            eccentricityVector = eccentricityVector.scalarmult(1.0/gravparam).subtract(position.unitVector())

        self.e = eccentricityVector.mag() 
        
        # Calculate Semi-latus rectum
        
        semilat = self.angmom*self.angmom/(gravparam)
        
        # Semimajor axis
        self.a = semilat/(1.0-self.e*self.e)
        
        # Inclination

        if (self.angmom > 0.0):
            self.i = np.arccos(angmomvec.z/self.angmom)

            if(angmomvec.z < 0.0):
                self.i = twopi - self.i    
        else:
            self.i = 0.0 
            

        #Calculate Longitude of the Ascending Node

        nplane = Vector3D(0.0,0.0,0.0)
        if (self.i == 0.0):
            self.longascend = 0.0 

            nplane.x = self.angmom 
            nplane.y = 0.0
            nplane.z = 0.0 
            nscalar = nplane.mag() 

    
        else:
            nplane.x = -angmomvec.z 
            nplane.y = angmomvec.y 
            nplane.z = 0.0 

            nscalar = nplane.mag() 
            self.longascend = np.arccos(nplane.x / nscalar) 

            if (nplane.y < 0.0):
                self.longascend = twopi - self.longascend 

        # Calculate true anomaly

        magpos = position.mag() 

        # If orbit circular, no inclination, then use the position vector itself

        if (self.e == 0.0 and self.i == 0.0):
            self.trueanom = np.arccos(position.x / magpos) 
            if (velocity.x < 0.0):
                self.trueanom = twopi - self.trueanom 

        # If orbit circular and inclination non-zero, then use the orbital plane vector
        elif (self.e == 0.0):
     
            ndotR = nplane.dot(position) 
            ndotR = ndotR / (magpos * nscalar) 

            ndotV = nplane.dot(velocity) 

            self.trueanom = np.arccos(ndotR) 

            if (ndotV > 0.0):
                self.trueanom = twopi - self.trueanom 
         
     
        # For non-circular orbits use the eccentricity vector
        else:
            edotR = eccentricityVector.dot(position) 
            edotR = edotR / (magpos * self.e) 

            rdotV = velocity.dot(position) 

            self.trueanom = np.arccos(edotR) 

            if (rdotV < 0.0):
                self.trueanom = twopi - self.trueanom 

        # Finally, calculate the longitude of periapsis - first calculate the argument of periapsis

        if (self.e != 0.0):
            edotn = eccentricityVector.dot(nplane) 
            edotn = edotn / (nscalar * self.e) 

            self.argper = np.arccos(edotn) 

            if (eccentricityVector.z < 0.0): self.argper = 2.0 * pi - self.argper 

        else:
            self.argper = 0.0         
    
    def calcVectorFromOrbit(self, G, totalmass):
        """Returns position and velocity vectors from orbital calculations"""
        # calculate distance from CoM using semimajor axis, eccentricity and true anomaly

        magpos = self.a * (1.0 - self.e * self.e) / (1.0+ self.e * np.arccos(self.trueanom)) 

        position = Vector3D(0.0,0.0,0.0)
        velocity = Vector3D(0.0,0.0,0.0)
        # Calculate position vector in orbital plane

        position.x = magpos * np.cos(self.trueanom) 
        position.y = magpos * np.sin(self.trueanom) 
        position.z = 0.0 

        # Calculate velocity vector in orbital plane */
        semiLatusRectum = abs(self.a * (1.0 - self.e * self.e)) 
        gravparam = G * totalmass 

        if (semiLatusRectum > 0.0 or semiLatusRectum<0.0):
            magvel = np.sqrt(gravparam / semiLatusRectum) 
        else:
            magvel = 0.0 

        velocity.x = -magvel * np.sin(self.trueanom) 
        velocity.y = magvel * (np.cos(self.trueanom) + self.e) 
        velocity.z = 0.0 

        # Begin rotations:
        # Firstly, Rotation around z axis by -self.argper */

        if(self.argper!=0.0):
            position.rotateZ(-1 * self.argper) 
            velocity.rotateZ(-1 * self.argper) 
     
        # Secondly, Rotate around x by -inclination */

        if(self.i !=0.0):
            position.rotateX(-1 * self.i) 
            velocity.rotateX(-1 * self.i) 
     
        # Lastly, Rotate around z by self.longascend */

        if(self.longascend !=0.0):
            position.rotateZ(-1 * self.longascend) 
            velocity.rotateZ(-1 * self.longascend) 
            
    def calcOrbitTrack(self, G, totalmass, npoints):
        '''Given an input body's orbital parameters, 
        calculates x and y coordinates for
        its orbit over N points'''
        
        if(self.e <=1.0):
            nu = np.linspace(0,2.0*np.pi, num=npoints)
        else:
            nu = np.linspace(-np.arccos(1.0/self.e), np.arccos(1.0/self.e), num=npoints)                
    
        semilat = self.angmom*self.angmom/(G*totalmass)
    
        r = semilat/(1.0+self.e*np.cos(nu))
    
        # Generate 3 numpy arrays based on this orbital information
            
        x = r*(np.cos(self.longascend)*np.cos(self.argper+nu) - np.sin(self.longascend)*np.sin(self.argper+nu)*np.cos(self.i))
        y = r*(np.sin(self.longascend)*np.cos(self.argper+nu) - np.cos(self.longascend)*np.sin(self.argper+nu)*np.cos(self.i))
        z = r*(np.sin(self.argper+nu)* np.sin(self.i))
    
        return x,y,z    
     
