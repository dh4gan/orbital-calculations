# Written 30/5/14 by dh4gan
# Conversion of orbital state vectors into orbital elements
# And vice versa

from vector import Vector3D
import numpy as np

pi = 3.141592654
twopi = 2.0*pi

class orbitalElements(object):
    """Set of orbital elements"""
    def __init__(self,a, e, i, longascend, argper, trueanom, position,velocity, G, totalMass):
        self.a = a
        self.e = e 
        self.rper = self.a*(1.0-e)
        self.i = i 
        self.longascend = longascend 
        self.argper = argper 
        self.trueanom = trueanom
        self.angmom = 0
        
        self.position = position
        self.velocity = velocity
        self.G = G
        self.totalMass = totalMass
        
    def __str__(self):
        s= 'a= %f \ne= %f \ni= %f \nlongascend= %f \nargper= %f \ntrueanom= %f \n ' % (self.a, self.e, self.i, self.longascend, self.argper, self.trueanom)
        return s
    
    def clone(self):
        return orbitalElements(self.a,self.e,self.i,self.longascend,self.argper,self.trueanom, self.position,self.velocity, self.G, self.totalMass)

    def calcOrbitFromVector(self):
        """Takes input state vectors and calculates orbits"""
        tiny = 1.0e-10
        # Calculate orbital angular momentum

        angmomvec = self.position.cross(self.velocity)
        self.angmom = angmomvec.mag()
        
        # Calculate Eccentricity Vector
        
        gravparam = self.G * self.totalMass 
        magpos = self.position.mag() 
        magvel = self.velocity.mag() 
        vdotr = self.velocity.dot(self.position) 
        
        if (magpos == 0.0):
            eccentricityVector = Vector3D(0.0,0.0,0.0) 
        else:
            eccentricityVector = self.position.scalarmult(magvel*magvel).subtract(self.velocity.scalarmult(vdotr))
            eccentricityVector = eccentricityVector.scalarmult(1.0/gravparam).subtract(self.position.unitVector())

        self.e = eccentricityVector.mag() 
        
        # Calculate Semi-latus rectum
        
        self.semilat = self.angmom*self.angmom/(gravparam)
        print "semilat::", self.semilat
        print "Angular Momentum:: ",angmomvec
        print "Eccentricity Vector::", eccentricityVector
        
        # Semimajor axis
        try:
            self.a = self.semilat/(1.0-self.e*self.e)
            self.rper = self.a*(1.0-self.e)
        except ZeroDivisionError: # For parabolic orbits
            self.a = np.inf
            self.rper = self.semilat
            
        # Inclination

        if (self.angmom > 0.0):
            self.i = np.arccos(angmomvec.z/self.angmom)

            if(angmomvec.z < 0.0):
                self.i = twopi - self.i    
        else:
            self.i = 0.0 
            

        # Two primary cases - inclined and non-inclined orbits
        
        # Do non-inclined orbits first
        nplane = Vector3D(0.0,0.0,0.0)
        if(self.i>tiny):
            
            # Longitude of the ascending node
            nplane.x = -angmomvec.z 
            nplane.y = angmomvec.y 
            nplane.z = 0.0 
            nplane = nplane.unitVector()
            
            self.longascend = np.arccos(nplane.x) 

            if (nplane.y < 0.0):
                self.longascend = twopi - self.longascend 
                
            if(self.e>tiny):
                print self.e
                # Argument of Periapsis
                edotn = eccentricityVector.dot(nplane) 
                edotn = edotn /self.e 

                self.argper = np.arccos(edotn) 
                if (eccentricityVector.z < 0.0): self.argper = 2.0 * pi - self.argper 
                
                # True Anomaly
                edotR = eccentricityVector.dot(self.position) 
                edotR = edotR / (magpos * self.e) 
                print "edotR", edotR, magpos, self.e
                rdotV = self.velocity.dot(self.position) 

                self.trueanom = np.arccos(edotR) 

                if (rdotV < 0.0):
                    self.trueanom = twopi - self.trueanom 
            # If orbits are circular...
            else:
                self.argper = 0.0
                
                ndotR = nplane.dot(self.position.unitVector()) 
                ndotV = nplane.dot(self.velocity.unitVector()) 

                self.trueanom = np.arccos(ndotR) 

                if (ndotV > tiny):
                    self.trueanom = twopi - self.trueanom 
                
                    
        # If orbit has zero-inclination..        
        else:
            print "Zero inclination ", self.i
            self.longascend = 0.0
                
            # non-circular orbits
            if(self.e>tiny):
                
                # Argument of Periapsis 
                self.argper = np.arctan2(eccentricityVector.y,eccentricityVector.x)
                print self.argper
                if(self.argper<0.0):
                    self.argper +=twopi
                print self.argper
                if(angmomvec.z <0.0):
                    self.argper = twopi - self.argper
                    
                # True Anomaly
                edotR = eccentricityVector.dot(self.position) 
                edotR = edotR / (magpos * self.e) 
                print "edotR", edotR, magpos, self.e
                rdotV = self.velocity.dot(self.position) 
                print "rdotV", rdotV
                self.trueanom = np.arccos(edotR) 
                print "true:",self.trueanom
                if (rdotV < 0.0):
                    self.trueanom = twopi - self.trueanom 
                print "true:", self.trueanom
            # If orbits are circular
            else:
                # Argument of Periapsis
                self.argper = 0.0
                # True Anomaly
                self.trueanom = np.arccos(self.position.x / magpos) 
                if (self.velocity.x > tiny):
                    self.trueanom = twopi - self.trueanom
                print "ARGH" 
            
    
    def calcVectorFromOrbit(self):
        """Returns position and velocity vectors from orbital calculations"""
        # calculate distance from CoM using semimajor axis, eccentricity and true anomaly

        magpos = self.semilat / (1.0+ self.e * np.cos(self.trueanom)) 

        self.position = Vector3D(0.0,0.0,0.0)
        self.velocity = Vector3D(0.0,0.0,0.0)
        
        # Calculate self.position vector in orbital plane
        self.position.x = magpos * np.cos(self.trueanom) 
        self.position.y = magpos * np.sin(self.trueanom) 
        self.position.z = 0.0 

        # Calculate self.velocity vector in orbital plane */
        gravparam = self.G * self.totalMass 

        try:
            magvel = np.sqrt(gravparam/self.semilat)
        except ZeroDivisionError:
            magvel = np.sqrt(2.0*gravparam/magpos)
            
        self.velocity.x = -magvel * np.sin(self.trueanom) 
        self.velocity.y = magvel * (np.cos(self.trueanom) + self.e) 
        self.velocity.z = 0.0 

        # Begin rotations:
        # Firstly, Rotation around z axis by -self.argper */

        if(self.argper>0.0):
            self.position.rotateZ(self.argper) 
            self.velocity.rotateZ(self.argper) 
     
        # Secondly, Rotate around x by -inclination */

        if(self.i >0.0):
            self.position.rotateX(self.i) 
            self.velocity.rotateX(self.i) 
     
        # Lastly, Rotate around z by self.longascend */

        if(self.longascend >0.0):
            self.position.rotateZ(self.longascend) 
            self.velocity.rotateZ(self.longascend) 
            
            
    def calcOrbitTrack(self, npoints):
        '''Given an input body's orbital parameters, 
        calculates x and y coordinates for
        its orbit over N points'''
        
        
        orbit = self.clone()
        orbit.calcOrbitFromVector()
        
        if(self.e <=1.0):
            nu = np.linspace(0,2.0*np.pi, num=npoints)
        else:
            nu = np.linspace(-np.arccos(1.0/self.e), np.arccos(1.0/self.e), num=npoints)                
    
        x = np.zeros(npoints)
        y = np.zeros(npoints)
        z = np.zeros(npoints)
        
        for i in range(npoints):
            orbit.trueanom = nu[i]
            orbit.calcVectorFromOrbit()
            x[i] = orbit.position.x
            y[i] = orbit.position.y
            z[i] = orbit.position.z
        
        #r = self.semilat/(1.0+self.e*np.cos(nu))
        #print "semilat: ",self.semilat
        
        # Generate 3 numpy arrays based on this orbital information
            
        #x = r*(np.cos(self.longascend)*np.cos(self.argper+nu) - np.sin(self.longascend)*np.sin(self.argper+nu)*np.cos(self.i))
        #y = r*(np.sin(self.longascend)*np.cos(self.argper+nu) - np.cos(self.longascend)*np.sin(self.argper+nu)*np.cos(self.i))
        #z = r*(np.sin(self.argper+nu)* np.sin(self.i))
    
        return x,y,z    
     
