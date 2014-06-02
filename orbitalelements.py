# Written 30/5/14 by dh4gan
# Conversion of orbital state vectors into orbital elements
# And vice versa

from vector import Vector3D
import numpy as np

pi = 3.141592654
twopi = 2.0*pi
tiny = 1.0e-10

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
        s= 'a= %e \ne= %e \ni= %e \nlongascend= %e \nargper= %e \ntrueanom= %e \n ' % (self.a, self.e, self.i, self.longascend, self.argper, self.trueanom)
        return s
    
    def clone(self):
        return orbitalElements(self.a,self.e,self.i,self.longascend,self.argper,self.trueanom, self.position,self.velocity, self.G, self.totalMass)

    def calcOrbitFromVector(self):
        """Takes input state vectors and calculates orbital elements"""
        
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
        
        etot = 0.5*magvel*magvel - gravparam/magpos
        
        print "Total energy is ", etot
        
        
        
        # Semimajor axis
        try:
            self.a = self.semilat/(1.0-self.e*self.e)
            self.rper = self.a*(1.0-self.e)
        except ZeroDivisionError: # For parabolic orbits
            self.a = np.inf
            self.rper = self.semilat/2.0
            
        # Inclination

        if (self.angmom > 0.0):
            self.i = np.arccos(angmomvec.z/self.angmom)

            if(angmomvec.z < 0.0):
                self.i = twopi - self.i    
        else:
            self.i = 0.0 
            

        # Longitude of the ascending node
        
        self.longascend = np.arctan2(-angmomvec.y, angmomvec.x)
        print self.longascend
        if(self.longascend < -tiny):
            self.longascend = self.longascend + twopi

        print self.longascend, np.cos(self.longascend), np.sin(self.longascend)
        nplane = Vector3D(np.cos(self.longascend), np.sin(self.longascend), 0.0)

        print "Node Vector::", nplane
      
            
        # True anomaly 
        if(self.e>tiny):
            edotR = eccentricityVector.dot(self.position) 
            edotR = edotR / (magpos * self.e) 
            rdotV = self.velocity.dot(self.position) 
            print "edotR, rdotV", edotR, rdotV

            self.trueanom = np.arccos(edotR) 

            if (rdotV < tiny):
                self.trueanom = twopi - self.trueanom 
                    
        else:
            ndotR = nplane.dot(self.position.unitVector()) 
            ndotV = nplane.dot(self.velocity.unitVector()) 

            self.trueanom = np.arccos(ndotR) 

            if (ndotV > tiny):
                self.trueanom = twopi - self.trueanom 

        # Argument of periapsis 
        
        average_v = np.sqrt(gravparam/magpos)
        print "Velocity vs average:", magvel, average_v
        
        if(self.e>tiny):
            ncrosse = nplane.cross(eccentricityVector.unitVector())
            ndote = nplane.dot(eccentricityVector.unitVector())
            self.argper = np.arccos(ndote)
            
            ncrosse = ncrosse.dot(angmomvec)
            print ncrosse
            if(ncrosse>0.0):
                print "FLIPPED"
                self.argper =twopi - self.argper

        else:
            self.argper = 0.0

    
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

        if(self.argper>tiny):
            self.position.rotateZ(-self.argper) 
            self.velocity.rotateZ(-self.argper) 
     
        # Secondly, Rotate around x by -inclination */

        if(self.i >tiny):
            self.position.rotateX(-self.i) 
            self.velocity.rotateX(-self.i) 
     
        # Lastly, Rotate around z by self.longascend */

        if(self.longascend >tiny):
            self.position.rotateZ(-self.longascend) 
            self.velocity.rotateZ(-self.longascend) 
            
            
    def calcOrbitTrack(self, npoints):
        '''Given an input body's orbital parameters, 
        calculates x and y coordinates for
        its orbit over N points'''
        
        
        orbit = self.clone()
 #       orbit.calcOrbitFromVector()
        orbit.semilat = orbit.a*(1.0-orbit.e*orbit.e)
        if(orbit.e <1.0):
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
     
