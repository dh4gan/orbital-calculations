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

        self.e = eccentricityVector.mag();
        
        # Calculate Semi-latus rectum
        
        semilat = angmom*angmom/(gravparam)
        
        # Semimajor axis
        self.a = semilat/(1.0-self.e*self.e)
        
           // Calculate Orbital Inclination

    if (magOrbitalAngularMomentum > 0.0)
    {
    inclination = acos(orbitalAngularMomentum.elements[2]
        / magOrbitalAngularMomentum);

    if(orbitalAngularMomentum.elements[2] < 0.0)
        {
        inclination = twopi - inclination;
        }
    }
    else
    {
    inclination = 0.0;
    }

    // Calculate Longitude of the Ascending Node

    if (inclination == 0.0)
    {

    longitudeAscendingNode = 0.0;

    nplane.elements[0] = magOrbitalAngularMomentum;
    nplane.elements[1] = 0.0;
    nplane.elements[2] = 0.0;
    nscalar = nplane.magVector();

    }
    else
    {

    nplane.elements[0] = -orbitalAngularMomentum.elements[2];
    nplane.elements[1] = orbitalAngularMomentum.elements[1];
    nplane.elements[2] = 0.0;

    nscalar = nplane.magVector();
    longitudeAscendingNode = acos(nplane.elements[0] / nscalar);

    if (nplane.elements[1] < 0.0)
        {
        longitudeAscendingNode = twopi - longitudeAscendingNode;
        }
    }

    // Calculate true anomaly

    magpos = position.magVector();

    // If orbit circular, no inclination, then use the position vector itself

    if (eccentricity == 0.0 and inclination == 0.0)
    {
    trueAnomaly = acos(position.elements[0] / magpos);
    if (velocity.elements[0] < 0.0)
        {
        trueAnomaly = twopi - trueAnomaly;
        }

    }

    // If orbit circular and inclination non-zero, then use the orbital plane vector
    else if (eccentricity == 0.0)
    {
    ndotR = nplane.dotProduct(position);
    ndotR = ndotR / (magpos * nscalar);

    ndotV = nplane.dotProduct(velocity);

    trueAnomaly = acos(ndotR);

    if (ndotV > 0.0)
        {
        trueAnomaly = twopi - trueAnomaly;
        }
    }
    // For non-circular orbits use the eccentricity vector
    else
    {
    edotR = eccentricityVector.dotProduct(position);
    edotR = edotR / (magpos * eccentricity);

    rdotV = velocity.dotProduct(position);

    trueAnomaly = acos(edotR);

    if (rdotV < 0.0)
        {
        trueAnomaly = twopi - trueAnomaly;
        }

    }

    // Finally, calculate the longitude of periapsis - first calculate the argument of periapsis

    if (eccentricity != 0.0)
    {
    edotn = eccentricityVector.dotProduct(nplane);
    edotn = edotn / (nscalar * eccentricity);

    argumentPeriapsis = acos(edotn);

    if (eccentricityVector.elements[2] < 0.0)
        {
        argumentPeriapsis = 2.0 * pi - argumentPeriapsis;
        }

    longitudePeriapsis = argumentPeriapsis + longitudeAscendingNode;
    }

    else
    {
    argumentPeriapsis = 0.0;
    longitudePeriapsis = 0.0;
    }

    }
        
    
    def calcVectorFromOrbit(self):
        """Returns position and velocity vectors from orbital calculations"""
        
        # Calculate distance from orbital focus
        
        # Calculate position in the orbital plane
        
        # Do rotations 
        
        void Body::calcVectorFromOrbit(double G, double totmass)
    {

    /* Author:dh4gan 1/8/13
     *
     * Calculates a Body's position and velocity, given its orbital elements (relative to CoM)
     * Uses separation and true anomaly to calculate position in frame coplanar with orbit
     * Then uses a series of rotations to give correct inclination and longitudes of periapsis and ascending nodes
     *
     */

    double magpos, magvel;
    double semiLatusRectum, gravparam;

    /* 1. calculate distance from CoM using semimajor axis, eccentricity and true anomaly*/

    magpos = semiMajorAxis * (1.0 - eccentricity * eccentricity) / (1.0
        + eccentricity * cos(trueAnomaly));

    /* 2. Calculate position vector in orbital plane */

    position.elements[0] = magpos * cos(trueAnomaly);
    position.elements[1] = magpos * sin(trueAnomaly);
    position.elements[2] = 0.0;

    /* 3. Calculate velocity vector in orbital plane */
    semiLatusRectum = fabs(semiMajorAxis * (1.0 - eccentricity * eccentricity));
    gravparam = G * totmass;

    if (semiLatusRectum != 0.0)
    {
    magvel = sqrt(gravparam / semiLatusRectum);
    }
    else
    {
    magvel = 0.0;
    }

    velocity.elements[0] = -magvel * sin(trueAnomaly);
    velocity.elements[1] = magvel * (cos(trueAnomaly) + eccentricity);
    velocity.elements[2] = 0.0;

    /* 4. Begin rotations:
     * Firstly, Rotation around z axis by -argumentPeriapsis */

    if(argumentPeriapsis!=0.0)
    {
    position.rotateZ(-1 * argumentPeriapsis);
    velocity.rotateZ(-1 * argumentPeriapsis);
    }

    /* Secondly, Rotate around x by -inclination */

    if(inclination !=0.0)
    {
    position.rotateX(-1 * inclination);
    velocity.rotateX(-1 * inclination);
    }

    /* Lastly, Rotate around z by longitudeAscendingNode */

    if(longitudeAscendingNode !=0.0)
    {
    position.rotateZ(-1 * longitudeAscendingNode);
    velocity.rotateZ(-1 * longitudeAscendingNode);
    }
