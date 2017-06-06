import orbitalelements as orb
import vector

# Test script - reads in position/velocity data
# computes orbits, then modifies a and recomputes new positions/velocities

G = orb.GmsolAUday
totalmass = 1.0

x = 1.537969e01
y = -2.591931e01
z = 1.7925877e-01

position = vector.Vector3D(x,y,z)

vx = 2.68067772e-3
vy =  1.6282417e-3
vz = -9.5159225e-5

velocity = vector.Vector3D(vx,vy,vz)

orbit = orb.orbitalElements(0.0,0.0,0.0,0.0,0.0,0.0,position,velocity,G,totalmass)
orbit.calcOrbitFromVector()

orbit.a = 15.0

orbit.calcVectorFromOrbit()

print orbit.position
print orbit.velocity

