import numpy as np
import sys
from numpy import linalg as LA
from matplotlib import pyplot as plt

r0 = np.array([1.0, 1.0]) # density 
u0 = np.array([0.0, 0.0]) # velocity
p0 = np.array([3.0, 1.0]) # pressure

gamma = 5.0/3.0
g = gamma
N = 200
h = 1.0/N
tau = 0.0001 # time step
T = 0.0 # general time
q = np.zeros((3,N))
q_new = np.zeros((3,N))

def get_A(u, hs):
    return np.array([[0.0,               1.0,       0.0],
                     [-(3.0-g)*0.5*u**2, (3.0-g)*u, g-1.0],
                     [u*(-hs + (g - 1.0)*0.5*u**2), hs - (g - 1.0)*u**2, g*u]])

for i in xrange(N/2):
    E = p0[0]/(r0[0]*(gamma - 1)) + 0.5*u0[0]**2
    hs = E + p0[0]/r0[0]
    q[:,i] = np.array([r0[0], r0[0]*u0[0], r0[0]*E])
#    F[:,i] = np.array([u0[0]*r0[0], r0[0]*u0[0]**2 + p0[0], r0[0]*u0[0]*hs])

for i in xrange(N/2, N):
    E = p0[1]/(r0[1]*(gamma - 1)) + 0.5*u0[1]**2
    hs = E + p0[1]/r0[1]
    q[:,i] = np.array([r0[1], r0[1]*u0[1], r0[1]*E])
#    F[:,i] = np.array([u0[1]*r0[1], r0[1]*u0[1]**2 + p0[1], r0[1]*u0[1]*hs])

while T < 0.26:
    
    # left boundary
    
    ul = q[1,0]/q[0,0]
    u  = q[1,0]/q[0,0]
    ur = q[1,1]/q[0,1]
    
    rl = q[0,0]
    r  = q[0,0]
    rr = q[0,1]
        
    El = q[2,0]/q[0,0]
    E  = q[2,0]/q[0,0]
    Er = q[2,1]/q[0,1]
        
    pl = (El - 0.5*ul**2)*rl*(g-1.0)
    p  = (E - 0.5*u**2)*r*(g-1.0)
    pr = (Er - 0.5*ur**2)*rr*(g-1.0)
        
    hsl = (El - 0.5*ul**2) + pl/rl + ul**2/2
    hs = (E - 0.5*u**2) + p/r + u**2/2
    hsr = (Er - 0.5*ur**2) + pr/rr + ur**2/2

    A   = get_A(u, hs)
    A_r = get_A(ur,hsr)

    a  = LA.norm(A, 2)
    ar = LA.norm(A_r, 2)
    nu = np.max([a, ar])

    F_l = np.array([ul*rl, rl*ul**2 + pl, rl*ul*hsl])
    F   = np.array([u*r, r*u**2 + p, r*u*hs])
    F_r = np.array([ur*rr, rr*ur**2 + pr, rr*ur*hsr])
        
    flx_l = 0.5*(F + F_l) - 0.5*nu*(q[:,0] - q[:,0])
    flx_r = 0.5*(F + F_r) - 0.5*nu*(q[:,1] - q[:,0])

    q_new[:,0] = q[:,0] - (tau/h)*(flx_r - flx_l)
    
    for i in xrange(1,N-1):

        ul = q[1,i-1]/q[0,i-1]
        u  = q[1,i]/q[0,i]
        ur = q[1,i+1]/q[0,i+1]

        rl = q[0,i-1]
        r  = q[0,i]
        rr = q[0,i+1]

        El = q[2,i-1]/q[0,i-1]
        E  = q[2,i]/q[0,i]
        Er = q[2,i+1]/q[0,i+1]

        pl = (El - 0.5*ul**2)*rl*(g - 1.0)
        p  = (E - 0.5*u**2)*r*(g - 1.0)
        pr = (Er - 0.5*ur**2)*rr*(g - 1.0)
    
        hsl = (El - 0.5*ul**2) + pl/rl + 0.5*ul**2
        hs  = (E - 0.5*u**2) + p/r + 0.5*u**2
        hsr = (Er - 0.5*ur**2) + pr/rr + 0.5*ur**2

        A   = get_A(u, hs)
        A_r = get_A(ur,hsr)

        a  = LA.norm(A,2)
        ar = LA.norm(A_r,2)
        nu = np.max([a, ar])

        F_l = np.array([ul*rl, rl*ul**2 + pl, rl*ul*hsl])
        F   = np.array([u*r, r*u**2 + p, r*u*hs])
        F_r = np.array([ur*rr, rr*ur**2 + pr, rr*ur*hsr])
        
        flx_l = 0.5*(F + F_l) - 0.5*nu*(q[:,i] - q[:,i-1])
        flx_r = 0.5*(F + F_r) - 0.5*nu*(q[:,i+1] - q[:,i])

        q_new[:,i] = q[:,i] - (tau/h)*(flx_r - flx_l)

        if tau*nu/h >=1:
            print tau*nu/h
            sys.exit(0)
    

    # right boundary

    ul = q[1,N-2]/q[0,N-2]
    u  = q[1,N-1]/q[0,N-1]
    ur = q[1,N-1]/q[0,N-1]
    
    rl = q[0,N-2]
    r  = q[0,N-1]
    rr = q[0,N-1]

    El = q[2,N-2]/q[0,N-2]
    E  = q[2,N-1]/q[0,N-1]
    Er = q[2,N-1]/q[0,N-1]
    
    pl = (El - 0.5*ul**2)*rl*(g-1)
    p  = (E - 0.5*u**2)*r*(g-1)
    pr = (Er - 0.5*ur**2)*rr*(g-1)
        
    hsl = (El - 0.5*ul**2) + pl/rl + ul**2/2
    hs = (E - 0.5*u**2) + p/r + u**2/2
    hsr = (Er - 0.5*ur**2) + pr/rr + ur**2/2

    A   = get_A(u, hs)
    A_r = get_A(ur,hsr)

    a  = LA.norm(A, 2)
    ar = LA.norm(A_r, 2)
    nu = np.max([a, ar])

    F_l = np.array([ul*rl, rl*ul**2 + pl, rl*ul*hsl])
    F   = np.array([u*r, r*u**2 + p, r*u*hs])
    F_r = np.array([ur*rr, rr*ur**2 + pr, rr*ur*hsr])
        
    flx_l = 0.5*(F + F_l) - 0.5*nu*(q[:,N-1] - q[:,N-2])
    flx_r = 0.5*(F + F_r) - 0.5*nu*(q[:,N-1] - q[:,N-1])

    q_new[:,N-1] = q[:,N-1] - (tau/h)*(flx_r - flx_l)

    q = np.copy(q_new)

    T += tau
    print T

l = np.linspace(0,1,N)
plt.subplot(211)
plt.plot(l, q[0])
plt.subplot(212)
plt.plot(l, q[1]/q[0])
plt.show()
