#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 15:15:59 2018

@author: jianinglu
"""

import math
from   matplotlib.pyplot import *   # http://matplotlib.sourceforge.net/api/pyplot_api.html
            	# http://numpy.scipy.org/
from   numpy import *
import  scipy as sc                   # http://scipy.org/
from scipy import constants  



#simulate the motion of 16 particles

#Set the initial positions of the particles to be evenly spaced in a square 
#domain of side length L

#define # of particles and the side lengths
N = 16
Lx = 4.0
Ly = 4.0
#Since each particle occupies a small cube,
#define dx and dy to be side lengths of the small cube, which are side length/sauqreroot of N

dx = Lx/sqrt(N)
dy = Ly/sqrt(N)

#x_grid and y_grid are arrays starting from half of the first small cube side length,
#the inidividual values go up by Lx AND ly, all the way to half small cube length
#from the end of the large cube
x_grid = arange(dx/2, Lx, dx)
y_grid = arange(dy/2, Ly, dy)

#make xx_grid, yy_grid, which are n by n matrix with the repetition of x_grid, y_grid,
#respectively
#they are the x, y coordinates of the particles
xx_grid, yy_grid = meshgrid(x_grid, y_grid)

#flatten the coordinate values
x_initial = xx_grid.flatten()
y_initial = yy_grid.flatten()

#print("initial x",x_initial, "initial y",y_initial)
#print(xx_grid)

#distance rij
def dis(r):
    return ((r[0]-r[1])**2+(r[2]-r[3])**2)**(1/2)


#Define Lennard-Jones potential,using the given equation
def p(r):
    sig=1
    ep=1
    return 4*ep*((sig/dis(r))**12-(sig/dis(r))**6)

#Define F, which is the gradient of Potential,
#in terms of the vector coordinate r consisting of the corordinates of both two particles
#F given by taking the gradient of V    
m=1
e=1
o=1
def F(r,t):
    x=r[0]
    x1=r[1]
    y=r[2]
    y1=r[3]
#    vx=r[4]
#    vx1=r[5]
#    vy=r[6]
#    vy1=r[7]
#    fx=vx
#    fx1=vx1
#    fy=vy
#    fy1=vy1
    fvx=(4*e/m)*(o**12*(-6)*((x1-x)**2+(y1-y)**2)**(-7)*2*(x1-x)*(-1)+o**6*3*((x1-x)**2+(y1-y)**2)**(-4)*2*(x1-x)*(-1))
    fvx1=(4*e/m)*(o**12*(-6)*((x1-x)**2+(y1-y)**2)**(-7)*2*(x1-x)+o**6*3*((x1-x)**2+(y1-y)**2)**(-4)*2*(x1-x))
    fvy=(4*e/m)*(o**12*(-6)*((x1-x)**2+(y1-y)**2)**(-7)*2*(y1-y)*(-1)+o**6*3*((x1-x)**2+(y1-y)**2)**(-4)*2*(y1-y)*(-1))
    fvy1=(4*e/m)*(o**12*(-6)*((x1-x)**2+(y1-y)**2)**(-7)*2*(y1-y)+o**6*3*((x1-x)**2+(y1-y)**2)**(-4)*2*(y1-y))
    
    
    fvx=-fvx
    fvx1=-fvx1
    fvy=-fvy
    fvy1=-fvy1
    
    return array([fvx,fvx1,fvy,fvy1],float)

def abF(r):
    x=r[0]
    x1=r[1]
    y=r[2]
    y1=r[3]
#    vx=r[4]
#    vx1=r[5]
#    vy=r[6]
#    vy1=r[7]
#    fx=vx
#    fx1=vx1
#    fy=vy
#    fy1=vy1
    fvx=(4*e/m)*(o**12*(-6)*((x1-x)**2+(y1-y)**2)**(-7)*2*(x1-x)*(-1)+o**6*3*((x1-x)**2+(y1-y)**2)**(-4)*2*(x1-x)*(-1))
    fvx1=(4*e/m)*(o**12*(-6)*((x1-x)**2+(y1-y)**2)**(-7)*2*(x1-x)+o**6*3*((x1-x)**2+(y1-y)**2)**(-4)*2*(x1-x))
    fvy=(4*e/m)*(o**12*(-6)*((x1-x)**2+(y1-y)**2)**(-7)*2*(y1-y)*(-1)+o**6*3*((x1-x)**2+(y1-y)**2)**(-4)*2*(y1-y)*(-1))
    fvy1=(4*e/m)*(o**12*(-6)*((x1-x)**2+(y1-y)**2)**(-7)*2*(y1-y)+o**6*3*((x1-x)**2+(y1-y)**2)**(-4)*2*(y1-y))
    
    
    fvx=-fvx
    fvx1=-fvx1
    fvy=-fvy
    fvy1=-fvy1
    
    return (fvx**2+fvy**2)**(1/2)
    

#g= dv/dt = F/m
def g(r,t):
    m=1
    return F(r,t)/m


    

#Define Verlet Method
def vm(r):
    #set initial values
    a=0.0
    b=10
    N=1000
    h=(b-a)/N
    tpoints=arange(a,b,h)
    xpoints=[]
    x1points=[]
    ypoints=[]
    y1points=[]


#i).i. r1 = [4, 4],r2 = [5.2, 4], 
    v=zeros(4)

    v_plus_half_h=v+0.5*h*  F(r,a)



#appending the new x, y values to the lists
    for t in tpoints:
        xpoints.append(r[0])
        x1points.append(r[1])
        ypoints.append(r[2])
        y1points.append(r[3])
#    vxpoints.append(r[4])
#    vx1points.append(r[5])
#    vypoints.append(r[6])
#    vy1points.append(r[7])
   
        r += h*v_plus_half_h
        k=h*F(r,t+h)
#    v_plus_h = v_plus_half_h+0.5*k
        v_plus_3_over_2_h=v_plus_half_h+k
        v_plus_half_h=v_plus_3_over_2_h
        
        
    k=array([xpoints,x1points,ypoints,y1points])
    return k


k=vm([0.5,2.5,0.5,0.5])
#g=vm([0.5,3.5,0.5,0.5])
#print(k)
#print(g[2])




  


#kinetic energy for one particles w.r.p. to another through time
def ke(r):
    #set initial values
    a=0.0
    b=10
    N=1000
    h=(b-a)/N
    tpoints=arange(a,b,h)
    xpoints=[]
    x1points=[]
    ypoints=[]
    y1points=[]
    ke=[]


#i).i. r1 = [4, 4],r2 = [5.2, 4], 
    v=zeros(4)

    v_plus_half_h=v+0.5*h*  F(r,a)
    



#appending the new x, y values to the lists
    for t in tpoints:
        ke.append((v[0]**2+v[2]**2)/2)
        #xpoints.append(r[0])
        #x1points.append(r[1])
        #ypoints.append(r[2])
        #y1points.append(r[3])
#    vxpoints.append(r[4])
#    vx1points.append(r[5])
#    vypoints.append(r[6])
#    vy1points.append(r[7])
   
        r += h*v_plus_half_h
        k=h*F(r,t+h)
        v = v_plus_half_h+0.5*k
        v_plus_3_over_2_h=v_plus_half_h+k
        v_plus_half_h=v_plus_3_over_2_h
        
    
    
    return ke


 





   
    
#distance rij
def dis(r):
    return ((r[0]-r[1])**2+(r[2]-r[3])**2)**(1/2)



        
#acceleration
def ac(r,t):

    return F(r,t)/m

x=x_initial 
y=y_initial 


#KE of whole sys
KE=zeros(1000)
r=zeros(4)
for i in range(16):
    for u in range(16):
        if u==i:
            pass
        else:
            r[0]=x[i]
            r[1]=x[u]
            r[2]=y[i]
            r[3]=y[u]
            KE+=ke(r)

#print(KE)


#Potential Energy of one by another thorugh time
def pe(r):
    #set initial values
    a=0.0
    b=10
    N=1000
    h=(b-a)/N
    tpoints=arange(a,b,h)

    pe=[]


#i).i. r1 = [4, 4],r2 = [5.2, 4], 
    v=zeros(4)

    v_plus_half_h=v+0.5*h*  F(r,a)
    



#appending the new x, y values to the lists
    for t in tpoints:
        pe.append(p(r)/2)
        #xpoints.append(r[0])
        #x1points.append(r[1])
        #ypoints.append(r[2])
        #y1points.append(r[3])
#    vxpoints.append(r[4])
#    vx1points.append(r[5])
#    vypoints.append(r[6])
#    vy1points.append(r[7])
   
        r += h*v_plus_half_h
        k=h*F(r,t+h)
        v = v_plus_half_h+0.5*k
        v_plus_3_over_2_h=v_plus_half_h+k
        v_plus_half_h=v_plus_3_over_2_h
        
    
    
    return pe

"""




"""


#PE
PE=zeros(1000)
r=zeros(4)
for i in range(16):
    for u in range(16):
        if u==i:
            pass
        else:
            r[0]=x[i]
            r[1]=x[u]
            r[2]=y[i]
            r[3]=y[u]
            PE+=pe(r)




E=PE+KE


tpoints=arange(0,10,0.01)                
        
plot(E,label="Energy Through time")
xlabel("time(s)")
ylabel("E(J)")
title("Energy of system vs. time")

show()






#plotting
#assigning the ini x and y to two l=elements of r 
#loop through len(x)
totx=[]
toty=[]
a=0.0
b=10
N=1000
t=arange(a,b,0.01)
pe=0
#allpe=zeros(1000)
for u in range(len(x)):
    

    #zero array of xpoints and y points
    xpoints=zeros(1000)
    ypoints=zeros(1000)
    
    #loop throough the elements except the above chosen one
    
    for p in range(16):
        
        r=zeros(4)
        r[0]=x[u]
        r[2]=y[u]
        
        if p==u:
            pass
        else:
            r[1]=x[p]
            r[3]=y[p]
            #print(r)
            a=vm(r)
        
            #pe=pe-dis(r)*abF(r)
        
            #print(a[0])
      
         #add the time evolving x y points given by each set of r together 
            
            
        
            
            xpoints=xpoints+a[0]/15
            ypoints=ypoints+a[2]/15
            
            #print("xpoints",xpoints)
            #print("ypoints",ypoints)
            
            
            
            
    
    #xpoints=xpoints/15
    #ypoints=ypoints/15  
    totx.append(xpoints)
    toty.append(ypoints)
            
    #print("xpoints(unit)",xpoints)
    #print("ypoints(unit)",ypoints)
    
    
    plot(xpoints,ypoints,label="{} th particle".format(u+1))          

      
xlabel("x direction(units)")
ylabel("y direction(units)")
title("orbit trajactory of 16 particles")

legend(loc='center left', bbox_to_anchor=(1, 0.5))
show()




