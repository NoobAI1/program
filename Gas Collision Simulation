from vpython import *
import random
import sys
import time
from itertools import product

start_time = time.time()

n = 1000
T = 300
dt = 0.00001
m = 4E-3/6E23
kb = 1.38E-23
vrms = ((3*kb*T)/(m))**(1/2)

scene = canvas(title="Ideal Gas Collision Simulation", width=600, height=600, range=1)

#scene
d = 1
container = curve(color=color.white, radius=0.01)
side1 = curve(color=color.white, radius=0.01)
side2 = curve(color=color.white, radius=0.01)
side3 = curve(color=color.white, radius=0.01)
container.append([vector(d,0,0), vector(0,d,0), vector(-d,0,0), vector(0,-d,0), vector(d,0,0)])
container.append([vector(d,0,2**(1/2)*d), vector(0,d,2**(1/2)*d), vector(-d,0,2**(1/2)*d), vector(0,-d,2**(1/2)*d), vector(d,0,2**(1/2)*d)])
side1.append([vector(0,d,0), vector(0,d,2**(1/2)*d)])
side2.append([vector(-d,0,0), vector(-d,0,2**(1/2)*d)])
side3.append([vector(0,-d,0), vector(0,-d,2**(1/2)*d)])

container.rotate(axis=vec(0,0,-1), angle=pi/4, origin=vec(0,0,0))
side1.rotate(axis=vec(0,0,-1), angle=pi/4, origin=vec(0,0,0))
side2.rotate(axis=vec(0,0,-1), angle=pi/4, origin=vec(0,0,0))
side3.rotate(axis=vec(0,0,-1), angle=pi/4, origin=vec(0,0,0))

#spheres
spheres = []
r = 0.02
for i in range(n):
    x = random.uniform(2**(-1/2)*d-r, -(2**(-1/2)*d-r))
    y = random.uniform(2**(-1/2)*d-r, -(2**(-1/2)*d-r))
    z = random.uniform(2**(1/2)*d-r, 0+r)
    spheres.append(sphere(pos=vector(x,y,z), color=color.red, radius=r))

velocities = []
for i in range(n):
    vx = vrms
    vy = 0
    vz = 0
    velocities.append((vector(vx,vy,vz)))

def col_v(v1, v2, r1, r2):
    v1_a = v1 + dot((v2-v1), (r1-r2))/mag2(r1-r2)*(r1-r2)
    v2_a = v2 + dot((v1-v2), (r2-r1))/mag2(r2-r1)*(r2-r1)
    return (v1_a, v2_a)

def col():
    grid = {}
    grid_size = 4 * r
    
    for idx in range(n):
        cell = (
            int(spheres[idx].pos.x / grid_size),
            int(spheres[idx].pos.y / grid_size),
            int(spheres[idx].pos.z / grid_size)
        )
        if cell not in grid:
            grid[cell] = []
        grid[cell].append(idx)
    
    pairs = []
    offsets = list(product([-1, 0, 1], repeat=3))  
    for cell in grid:
        for dx, dy, dz in offsets:
            neighbor_cell = (cell[0]+dx, cell[1]+dy, cell[2]+dz)
            if neighbor_cell in grid:  
                for i in grid[cell]:
                    for j in grid[neighbor_cell]:
                        if i < j:  
                            if mag2(spheres[i].pos - spheres[j].pos) <= (2*r)**2:
                                pairs.append((i, j))
    return pairs
           
def col_wall_v():
    for i in spheres:
        if i.pos.x >= 2**(-1/2)*d-r or i.pos.x <= -(2**(-1/2)*d-r):
            velocities[spheres.index(i)].x = -velocities[spheres.index(i)].x
        if i.pos.y >= 2**(-1/2)*d-r or i.pos.y <= -(2**(-1/2)*d-r):
            velocities[spheres.index(i)].y = -velocities[spheres.index(i)].y
        if i.pos.z >= 2**(1/2)*d-r or i.pos.z <= 0+r:
            velocities[spheres.index(i)].z = -velocities[spheres.index(i)].z   

g1 = graph(title='Rate Distribution',ytitle='N',xtitle='v(m/s)', xmax = 5000, ymax = n/20)
gc = gcurve()
dv = 50
gc.delete()
for v in arange(0, 5000, dv):
    gc.plot(v, ((2/pi)*((m/(kb*T))**3))**0.5*(v**2)*exp(((-m)*(v**2))/(2*kb*T))*dv*n)

gb = gvbars(delta=dv)
bars = []
for i in range(10000//dv):
    if i == int(vrms/dv):
        bars.append([dv*(i+0.5), n])
    else:
        bars.append([dv*(i+0.5), 0])

while True:
    rate(300)
    end_time = time.time()
    exe_time = end_time - start_time
    if exe_time > 60:
        sys.exit()

    for i in range(n):
        spheres[i].pos += velocities[i]*dt

    col_list = col()
    if col_list != []:
        for i in range(len(col_list)):
            v1 = velocities[col_list[i][0]]
            v2 = velocities[col_list[i][1]]
            v1_a, v2_a= col_v(v1, v2, spheres[col_list[i][0]].pos, spheres[col_list[i][1]].pos)
            if int(mag(v1)/dv) < len(bars) and int(mag(v1_a)/dv) < len(bars):
                bars[int(mag(v1)/dv)][1] -= 1
                bars[int(mag(v1_a)/dv)][1] += 1
            if int(mag(v2)/dv) < len(bars) and int(mag(v2_a)/dv) < len(bars):    
                bars[int(mag(v2)/dv)][1] -= 1
                bars[int(mag(v2_a)/dv)][1] += 1
            velocities[col_list[i][0]] = v1_a
            velocities[col_list[i][1]] = v2_a
        col_list.clear()

    col_wall_v()

    gb.data = bars
