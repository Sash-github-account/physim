#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
#import time
import numpy as np
import math
#from constraint import *
#from physim.py import collision

E0= 8.854e-12
pi= 3.14
Mu0 = 4e-7*pi
G = 6.674e-11

#def RK_algo():

#NOTE: inaccurate algorithm (error prop. to timeres)


def lagrange_interpole(x, f_x, x_comp):

    P_x = 0
    for k in range(len(x)):
        Lambda = 1
        for l in range(len(x)):
            if  k!= l:
                Lambda = Lambda*(x_comp - x[l])/(x[k] - x[l])
        P_x += Lambda*f_x[k]
    return P_x

#def Neville_algo(x, f_x, x_comp):


def cross(a,b):
    p = []
    for i in range(len(a)):
        if i == 0:
            p.append(a[1]*b[2]-a[2]*b[1])
        elif i ==1:
            p.append(a[2]*b[0]-b[2]*a[0])
        elif i== 2:
            p.append(a[0]*b[1]-b[0]*a[1])
    return p



def grav_poten(sys,r):
    GP =0
    for i in range(len(sys)):
        d = math.sqrt((sys[i].pos[0] - r[0])**2+(sys[i].pos[1] - r[1])**2+(sys[i].pos[2] - r[2])**2)
        GP += (G*sys[i].mass/d)
    return GP



def stat_elec_poten(sys,r):
    EP =0
    for i in range(len(sys)):
        d = math.sqrt((sys[i].pos[0] - r[0])**2+(sys[i].pos[1] - r[1])**2+(sys[i].pos[2] - r[2])**2)
        EP += ((1/4*pi*E0)*sys[i].charge/d)
    return EP




def grav_field(sys,r):
    GF = [0,0,0]
    for i in range(len(sys)):
        d = math.sqrt((sys[i].pos[0] - r[0])**2+(sys[i].pos[1] - r[1])**2+(sys[i].pos[2] - r[2])**2)
        if d!=0:
            d_vec = [sys[i].pos[0] -r[0], sys[i].pos[1]-r[1],sys[i].pos[2]-r[2]]
            for j in range(3):
                GF[j] += G* sys[i].mass*d_vec[j]/d**3
    return GF


def mag_field(sys,r):
    B=[0,0,0]
    for i in range(len(sys)):
        d = math.sqrt((sys[i].pos[0] - r[0])**2+(sys[i].pos[1] - r[1])**2+(sys[i].pos[2] - r[2])**2)
        if d != 0:
            d_vec = [ r[0]-sys[i].pos[0] ,r[1]-sys[i].pos[1],r[2]-sys[i].pos[2]]
            v_x_d =cross(sys[i].vel,d_vec)
            for j in range(3):
                B[j] += (Mu0/4/pi)*sys[i].charge*v_x_d[j]/d**3
    return B



def elec_field(sys,r):
    E=[0,0,0]
    for i in range(len(sys)):
        d = math.sqrt((sys[i].pos[0] - r[0])**2+(sys[i].pos[1] - r[1])**2+(sys[i].pos[2] - r[2])**2)
        if d !=0:
            d_vec = [ r[0]-sys[i].pos[0] ,r[1]-sys[i].pos[1],r[2]-sys[i].pos[2]]
            for j in range(3):
                E[j]+=(1/(4*pi*E0))*sys[i].charge*d_vec[j]/d**3
    return E


def collision(a,b):
    anga = math.acos((((b.pos[0] - a.pos[0])*a.vel[0])+((b.pos[1] - a.pos[1])*a.vel[1])+((b.pos[2] - a.pos[2])*a.vel[2]))/math.sqrt(((b.pos[0] - a.pos[0])**2+(b.pos[1] - a.pos[1])**2+(b.pos[2] - a.pos[2])**2)*((a.vel[0])**2+(a.vel[1])**2+(a.vel[2])**2)))
    angb = math.acos((((a.pos[0] - b.pos[0])*b.vel[0])+((a.pos[1] - b.pos[1])*b.vel[1])+((a.pos[2] - b.pos[2])*b.vel[2]))/math.sqrt(((a.pos[0] - b.pos[0])**2+(a.pos[1] - b.pos[1])**2+(a.pos[2] - b.pos[2])**2)*((b.vel[0])**2+(b.vel[1])**2+(b.vel[2])**2)))

    file1.write(str( anga)+'\t'+str(angb) +'\t'+str(math.cos(anga))+'\t'+str(math.cos(angb))+'\n')

    #angaxz = math.acos((((b.pos[0] - a.pos[0])*a.vel[0])+((b.pos[2] - a.pos[2])*a.vel[2]))/math.sqrt(((b.pos[0] - a.pos[0])**2+(b.pos[2] - a.pos[2])**2)*((a.vel[0])**2+(a.vel[2])**2)))
    #angaxy = math.acos((((b.pos[0] - a.pos[0])*a.vel[0])+((b.pos[1] - a.pos[1])*a.vel[1]))/math.sqrt(((b.pos[0] - a.pos[0])**2+(b.pos[1] - a.pos[1])**2)*((a.vel[0])**2+(a.vel[1])**2)))
    #angbxz = math.acos((((a.pos[0] - b.pos[0])*b.vel[0])+((a.pos[2] - b.pos[2])*b.vel[2]))/math.sqrt(((a.pos[0] - b.pos[0])**2+(a.pos[2] - b.pos[2])**2)*((b.vel[0])**2+(b.vel[2])**2)))
    #angbxy = math.acos((((a.pos[0] - b.pos[0])*b.vel[0])+((a.pos[1] - b.pos[1])*b.vel[1]))/math.sqrt(((a.pos[0] - b.pos[0])**2+(a.pos[1] - b.pos[1])**2)*((b.vel[0])**2+(b.vel[1])**2)))

    vcena = math.sqrt(a.vel[0]**2+a.vel[1]**2+a.vel[2]**2)*math.cos(anga)
    vcenb = math.sqrt(b.vel[0]**2+b.vel[1]**2+b.vel[2]**2)*math.cos(angb)

    vcenax = vcena*(b.pos[0]-a.pos[0])/math.sqrt((b.pos[0]-a.pos[0])**2+(b.pos[1]-a.pos[1])**2+(b.pos[2]-a.pos[2])**2)
    vcenay = vcena*(b.pos[1]-a.pos[1])/math.sqrt((b.pos[0]-a.pos[0])**2+(b.pos[1]-a.pos[1])**2+(b.pos[2]-a.pos[2])**2)
    vcenaz = vcena*(b.pos[2]-a.pos[2])/math.sqrt((b.pos[0]-a.pos[0])**2+(b.pos[1]-a.pos[1])**2+(b.pos[2]-a.pos[2])**2)

    vcenbx = vcenb*(a.pos[0]-b.pos[0])/math.sqrt((b.pos[0]-a.pos[0])**2+(b.pos[1]-a.pos[1])**2+(b.pos[2]-a.pos[2])**2)
    vcenby = vcenb*(a.pos[1]-b.pos[1])/math.sqrt((b.pos[0]-a.pos[0])**2+(b.pos[1]-a.pos[1])**2+(b.pos[2]-a.pos[2])**2)
    vcenbz = vcenb*(a.pos[2]-b.pos[2])/math.sqrt((b.pos[0]-a.pos[0])**2+(b.pos[1]-a.pos[1])**2+(b.pos[2]-a.pos[2])**2)

    file1.write(str( vcena)+'\t'+str( vcenax)+'\t'+str(vcenay)+'\t'+str(vcenaz)+'\n')
    file1.write(str( vcenb)+'\t'+str( vcenbx)+'\t'+str(vcenby)+'\t'+str(vcenbz)+'\n')


    vnormax = a.vel[0]-vcenax
    vnormay = a.vel[1]-vcenay
    vnormaz = a.vel[2]-vcenaz

    vnormbx = b.vel[0]-vcenbx
    vnormby = b.vel[1]-vcenby
    vnormbz = b.vel[2]-vcenbz

    vnorma = math.sqrt(vnormax**2+vnormay**2+vnormaz**2)
    vnormb = math.sqrt(vnormbz**2+vnormby**2+vnormbz**2)

    file1.write(str( vnorma)+'\t'+str( vnormax)+'\t'+str(vnormay)+'\t'+str(vnormaz)+'\n')
    file1.write(str( vnormb)+'\t'+str( vnormbx)+'\t'+str(vnormby)+'\t'+str(vnormbz)+'\n')


    tempx = vcenax
    vcenax = vcenbx*b.mass/a.mass
    vcenbx = tempx*a.mass/b.mass

    tempy = vcenay
    vcenay = vcenby*b.mass/a.mass
    vcenby = tempy*a.mass/b.mass

    tempz = vcenaz
    vcenaz = vcenbz*b.mass/a.mass
    vcenbz = tempz*a.mass/b.mass

    file1.write(str( vcena)+'\t'+str( vcenax)+'\t'+str(vcenay)+'\t'+str(vcenaz)+'\n')
    file1.write(str( vcenb)+'\t'+str( vcenbx)+'\t'+str(vcenby)+'\t'+str(vcenbz)+'\n')


    a.vel[0] = vnormax +vcenax
    a.vel[1] = vnormay +vcenay
    a.vel[2] = vnormaz +vcenaz

    file1.write(str( a.vel)+'\n')
    b.vel[0] = vnormbx +vcenbx
    b.vel[1] = vnormby +vcenby
    b.vel[2] = vnormbz +vcenbz
    file1.write( str(b.vel)+'\n')

    return [a,b]



class Emag_particle:
      def __init__(self):
          self.mass= 1
          self.radius =1
          self.pos = [0,0,0]
          self.vel = [0,0,0]
          self.acc = [0,0,0]
          self.KE = 0.5*self.mass*(self.vel[0]**2+self.vel[1]**2+self.vel[2]**2)
          self.charge = 1

      def data(self):
          self.mass =input('enter mass: ')
          self.radius = input('enter radius: ')
          self.charge = input('enter charge: ')
          self.pos = input('enter initial position vector: ')
          self.vel = input('enter initial velocity vector: ')
          self.KE = 0.5*self.mass*(self.vel[0]**2+self.vel[1]**2+self.vel[2]**2)

class obj_rotation(Emag_particle):
    def __init__(self):
        Emag_particle.__init__(self)
        #self.number_part = 1
        self.pos_cm =[0,0,0]
        self.dist_frm_centre = [0,0,0]
        self.ang_vel = [0,0,0]
        self.torque = [0,0,0]

        def plot_object():

            fig = plt.figure()
            ax = p3.Axes3D(fig)
            #line, = ax.plot([],[],[])
            ax.set_xlim(min([min(xplot[i]) for i in range(num)]),max([max(xplot[i]) for i in range(num)]))
            ax.set_xlabel('X')
            ax.set_ylim(min([min(yplot[i]) for i in range(num)]),max([max(yplot[i]) for i in range(num)]))
            ax.set_ylabel('Y')
            ax.set_zlim(min([min(zplot[i]) for i in range(num)]),max([max(zplot[i]) for i in range(num)]))
            ax.set_zlabel('Z')

            for i in range(num):

                ax.plot(xplot[i], yplot[i], zplot[i])
                ax1.plot(xplot[i], yplot[i], zplot[i])
                ax.legend()

                print systm[i]
                #time.sleep(10)
                #plt.close()

            plt.show(fig)







systm = []
file2 = open('input_set.txt','r')
#num = input('enter number of particles: ')
file2_lines = file2.readlines()
num = int(file2_lines[0])
print num
#read_choice =int(input('1.enter manually 2.read from input file:'))
read_choice = int(file2_lines[1])
endtime = float(file2_lines[3])
timeres = float(file2_lines[4])
throttle_ip = int(file2_lines[5])
print read_choice
if read_choice ==2:
    #file_name = str(input('enter input data file name: '))
    file_name = str(file2_lines[2].rstrip())
    print file_name
    file1= open(file_name,'r')
    read_lines = file1.readlines()
    data = [[read_lines[i+(j*5)] for i in range(5)] for j in range(num)]
    print data
    file2.close()
    file1.close()

for i in range(num):
    if(read_choice == 1):
        ins = Emag_particle()
        ins.data()
        systm.append(ins)
    elif read_choice ==2:
        ins = Emag_particle()
        ins.mass = float(data[i][0])
        #print ins.mass
        ins.radius = float(data[i][1])
        ins.charge = float(data[i][2])
        ins.pos = data[i][3].split(',')
        ins.pos = [float(ins.pos[j]) for j in range(3)]
        ins.vel = data[i][4].split(',')
        ins.vel = [float(ins.vel[j]) for j in range(3)]
        ins.KE = 0.5*ins.mass*(ins.vel[0]**2+ins.vel[1]**2+ins.vel[2]**2)
        systm.append(ins)



#print systm[1].mass, systm[1].pos, systm[1].vel


timesteps = int(endtime/timeres)
print timesteps
xplot =[]
yplot =[]
zplot =[]

for i in range(num):
    xplot.append([])
    yplot.append([])
    zplot.append([])

timeplot = []

col_matrix = []
col_matrix=[[False for i in range(num)] for j in range(num)]

thrle_cnt =0
file1 = open('partdata.txt','w')
for  i in range(timesteps):
    timeplot.append(i*timeres)
    print 'time: '+ str(i*timeres)
    file1.write('time: '+ str(i*timeres)+'\n')
    PE = 0
    TKE = 0
    TE = 0

    if  thrle_cnt == throttle_ip :
        thrle_cnt = 0
        for i in range(num):
            xplot[i].append(systm[i].pos[0])
            yplot[i].append(systm[i].pos[1])
            zplot[i].append(systm[i].pos[2])
    thrle_cnt += 1




file1.close()


x=[2.3,2.7,2.9,3.2,3.5,3.7]
f_x =[6.38,13.6,18.7,28.26,40.4, 50]
P = lagrange_interpole(x,f_x,3.0)
print P

fig = plt.figure()
ax = p3.Axes3D(fig)
#line, = ax.plot([],[],[])
ax.set_xlim(min([min(xplot[i]) for i in range(num)]),max([max(xplot[i]) for i in range(num)]))
ax.set_xlabel('X')
ax.set_ylim(min([min(yplot[i]) for i in range(num)]),max([max(yplot[i]) for i in range(num)]))
ax.set_ylabel('Y')
ax.set_zlim(min([min(zplot[i]) for i in range(num)]),max([max(zplot[i]) for i in range(num)]))
ax.set_zlabel('Z')

for i in range(num):
    fig1 = plt.figure()
    #ax = Axes3D(fig)
    ax1 = p3.Axes3D(fig1)
    ax1.set_xlim(min(xplot[i]),max(xplot[i]))
    ax1.set_xlabel('X')
    ax1.set_ylim(min(yplot[i]),max(yplot[i]))
    ax1.set_ylabel('Y')
    ax1.set_zlim(min(zplot[i]),max(zplot[i]))
    ax1.set_zlabel('Z')

    ax.plot(xplot[i], yplot[i], zplot[i])
    ax1.plot(xplot[i], yplot[i], zplot[i])
    ax.legend()

    print systm[i]
    #time.sleep(10)
    #plt.close()

plt.show(fig1)
plt.show(fig)





fig3 = plt.figure()
ax3 = p3.Axes3D(fig3)
colors = plt.cm.jet(np.linspace(0, 1, num))

# set up lines and points
lines = sum([ax3.plot([], [], [], '-', c=c)
             for c in colors], [])
pts = sum([ax3.plot([], [], [], 'o', c=c)
           for c in colors], [])

ax3.set_xlim(min([min(xplot[i]) for i in range(num)]),max([max(xplot[i]) for i in range(num)]))
ax3.set_xlabel('X')
ax3.set_ylim(min([min(yplot[i]) for i in range(num)]),max([max(yplot[i]) for i in range(num)]))
ax3.set_ylabel('Y')
ax3.set_zlim(min([min(zplot[i]) for i in range(num)]),max([max(zplot[i]) for i in range(num)]))
ax3.set_zlabel('Z')


def init():
    for line, pt in zip(lines, pts):
        line.set_data([], [])
        line.set_3d_properties([])

        pt.set_data([], [])
        pt.set_3d_properties([])
    return lines + pts

def animate(i):

    for line, pt, xi, yi, zi in zip(lines, pts, xplot, yplot, zplot):
        print 'entered'
        x, y, z = xi[:i], yi[:i], zi[:i]
        line.set_data(x, y)
        line.set_3d_properties(z)

        pt.set_data(x[-1:], y[-1:])
        pt.set_3d_properties(z[-1:])

    #ax.view_init(30, 0.3 * i)
    #fig.canvas.draw()
    return lines + pts


# animation function.  This will be called sequentially with the frame number
anim = animation.FuncAnimation(fig3, animate, init_func=init, frames=(len(range(timesteps))/throttle_ip), interval=1, blit=False)

plt.show()
