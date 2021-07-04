#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import math

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




class particle:
    def __init__(self):
        self.mass= 1
        self.radius = 1
        self.pos = [0,0,0]
        self.vel = [0,0,0]
        self.acc =[0,0,0]
        #self.PE = 0
        self.KE = 0.5*self.mass*(self.vel[0]**2+self.vel[1]**2+self.vel[2]**2)
        #self.TE = 0

    def data(self):
        self.mass = input('enter mass : ')
        self.radius = input('enter radius : ')
        self.pos = input('enter initial position vector : ')
        self.vel = input('enter initial velocity vector: ')
        self.KE = 0.5*self.mass*(self.vel[0]**2+self.vel[1]**2+self.vel[2]**2)
        #self.acc = acc







systm = []

G = 6.674e-11
num = input('enter number of particles: ')
for i in range(num):
    ins = particle()
    ins.data()
    systm.append(ins)

print systm[1].mass, systm[1].pos, systm[1].vel

endtime = input ('enter time for sim end : ')
timeres = input ('enter time res : ')

timesteps = endtime/timeres
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

file1 = open('partdata.txt','w')
for  i in range(timesteps):
    timeplot.append(i*timeres)
    print 'time: '+ str(i*timeres)
    file1.write('time: '+ str(i*timeres)+'\n')
    PE = 0
    TKE = 0
    TE = 0
    for i in range(num):
        xplot[i].append(systm[i].pos[0])
        yplot[i].append(systm[i].pos[1])
        zplot[i].append(systm[i].pos[2])

    for part  in range(num):

        TKE += systm[part].KE
        for accord in range(3):
            #print accord
            for ot_part in range(num):
                #print ot_part
                if(ot_part != part):
                    r = math.sqrt((systm[part].pos[0] - systm[ot_part].pos[0])**2 + (systm[part].pos[1] - systm[ot_part].pos[1])**2 +(systm[part].pos[2] - systm[ot_part].pos[2])**2 )
                    #print r


                    systm[part].acc[accord] = -(systm[part].acc[accord] + (G*systm[ot_part].mass*(systm[part].pos[accord] - systm[ot_part].pos[accord])/(r)**3))
                    #print r,part, ot_part, systm[part].acc[accord],systm[ot_part].acc[accord],systm[part].pos[accord] - systm[ot_part].pos[accord]
                    PE -= G*(systm[part].mass*systm[ot_part].mass)/2*r




    #q_for_col=[]
    TE = TKE+PE
    for part  in range(num):

        for ot_part in range(num):
                #print ot_part
            if(ot_part != part):
                r = math.sqrt((systm[part].pos[0] - systm[ot_part].pos[0])**2 + (systm[part].pos[1] - systm[ot_part].pos[1])**2 +(systm[part].pos[2] - systm[ot_part].pos[2])**2 )
                    #print r
                if r <= systm[part].radius+systm[ot_part].radius  and col_matrix[part][ot_part] ==False:
                    print('collision')
                    file1.write('collision\n' + str(part)+'\t'+str(ot_part)+'\n')
                    temp_col =collision(systm[part],systm[ot_part])
                    systm[part] = temp_col[0]
                    systm[ot_part] = temp_col[1]
                    #q_for_col.append((part,ot_part))
                    col_matrix[part][ot_part] = True
                    col_matrix[ot_part][part] = True

                if col_matrix[part][ot_part] ==True and r > systm[part].radius+systm[ot_part].radius:
                    col_matrix[part][ot_part] =False
                    col_matrix[ot_part][part] = False
                   # sys.exit(0)
                        #systm[part] = coll_arr[0]
                        #systm[ot_part] = coll_arr[1]

        print  'pos of '+ str(part) +' = ',systm[part].pos, 'vel'+str(part)+'= ',systm[part].vel  ,' acc'+str(part)+'= ',systm[part].acc
        file1.write('pos of '+ str(part) +' = '+str(systm[part].pos)+ 'vel'+str(part)+'= '+str(systm[part].vel)  +' acc'+str(part)+'= '+str(systm[part].acc)+'\t'+str(systm[part].KE)+'\n')



    file1.write('TKE: '+str(TKE)+'\t'+'PE: '+str(PE)+'\t'+'TE: '+str(TE)+'\n')

    for par in range(num):
        #print systm[par].acc
        for coord in range(3):
            systm[par].pos[coord] = systm[par].pos[coord] + (systm[par].vel[coord]*timeres)
            systm[par].vel[coord] = systm[par].vel[coord] + (systm[par].acc[coord]*timeres)
            systm[par].acc[coord] = 0
            systm[par].KE = 0.5*systm[par].mass*(systm[par].vel[0]**2+systm[par].vel[1]**2+systm[par].vel[2]**2)
            #systm[par].TE =

file1.close()


fig = plt.figure()
ax = Axes3D(fig)
for i in range(num):

    ax.plot(xplot[i], yplot[i], zplot[i])
    ax.legend()

    print systm[i]
    #time.sleep(10)
    #plt.close()

plt.show(fig)
