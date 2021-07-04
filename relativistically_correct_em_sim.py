#!/usr/bin/env python
#---------#
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import time
import math
import numpy as np
from IPython.display import HTML
import os, psutil
#---------#



#---------#
E0= 8.854e-12
pi= np.pi
Mu0 = 4e-7*pi
G = 6.674e-11
C_EM = 3e8
#---------#

#def RK_algo():

#---------#
def lorentz_factor(particle):
    vel_mag = math.sqrt((particle.vel[0])**2 + (particle.vel[1])**2 + (particle.vel[2])**2)
    print("mag vel: " + str(vel_mag))
    gmma = 1/(math.sqrt(1-(vel_mag**2/C_EM**2)))
    return gmma
#---------#



#---------#
def constrt_4_vel(particle):
    gmma = lorentz_factor(particle)
    print("lorentz factor: " + str(gmma))
    four_vel = []
    for i in range(4):
       if(i == 0): four_vel.append(gmma*C_EM)
       else: four_vel.append(gmma*particle.vel[i-1])
    print("four vel: " + str(four_vel))
    return np.array(four_vel)
#---------#




#---------#
def constrt_EM_tensor(E, B, particle):
    row1 = [0,            (E[0]/C_EM), (E[1]/C_EM), (E[2]/C_EM)]
    row2 = [-(E[0]/C_EM),           0,        B[2],       -B[1]]
    row3 = [-(E[1]/C_EM),       -B[2],           0,        B[0]]
    row4 = [-(E[2]/C_EM),        B[1],       -B[0],           0]
    F_tensor = np.array([row1, row2, row3, row4])
    #print(" field tensor: " + str(F_tensor))
    return F_tensor
#---------#



#---------#
def minkowsky_force(F_tensor, four_vel, particle):
    mink_force = [0,0,0,0]
    for i in range(4):
        for j in range(4):
            mink_force[i] = mink_force[i] + (four_vel[i] * F_tensor[j][i])
    #mink_force = particle.charge * np.matmul(four_vel, F_tensor)
    #print("Minkowsky force: " + str(mink_force))
    return mink_force
#---------#



#---------#
def lagrange_interpole(x, f_x, x_comp):

    P_x = 0
    for k in range(len(x)):
        Lambda = 1
        for l in range(len(x)):
            if  k!= l:
                Lambda = Lambda*(x_comp - x[l])/(x[k] - x[l])
        P_x += Lambda*f_x[k]
    return P_x
#---------#

#def Neville_algo(x, f_x, x_comp):

#---------#
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
#---------#



#---------#
def grav_poten(sys,r):
    GP =0
    for i in range(len(sys)):
        d = math.sqrt((sys[i].pos[0] - r[0])**2+(sys[i].pos[1] - r[1])**2+(sys[i].pos[2] - r[2])**2)
        if(d!=0): GP += (G*sys[i].mass/d)
    return GP
#---------#



#---------#
def stat_elec_poten(sys,r):
    EP =0
    for i in range(len(sys)):
        d = math.sqrt((sys[i].pos[0] - r[0])**2+(sys[i].pos[1] - r[1])**2+(sys[i].pos[2] - r[2])**2)
        if(d!=0): EP += ((1/4*pi*E0)*sys[i].charge/d)
    return EP
#---------#




#---------#
def grav_field(sys,r):
    GF = [0,0,0]
    for i in range(len(sys)):
        d = math.sqrt((sys[i].pos[0] - r[0])**2+(sys[i].pos[1] - r[1])**2+(sys[i].pos[2] - r[2])**2)
        if d!=0:
            d_vec = [sys[i].pos[0] -r[0], sys[i].pos[1]-r[1],sys[i].pos[2]-r[2]]
            for j in range(3):
                GF[j] += G* sys[i].mass*d_vec[j]/d**3
    return GF
#---------#



#---------#
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
#---------#



#---------#
def elec_field(sys,r):
    E=[0,0,0]
    for i in range(len(sys)):
        d = math.sqrt((sys[i].pos[0] - r[0])**2+(sys[i].pos[1] - r[1])**2+(sys[i].pos[2] - r[2])**2)
        print("d :" + str(d))
        if d !=0:
            d_vec = [ r[0]-sys[i].pos[0] ,r[1]-sys[i].pos[1],r[2]-sys[i].pos[2]]
            for j in range(3):
                E[j]+=(1/(4*pi*E0))*sys[i].charge * d_vec[j] / d**3
    return E
#---------#



#---------#
def collision(a,b, file1):
    anga = math.acos((((b.pos[0] - a.pos[0])*a.vel[0])+((b.pos[1] - a.pos[1])*a.vel[1])+((b.pos[2] - a.pos[2])*a.vel[2]))/math.sqrt(((b.pos[0] - a.pos[0])**2+(b.pos[1] - a.pos[1])**2+(b.pos[2] - a.pos[2])**2)*((a.vel[0])**2+(a.vel[1])**2+(a.vel[2])**2)))
    angb = math.acos((((a.pos[0] - b.pos[0])*b.vel[0])+((a.pos[1] - b.pos[1])*b.vel[1])+((a.pos[2] - b.pos[2])*b.vel[2]))/math.sqrt(((a.pos[0] - b.pos[0])**2+(a.pos[1] - b.pos[1])**2+(a.pos[2] - b.pos[2])**2)*((b.vel[0])**2+(b.vel[1])**2+(b.vel[2])**2)))

    file1.write(str( anga)+'\t'+str(angb) +'\t'+str(math.cos(anga))+'\t'+str(math.cos(angb))+'\n')

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
#---------#



#---------#
class Emag_particle:
    def __init__(self):
          self.mass= 1
          self.radius =1
          self.pos = [0,0,0]
          self.prev_pos = [0,0,0]
          self.vel = [0,0,0]
          self.vel_mag = math.sqrt((self.vel[0])**2 + (self.vel[1])**2 + (self.vel[2])**2)
          self.acc = [0,0,0]
          self.KE = 0.5*self.mass*(self.vel[0]**2+self.vel[1]**2+self.vel[2]**2)
          self.charge = 1

    def data(self):
        self.charge = float(input('enter charge: '))
        self.mass = float(input('enter mass : '))
        self.radius = float(input('enter radius : '))
        self.KE = 0.5*self.mass*(self.vel[0]**2+self.vel[1]**2+self.vel[2]**2)
        for i in range(3):
          self.pos[i] = float(input('enter initial ' + str(i) +' position  : '))
          self.vel[i] = float(input('enter initial ' + str(i) +' velocity : '))
        self.vel_mag = math.sqrt((self.vel[0])**2 + (self.vel[1])**2 + (self.vel[2])**2)

#---------#







#---------#
def setup_system_to_simulate ():
    systm = []
    file2 = open('input_set.txt','r')
    #num = input('enter number of particles: ')
    file2_lines = file2.readlines()
    num = int(file2_lines[0])
    #print num
    #read_choice =int(input('1.enter manually 2.read from input file:'))
    read_choice = int(file2_lines[1])
    endtime = float(file2_lines[2])
    timeres = float(file2_lines[3])
    throttle_ip = float(file2_lines[4])
    file_name = str(file2_lines[5].rstrip())
    file2.close()
    print(read_choice)
    if read_choice ==2:
        #file_name = str(input('enter input data file name: '))
        #file_name = str(file2_lines[2].rstrip())
        print(file_name)
        file1= open(file_name,'r')
        read_lines = file1.readlines()
        data = [[read_lines[i+(j*5)] for i in range(5)] for j in range(num)]
        print(data)
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
            #ins.prev_pos = ins.pos
            ins.vel = data[i][4].split(',')
            ins.vel = [float(ins.vel[j]) for j in range(3)]
            ins.KE = 0.5*ins.mass*(ins.vel[0]**2+ins.vel[1]**2+ins.vel[2]**2)
            systm.append(ins)

    timesteps = endtime/timeres
    print( timesteps)
    return systm, int(num), endtime, timeres, int(timesteps), throttle_ip
#---------#



#---------#
def setup_x_y_z_t_plots(num):
    xplot =[]
    yplot =[]
    zplot =[]
    timeplot = []
    for i in range(num):
        xplot.append([])
        yplot.append([])
        zplot.append([])
    return xplot, yplot, zplot, timeplot
#---------#



#---------#
def setup_collision_matrix(num):
    col_matrix = []
    col_matrix=[[False for i in range(num)] for j in range(num)]
    return col_matrix
#---------#



#---------#
def init_prev_pos(num, systm, timeres):
    for partcle in range (num):
        for cod in range(3):
            systm[partcle].prev_pos[cod] = systm[partcle].pos[cod]
            systm[partcle].pos[cod] = systm[partcle].pos[cod]+(systm[partcle].vel[cod]*timeres)
        print(str(systm[partcle].prev_pos)+str(systm[partcle].pos))
    return systm
#---------#



#---------#
def update_xyzt_plot(i, timeres, thrle_cnt, throttle_ip, num, xplot, yplot, zplot, timeplot, systm):
    if  thrle_cnt == throttle_ip :
        thrle_cnt = 0
        for i in range(num):
            xplot[i].append(systm[i].pos[0])
            yplot[i].append(systm[i].pos[1])
            zplot[i].append(systm[i].pos[2])
    thrle_cnt += 1
    timeplot.append(i*timeres)
    print( 'time: '+ str(i*timeres))
    return xplot, yplot, zplot, timeplot, thrle_cnt
#---------#



#---------#
def solve_G_accl(systm, TKE, PE, num):
    for part  in range(num):
        TKE += systm[part].KE
        for accord in range(3):
            for ot_part in range(num):
                if(ot_part != part):
                    r = math.sqrt((systm[part].pos[0] - systm[ot_part].pos[0])**2 + (systm[part].pos[1] - systm[ot_part].pos[1])**2 +(systm[part].pos[2] - systm[ot_part].pos[2])**2 )                    
                    systm[part].acc[accord] = -(systm[part].acc[accord] + (G*systm[ot_part].mass*(systm[part].pos[accord] - systm[ot_part].pos[accord])/(r)**3))
                    PE -= G*(systm[part].mass*systm[ot_part].mass)/2*r
    return TKE, PE, systm
#---------#



#---------#
def solve_complete_classical_accelaration(num, TKE, PE, systm):
    E_PE = 0
    G_PE = 0
    for part1 in range(num):
        TKE += systm[part1].KE
        B = mag_field(systm,systm[part1].pos)
        E = elec_field(systm, systm[part1].pos)
        print("Elec accel: "+ str(part1) + " : " + str(systm[part1].charge * np.array(E) / systm[part1].mass))
        GF = grav_field(systm,systm[part1].pos)
        print("Grav acc: " + str(part1) + ": " + str(GF))
        v_x_B = cross(systm[part1].vel,B)
        E_PE += stat_elec_poten(systm, systm[part1].pos)
        G_PE += grav_poten(systm, systm[part1].pos)
        PE += E_PE + G_PE
        print("EMAG force: " + str(part1) + ": " + str(np.array(E) + np.array(v_x_B)))
        for cord1 in range(3):
            systm[part1].acc[cord1] = ((systm[part1].charge/systm[part1].mass)*(E[cord1]+v_x_B[cord1])) + GF[cord1]
    return systm, TKE, PE
#---------#


#---------#
def solve_sr_correct_accelaration(num, TKE, PE, systm):
    E_PE = 0
    G_PE = 0
    for part1 in range(num):
        TKE += systm[part1].KE
        B = mag_field(systm,systm[part1].pos)
        E = elec_field(systm, systm[part1].pos)  
        print("Elec accel: "+ str(part1) + " : " + str(systm[part1].charge * np.array(E) / systm[part1].mass))
        GF = grav_field(systm,systm[part1].pos)
        print("Grav acc: " + str(part1) + ": " + str(GF))
        v_x_B = cross(systm[part1].vel,B)
        E_PE += stat_elec_poten(systm, systm[part1].pos)
        G_PE += grav_poten(systm, systm[part1].pos)
        PE += E_PE + G_PE
        gmma = lorentz_factor(systm[part1])
        mink_force = gmma * (np.array(E) + np.array(v_x_B))
        print("Mink force: " + str(part1) + ": " + str(mink_force))
        for cord1 in range(3):
            systm[part1].acc[cord1] = ((systm[part1].charge/systm[part1].mass)*(mink_force[cord1])) + GF[cord1]
    return systm, TKE, PE    
#---------#


# #---------#
# def solve_sr_correct_accelaration(num, TKE, PE, systm):
#     E_PE = 0
#     G_PE = 0
#     for part1 in range(num):
#         TKE += systm[part1].KE
#         B = mag_field(systm,systm[part1].pos)
#         E = elec_field(systm, systm[part1].pos)  
#         print("Elec accel: "+ str(part1) + " : " + str(systm[part1].charge * np.array(E) / systm[part1].mass))
#         GF = grav_field(systm,systm[part1].pos)
#         print("Grav acc: " + str(part1) + ": " + str(GF))
#         E_PE += stat_elec_poten(systm, systm[part1].pos)
#         G_PE += grav_poten(systm, systm[part1].pos)
#         PE += E_PE + G_PE
#         four_vel = constrt_4_vel(systm[part1])
#         F_tensor = constrt_EM_tensor(E, B, systm[part1])
#         mink_force = minkowsky_force(F_tensor, four_vel, systm[part1])
#         print("Mink force: " + str(part1) + ": " + str(mink_force))
#         for cord1 in range(3):
#             systm[part1].acc[cord1] = ((systm[part1].charge/systm[part1].mass)*(mink_force[cord1+1])) + GF[cord1]
#     return systm, TKE, PE    
# #---------#




#---------#
def solve_collision(systm, num, col_matrix, file1):
    for part  in range(num):
        for ot_part in range(num):
            if(ot_part != part):
                r = math.sqrt((systm[part].pos[0] - systm[ot_part].pos[0])**2 + (systm[part].pos[1] - systm[ot_part].pos[1])**2 +(systm[part].pos[2] - systm[ot_part].pos[2])**2 )
                if r <= systm[part].radius+systm[ot_part].radius  and col_matrix[part][ot_part] ==False:
                    print(('collision'))
                    file1.write('collision\n' + str(part)+'\t'+str(ot_part)+'\n')
                    temp_col =collision(systm[part],systm[ot_part], file1)
                    systm[part] = temp_col[0]
                    systm[ot_part] = temp_col[1]
                    col_matrix[part][ot_part] = True
                    col_matrix[ot_part][part] = True
                if col_matrix[part][ot_part] ==True and r > systm[part].radius+systm[ot_part].radius:
                    col_matrix[part][ot_part] =False
                    col_matrix[ot_part][part] = False
  
        print(  'pos of '+ str(part) +' = ',systm[part].pos, 'vel'+str(part)+'= ',systm[part].vel  ,' acc'+str(part)+'= ',systm[part].acc)
    return col_matrix, systm
#---------#



#---------#    
def update_sys_state(systm, timeres, num, temp_prevpos):    
    for par in range(num):
        for coord in range(3):
            temp_prevpos[par][coord] = systm[par].pos[coord]
            systm[par].pos[coord] = (2*systm[par].pos[coord]) - systm[par].prev_pos[coord]+ (timeres**2*systm[par].acc[coord])
            systm[par].vel[coord] = (systm[par].pos[coord] - systm[par].prev_pos[coord])/(2*timeres)
            systm[par].prev_pos[coord] = temp_prevpos[par][coord]

        print(str(systm[par].vel))    
    return systm
#---------#





#---------#
def setup_system():
    
    systm, num, endtime, timeres, timesteps, throttle_ip = setup_system_to_simulate()  
        
    xplot, yplot, zplot, timeplot = setup_x_y_z_t_plots(num)
    
    col_matrix = setup_collision_matrix(num)
    
    systm = init_prev_pos(num, systm, timeres)
    
    temp_prevpos = [[0 for i in range(3)]for j in range(num)]

    return  temp_prevpos, systm, num, endtime, timeres, timesteps, throttle_ip, xplot, yplot, zplot, timeplot,col_matrix 
#---------#
 


#---------#           
def run_sim(temp_prevpos, file1, timesteps, throttle_ip, num, col_matrix, xplot, yplot, zplot, timeplot, timeres, systm, endtime):   
    thrle_cnt =0
    for  i in range(timesteps):
    
        file1.write('time: '+ str(i*timeres)+'\n')
        PE = 0
        TKE = 0
        TE = 0
        
        xplot, yplot, zplot, timeplot, thrle_cnt = update_xyzt_plot( i, timeres, thrle_cnt, throttle_ip, num, xplot, yplot, zplot, timeplot, systm)
        
        
        #systm, TKE, PE = solve_complete_classical_accelaration(num, TKE, PE, systm)
        
        systm, TKE, PE = solve_sr_correct_accelaration(num, TKE, PE, systm)
        
        col_matrix, systm = solve_collision(systm, num, col_matrix, file1)
    
        file1.write('TKE: '+str(TKE)+'\t'+'PE: '+str(PE)+'\t'+'TE: '+str(TE)+'\n')
    
        systm = update_sys_state(systm, timeres, num, temp_prevpos)
        
        
    print ('throttle:', thrle_cnt)
    return xplot, yplot, zplot
#---------#
    



#---------#
def run_sim_w_log_file(temp_prevpos, timesteps, throttle_ip, num, col_matrix, xplot, yplot, zplot, timeplot, timeres, systm, endtime):   
    file1 = open('partdata.txt','w')
    xplot, yplot, zplot = run_sim(temp_prevpos, file1, timesteps, throttle_ip, num, col_matrix, xplot, yplot, zplot, timeplot, timeres, systm, endtime)  
    file1.close()
    return xplot, yplot, zplot
#---------#




#---------#
def gen_1_plot_for_all_particles(num,  xplot, yplot, zplot):
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
  
    for i in range(num):
        ax.plot(xplot[i], yplot[i], zplot[i])
        
    plt.show()
#---------#
 

#---------#    
def init(lines, pts):
    for line, pt in zip(lines, pts):
        line.set_data([], [])
        line.set_3d_properties([])

        pt.set_data([], [])
        pt.set_3d_properties([])
    return lines + pts
#---------#



#---------#
def animate(i, ax3, xplot, yplot, zplot):
    # set up lines and points
    colors = plt.cm.jet(np.linspace(0, 1, num))
    lines = sum([ax3.plot([], [], [], '-', c=c)  for c in colors], [])
    pts = sum([ax3.plot([], [], [], 'o', c=c)  for c in colors], [])  
    for line, pt, xi, yi, zi in zip(lines, pts, xplot, yplot, zplot):
        #print 'entered'
        x, y, z = xi[:i], yi[:i], zi[:i]
        line.set_data(x, y)
        line.set_3d_properties(z)

        pt.set_data(x[-1:], y[-1:])
        pt.set_3d_properties(z[-1:])

    #ax.view_init(30, 0.3 * i)
    #fig.canvas.draw()
    return lines + pts
#---------#


   
#---------#
def traj_animation( throttle_ip, timesteps, xplot, yplot, zplot):
        
    fig3 = plt.figure()
    ax3 = p3.Axes3D(fig3)
    
  
    ax3.set_xlim(min([min(xplot[i]) for i in range(num)]),max([max(xplot[i]) for i in range(num)]))
    ax3.set_xlabel('X')
    ax3.set_ylim(min([min(yplot[i]) for i in range(num)]),max([max(yplot[i]) for i in range(num)]))
    ax3.set_ylabel('Y')
    ax3.set_zlim(min([min(zplot[i]) for i in range(num)]),max([max(zplot[i]) for i in range(num)]))
    ax3.set_zlabel('Z')
        
    # animation function.  This will be called sequentially with the frame number
    animation.FuncAnimation(fig3, animate, init_func=init, frames= int(timesteps/throttle_ip), interval=1, blit=False)   
    plt.show()
#---------#
   
    
    
#---------#   
if __name__ == '__main__':

    temp_prevpos, systm, num, endtime, timeres, timesteps, throttle_ip, xplot, yplot, zplot, timeplot, col_matrix = setup_system()

    start_time = time.time()
    xplot, yplot, zplot = run_sim_w_log_file(temp_prevpos, timesteps, throttle_ip, num, col_matrix, xplot, yplot, zplot, timeplot, timeres, systm, endtime)
    end_time = time.time()
    
    print(end_time - start_time)
     
    gen_1_plot_for_all_particles(num,  xplot, yplot, zplot)
    
    process = psutil.Process(os.getpid())
    print('Proc mem usage: '+ str(process.memory_info().rss/1024**2) + ' MBytes')  
    file_stats = os.stat(__file__)
    print(file_stats)
    print(f'File Size in KiloBytes is {file_stats.st_size / 1024}')    
    #traj_animation( throttle_ip, timesteps, xplot, yplot, zplot)
#---------#  
 
    
 
    