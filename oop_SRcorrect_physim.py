#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 12:11:34 2021

@author: sashwathnalinkanth

references: 
    1. https://people.nscl.msu.edu/~witek/Classes/PHY432/Potentials%20and%20Fields/Potentials2.pdf
    2. https://www.physics.princeton.edu//~mcdonald/examples/orbitdecay.pdf
    3. https://en.wikipedia.org/wiki/Abraham%E2%80%93Lorentz_force
    4. http://optics.hanyang.ac.kr/~shsong/Chapter%2011.%20Griffiths-Radiation-Point%20charge%20radiation.pdf
    5. INTRODUCTION TO ELECTRODYNAMICS, 4TH EDITION DAVID J. GRIFFITHS
    6. Introduction to Classical Mechanics: With Problems and Solutions Textbook by David J. Morin
    7. http://www.bu.edu/simulation/publications/dcole/PDF/DCColePhysRev4thWeb.pdf
    8. Singal, A. (2017). Radiation reaction from electromagnetic fields in the neighborhood of a point charge. The American Journal of Physics, 85(3), 202-206.


"""

#---------#
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import time
import math
import numpy as np
import os, psutil
#---------#

#---------#
E0= 8.854e-12
pi= np.pi
Mu0 = 4e-7*pi
G = 6.674e-11
C_EM = 3e8
#---------#
 
#---------#
class Emag_particle:
    
    #---------#
    E0= 8.854e-12
    pi= np.pi
    Mu0 = 4e-7*pi
    G = 6.674e-11
    C_EM = 3e8
    #---------#
    
    #---------#
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
          self.E_PE = 0
          self.G_PE = 0
          self.B = [0,0,0]
          self.E = [0,0,0]
          self.GF = [0,0,0]
          self.v_x_B = [0,0,0]
          self.gmma = 0
          self.four_vel = [0,0,0,0]
          self.F_tensor = np.array([])
          self.mink_force = [0,0,0,0]
          self.acc_hist = [[0,0,0]]
          self.a_dot = [0,0,0]
          self.prev_acc = [0,0,0]
          self.F_rad_reactn = [0,0,0]
    #---------#

    #---------#
    def get_usr_data(self):
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
    def upd_fields(self, systm):
        self.B = systm.mag_field(self.pos)
        self.E = systm.elec_field(self.pos)
        self.GF = systm.grav_field(self.pos)
        self.v_x_B = systm.cross(self.vel,self.B)
    #---------#
    
    #---------#
    def upd_retd_fields(self, systm, retd_E, retd_B):
        self.E , self.B = retd_E, retd_B
        self.GF = systm.grav_field(self.pos)
        self.v_x_B = systm.cross(self.vel,self.B)
    #---------#  
        
    #---------#
    def init_prev_pos(self, timeres):
        for cod in range(3):
            self.prev_pos[cod] = self.pos[cod]
            self.pos[cod] = self.pos[cod]+(self.vel[cod]*timeres)
        print(str(self.prev_pos)+str(self.pos))
    #---------#
    
    #---------#
    def solveNupd_complete_classical_accelaration(self, systm):
        self.upd_fields(self, systm)
        self.acc = (self.charge/self.mass)*(self.E + self.v_x_B) + self.GF
    #---------#

    #---------#    
    def upd_E_G_potential_energy(self, systm):
        self.E_PE = systm.stat_elec_poten(systm, self.pos)
        self.G_PE = systm.grav_poten(systm, self.pos)
    #---------#
       
    #---------#
    def upd_lorentz_factor(self):
        self.gmma = 1/(math.sqrt(1-(self.vel_mag**2/self.C_EM**2)))
    #---------#
       
    #---------#
    def constrt_4_vel(self):
        self.upd_lorentz_factor()
        for i in range(4):
           if(i == 0): self.four_vel.append(self.gmma * self.C_EM)
           else: self.four_vel.append(self.gmma * self.vel[i-1])
        self.four_vel = np.array(self.four_vel)
    #---------#
       
    #---------#
    def constrt_EM_tensor(self):
        row1 = [0,                 (self.E[0]/self.C_EM), (self.E[1]/self.C_EM), (self.E[2]/self.C_EM)]
        row2 = [-(self.E[0]/self.C_EM),                0,        self.B[2],       -self.B[1]]
        row3 = [-(self.E[1]/self.C_EM),       -self.B[2],                0,        self.B[0]]
        row4 = [-(self.E[2]/self.C_EM),        self.B[1],       -self.B[0],                0]
        self.F_tensor = np.array([row1, row2, row3, row4])
    #---------#
       
    #---------#
    def minkowsky_force(self):
        self.constrt_4_vel(self)
        self.constrt_EM_tensor(self)
        self.mink_force = self.charge * np.matmul(self.four_vel, self.F_tensor)
    #---------#
    
    #---------#
    def solve_G_accl(self, systm):
        self.GF = systm.grav_field(systm,self.pos)
        self.acc = self.GF
    #---------#
    
    #---------#
    def solve_rad_reactn(self, timeres):
        N = len(self.acc_hist)
        if(N > 2):
            #self.a_dot = ( np.array(self.acc_hist[N-1]) - np.array(self.acc_hist[N-2])) / timeres
            self.a_dot = (np.array(self.acc) - np.array(self.prev_acc) ) / timeres
        else:
            self.a_dot = [0,0,0]
        print("Acc dot: " + str(np.array(self.acc)) + str(np.array(self.prev_acc)))
        self.F_rad_reactn = (Mu0 * self.charge**2) * np.array(self.a_dot) / (6 * pi * C_EM)
    #---------#
            
    #---------#
    def solve_sr_correct_accelaration(self, systm):
        temp_acc = [0,0,0]
        self.upd_fields(systm) 
        self.upd_lorentz_factor()
        self.mink_force = self.gmma * np.array(self.E + self.v_x_B)
        temp_prev_acc = self.acc
        for cord1 in range(3):           
            #self.acc[cord1] = ((self.charge/self.mass)*(self.mink_force[cord1])) + self.GF[cord1] - (self.F_rad_reactn[cord1] / self.mass)
            temp_acc[cord1] = ((self.charge/self.mass)*(self.mink_force[cord1])) + self.GF[cord1] 
        self.prev_acc = temp_prev_acc
        self.acc_hist.append(temp_prev_acc)
        self.acc = temp_acc
        print("acc to append (prev) : "+ str(self) + " : " + str(self.prev_acc))
        print("acc after appnd (cur): "+ str(self) + " : " + str(self.acc))
    #---------#
    
    #---------#
    def solve_sr_correct_accelaration_w_rad_react(self, systm):
        temp_acc = [0,0,0]
        self.upd_fields(systm) 
        self.upd_lorentz_factor()
        self.mink_force = self.gmma * np.array(self.E + self.v_x_B)
        temp_prev_acc = self.acc
        for cord1 in range(3):           
            #self.acc[cord1] = ((self.charge/self.mass)*(self.mink_force[cord1])) + self.GF[cord1] - (self.F_rad_reactn[cord1] / self.mass)
            temp_acc[cord1] = ((self.charge/self.mass)*(self.mink_force[cord1])) + self.GF[cord1] + (self.F_rad_reactn[cord1] / self.mass)
        self.prev_acc = temp_prev_acc
        self.acc_hist.append(temp_prev_acc)
        self.acc = temp_acc
        print("acc to append (prev) : "+ str(self) + " : " + str(self.prev_acc))
        print("acc after appnd (cur): "+ str(self) + " : " + str(self.acc))
    #---------# 
    
    #---------#    
    def update_particle_state(self, timeres):  
        temp_prevpos = [0,0,0]
        self.solve_rad_reactn(timeres)
        for coord in range(3):
            temp_prevpos[coord] = self.pos[coord]
            self.pos[coord] = (2 * self.pos[coord]) - self.prev_pos[coord]+ (timeres**2 * self.acc[coord])
            self.vel[coord] = (self.pos[coord] - self.prev_pos[coord])/(2 * timeres)
            self.prev_pos[coord] = temp_prevpos[coord]

        print(str(self.vel))    
    #---------#
     
#---------#




#---------#
class Physical_System_of_particles:
    
    #---------#
    systm = [] # A list objects of type Emag_particle()
    num = 0
    col_matrix = []
    ret_t = 0
    retd_r = [0,0,0]
    retd_r_mag = 0
    retd_r_cap = [0,0,0]
    u_vec = [0,0,0]
    #---------#
    
    #---------#
    def setup_system_of_particles_to_simulate(self, sim_obj):
        print('choice: ' + str(sim_obj.read_choice))
        self.num = sim_obj.num
        self.col_matrix = [[False for i in range(self.num)] for j in range(self.num)]
        if sim_obj.read_choice ==2:            
            print("Reading input file: " + sim_obj.file_name)
            file1= open(sim_obj.file_name,'r')
            read_lines = file1.readlines()
            data = [[read_lines[i+(j*5)] for i in range(5)] for j in range(sim_obj.num)]
            print(data)
            file1.close() 
            print("Reading Complete. Closing input file: " + sim_obj.file_name)
        for i in range(sim_obj.num):
            if(sim_obj.read_choice == 1):
                ins = Emag_particle()
                ins.get_usr_data()
                ins.init_prev_pos(sim_obj.timeres)
                self.systm.append(ins)
            elif sim_obj.read_choice ==2:
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
                ins.init_prev_pos(sim_obj.timeres)
                self.systm.append(ins)    
    #---------#
    
    #---------#    
    def update_sys_state(self, timeres):    
        for par in self.systm:
             par.update_particle_state(timeres)
    #---------#
    
    #---------#
    def solve_retrd_time(self, par, r, t):
        r_mag = math.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
        t1 = (C_EM**2 * t) - (self.dot(par.vel, r))
        if (t1**2 + (C_EM**2 -par.vel_mag**2)*(r_mag**2 - (C_EM**2 * t**2))) >= 0:
            t2 = math.sqrt(t1**2 + (C_EM**2 -par.vel_mag**2)*(r_mag**2 - (C_EM**2 * t**2)))
        else:
            t2 = 0           
        t3 = C_EM**2 - par.vel_mag**2
        self.ret_t =  (t1 -t2)/t3
    #---------#
    
    #---------#
    def solve_ret_rNcap(self, par, r):
        self.retd_r[0] = r[0] - par.vel[0]*self.ret_t
        self.retd_r[1] = r[1] - par.vel[1]*self.ret_t
        self.retd_r[2] = r[2] - par.vel[2]*self.ret_t
        self.ret_r_mag = self.vec_mag(self.retd_r)
        if(self.ret_r_mag == 0):
            self.retd_r_cap = [0,0,0]
        else:
            self.retd_r_cap[0] = self.retd_r[0]/self.ret_r_mag
            self.retd_r_cap[1] = self.retd_r[1]/self.ret_r_mag
            self.retd_r_cap[2] = self.retd_r[2]/self.ret_r_mag
    #---------#

    #---------#    
    def solve_u_vec(self, par):
        self.u_vec[0] = par.C_EM*self.retd_r_cap[0] - par.vel[0]
        self.u_vec[1] = par.C_EM*self.retd_r_cap[1] - par.vel[1]
        self.u_vec[2] = par.C_EM*self.retd_r_cap[2] - par.vel[2]
    #---------#
    
    #---------#    
    def solve_retd_acc(self, par):
        print(self.ret_t)
        if(self.ret_t > 0): 
            self.retd_acc = par.acc_hist[int(self.ret_t)]
        else:
            self.retd_acc = [0,0,0]
    #---------#  
    
    #---------#  
    def solve_retd_elec_field(self, par, r, t):
        self.solve_retrd_time(par, r, t)
        self.solve_ret_rNcap(par, r)
        self.solve_u_vec(par)
        self.solve_retd_acc(par)
        print("ret acc:" + str(self.retd_acc))
        c2_v2xu = (C_EM**2 - par.vel_mag**2) * np.array(self.u_vec)
        uXa = self.cross(self.u_vec, self.retd_acc)
        rXuXa = self.cross(self.retd_r, uXa)
        r_dot_u_3 = (self.dot(self.retd_r, self.u_vec))**3
        charge_scalar = par.charge/(4*pi*E0)
        if(r_dot_u_3 != 0):
            r_mag_by_rdotu3 = self.ret_r_mag / r_dot_u_3
        else:
            r_mag_by_rdotu3 = 0
        retd_E = charge_scalar * r_mag_by_rdotu3 * (c2_v2xu + rXuXa)
        return retd_E
    #---------#  

    #---------#  
    def get_retd_elec_mag_field(self, par, r, t):
        retd_E = self.solve_retd_elec_field(par, r, t)
        retd_B = np.array(self.cross(self.retd_r_cap, retd_E))/C_EM
        return retd_E, retd_B
    #---------#  
        
    #---------#  
    def solve_sr_acc_sys_particles_w_retd_fields(self, t):
        for particle1 in self.systm:
            temp_E = [0,0,0]
            temp_B = [0,0,0]
            for other_particle in self.systm:
                if(particle1 != other_particle): 
                    retd_E, retd_B = self.get_retd_elec_mag_field(other_particle, particle1.pos, t)
                    for i in range(3):
                        temp_E[i] = temp_E[i] + retd_E[i]
                        temp_B[i] = temp_B[i] + retd_B[i]
            particle1.upd_retd_fields(self, temp_E, temp_B)
            particle1.solve_sr_correct_accelaration_w_rad_react(self)
    #---------#  
    
    #---------#
    def detect_collision(self, file1):
        for particle, part  in zip(self.systm, range(self.num)):
            for other_particle, ot_part in zip(self.systm, range(self.num)):
                if(other_particle != particle):
                    r = math.sqrt((particle.pos[0] - other_particle.pos[0])**2 + (particle.pos[1] - other_particle.pos[1])**2 +(particle.pos[2] - other_particle.pos[2])**2 )
                    if r <= particle.radius + other_particle.radius  and self.col_matrix[part][ot_part] == False:
                        print(('collision'))
                        file1.write('collision\n' + str(part)+'\t'+str(ot_part)+'\n')
                        temp_col = self.collision(particle, other_particle, file1)
                        self.systm[part] = temp_col[0]
                        self.systm[ot_part] = temp_col[1]
                        self.col_matrix[part][ot_part] = True
                        self.col_matrix[ot_part][part] = True
                    if self.col_matrix[part][ot_part] ==True and r > particle.radius + other_particle.radius:
                        self.col_matrix[part][ot_part] =False
                        self.col_matrix[ot_part][part] = False     
            print(  'pos of ' + str(part) + ' = ', particle.pos, 'vel' + str(part) + '= ', particle.vel  ,' acc' + str(part) + '= ', particle.acc)
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
       
    #---------#
    def cross(self, a, b):
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
    def dot(self, a, b):
        p = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]               
        return p
    #---------#  
    
    #---------#
    def vec_mag(self, a):
        p = math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)            
        return p
    #---------#      
     
    #---------#
    def grav_poten(self, r):
        GP =0
        for particle in self.systm:
            d = math.sqrt((particle.pos[0] - r[0])**2+(particle.pos[1] - r[1])**2+(particle.pos[2] - r[2])**2)
            if(d!=0): GP += (G*particle.mass/d)
        return GP
    #---------#
       
    #---------#
    def stat_elec_poten(self, r):
        EP =0
        for particle in self.systm:
            d = math.sqrt((particle.pos[0] - r[0])**2+(particle.pos[1] - r[1])**2+(particle.pos[2] - r[2])**2)
            if(d!=0): EP += ((1/4*pi*E0)*particle.charge/d)
        return EP
    #---------#
       
    #---------#
    def grav_field(self, r):
        GF = [0,0,0]
        for particle in self.systm:
            d = math.sqrt((particle.pos[0] - r[0])**2+(particle.pos[1] - r[1])**2+(particle.pos[2] - r[2])**2)
            if d!=0:
                d_vec = [particle.pos[0] -r[0], particle.pos[1]-r[1],particle.pos[2]-r[2]]
                for j in range(3):
                    GF[j] += G* particle.mass*d_vec[j]/d**3
        return GF
    #---------#
       
    #---------#
    def mag_field(self, r):
        B=[0,0,0]
        for particle in self.systm:
            d = math.sqrt((particle.pos[0] - r[0])**2+(particle.pos[1] - r[1])**2+(particle.pos[2] - r[2])**2)
            if d != 0:
                d_vec = [ r[0]-particle.pos[0] ,r[1]-particle.pos[1],r[2]-particle.pos[2]]
                v_x_d = self.cross(particle.vel, d_vec)
                for j in range(3):
                    B[j] += (Mu0/4/pi)*particle.charge*v_x_d[j]/d**3
        return B
    #---------#
        
    #---------#
    def elec_field(self, r):
        E=[0,0,0]
        for particle in self.systm:
            d = math.sqrt((particle.pos[0] - r[0])**2+(particle.pos[1] - r[1])**2+(particle.pos[2] - r[2])**2)
            if d !=0:
                d_vec = [ r[0]-particle.pos[0] ,r[1]-particle.pos[1],r[2]-particle.pos[2]]
                for j in range(3):
                    E[j]+=(1/(4*pi*E0))*particle.charge*d_vec[j]/d**3
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
    def solve_sr_acc_sys_particles(self):
        for particle in self.systm:
            particle.solve_sr_correct_accelaration(self)
    #---------#

#---------#





#---------#
class Simulation_object():
    
    #---------#
    endtime = 0 
    timeres = 0
    timesteps = 0
    cur_timestep = 0
    throttle_ip = 0
    xplot =[]
    yplot =[]
    zplot =[]
    timeplot = []   
    thrle_cnt = 0
    num = 0
    PE = 0
    TKE = 0
    TE = 0
    #---------#

    #---------#
    def read_sim_param_input_file(self, input_file):
            file2 = open(input_file,'r')
            file2_lines = file2.readlines()
            self.num = int(file2_lines[0])
            #-----1.enter manually 2.read from input file-----#
            self.read_choice = int(file2_lines[1])
            self.endtime = float(file2_lines[2])
            self.timeres = float(file2_lines[3])
            self.throttle_ip = float(file2_lines[4])
            self.file_name = str(file2_lines[5].rstrip())
            file2.close()
            self.timesteps = int(self.endtime/self.timeres)
            print(self.timesteps)
    #---------#
    
    #---------#
    def setup_x_y_z_t_plots(self, systm):
        print(systm.num)
        for i in range(systm.num):
            self.xplot.append([])
            self.yplot.append([])
            self.zplot.append([])
        print("Printint:" + str((self.xplot)))
    #---------#
    
    #---------#
    def update_xyzt_plot(self, system):
        if  self.thrle_cnt == self.throttle_ip :
            self.thrle_cnt = 0
            print(range(len(system.systm)))
            print(len(self.xplot))
            for particle in range(system.num):
                self.xplot[particle].append(system.systm[particle].pos[0])
                self.yplot[particle].append(system.systm[particle].pos[1])
                self.zplot[particle].append(system.systm[particle].pos[2])
        self.thrle_cnt += 1
        self.timeplot.append(self.cur_timestep * self.timeres)
        print( 'time: '+ str(self.cur_timestep * self.timeres))
    #---------#
       
    #---------#           
    def run_n_particle_sim(self, system, file1):   
        for  i in range(self.timesteps):
        
            self.cur_timestep = i    
        
            if (file1 != None): file1.write('time: '+ str(self.cur_timestep * self.timeres)+'\n')
            
            self.update_xyzt_plot(system)
            
            
            #system.solve_sr_acc_sys_particles()
            system.solve_sr_acc_sys_particles_w_retd_fields(self.cur_timestep)
            
            
            system.detect_collision(file1)
        
            if (file1 != None): file1.write('TKE: '+str(self.TKE)+'\t'+'PE: '+str(self.PE)+'\t'+'TE: '+str(self.TE)+'\n')
        
            system.update_sys_state(self.timeres)
            
            
        print ('throttle:', self.thrle_cnt)
    #---------#
           
    #---------#
    def run_n_particle_sim_w_log_file(self, systm):   
        file1 = open('partdata.txt','w')
        self.run_n_particle_sim(systm, file1)  
        file1.close()
    #---------#
    
    #---------#
    def gen_1_plot_for_all_particles(self):
        fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
      
        for i in range(self.num):
            ax.plot(self.xplot[i], self.yplot[i], self.zplot[i])
            
        plt.show()
    #---------#
           
    #---------#
    def anim_plot_for_all_particles(self):
        plt.ion()
        fig = plt.figure() 
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z') 

        for i in range(self.num):
            sc = ax.scatter(self.xplot[i], self.yplot[i], self.zplot[i])           
            fig.show()
            for plt_till in range(len(self.xplot[i])):
                plt.pause(0.001)
                #ax.plot(self.xplot[i][:plt_till], self.yplot[i][:plt_till], self.zplot[i][:plt_till])           
                sc._offsets3d = (self.xplot[i][:plt_till], self.yplot[i][:plt_till], self.zplot[i][:plt_till])
                plt.draw()
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
    def animate(self, i, ax3, xplot, yplot, zplot):
        # set up lines and points
        colors = plt.cm.jet(np.linspace(0, 1, self.num))
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
    def traj_animation(self, throttle_ip, timesteps, xplot, yplot, zplot):
            
        fig3 = plt.figure()
        ax3 = p3.Axes3D(fig3)
        
      
        ax3.set_xlim(min([min(xplot[i]) for i in range(self.num)]),max([max(xplot[i]) for i in range(self.num)]))
        ax3.set_xlabel('X')
        ax3.set_ylim(min([min(yplot[i]) for i in range(self.num)]),max([max(yplot[i]) for i in range(self.num)]))
        ax3.set_ylabel('Y')
        ax3.set_zlim(min([min(zplot[i]) for i in range(self.num)]),max([max(zplot[i]) for i in range(self.num)]))
        ax3.set_zlabel('Z')
            
        # animation function.  This will be called sequentially with the frame number
        animation.FuncAnimation(fig3, self.animate, init_func= self.init, frames= int(timesteps/throttle_ip), interval=1, blit=False)   
        plt.show()
    #---------#
    
#---------#      
        

    
#---------#   
if __name__ == '__main__':
    #-----Start Sim------#
    start_time = time.time()

    #------Setup Sim and System objects-----#
    sim = Simulation_object()
    system = Physical_System_of_particles()

    #------Create Sim and System params-----#    
    sim.read_sim_param_input_file("input_set.txt")
    system.setup_system_of_particles_to_simulate(sim)
    sim.setup_x_y_z_t_plots(system)
        
    #-------Run particle trajectory simulation-----#
    sim.run_n_particle_sim_w_log_file(system)
    
    #------Plot trajectories----------#
    sim.gen_1_plot_for_all_particles()
    
    #--------Animated plot----------#
    #sim.anim_plot_for_all_particles()

    #-----End of Sim--------#    
    end_time = time.time()    
    print( 'Tot run time: ' + str(end_time - start_time))
       
    #------Calculate memory Usage------#
    process = psutil.Process(os.getpid())
    print('Proc mem usage: '+ str(process.memory_info().rss/1024**2) + ' MBytes')  
    file_stats = os.stat(__file__)
    print(file_stats)
    print(f'File Size in KiloBytes is {file_stats.st_size / 1024}')
#---------#  
 
    
 
    