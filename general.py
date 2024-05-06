import math
import fluids as fds
import ht as ht
import heat_transfer as bht

from CoolProp.CoolProp import PropsSI

class duct:
    def __init__(self, type="tube", diameter=0, width=0, height=0,length=0,thickness=0,roughness=0,fluid="water",mdot=0,T_fluid_m=25+273.15,P_m=1.5):
        self.type = type
        self.Di = diameter
        self.W = width
        self.H = height
        self.L = length
        self.t = thickness

        if self.type == "rectangle":
            self.csA = self.W*self.H

            self.HD = (2*self.W*self.H)/(self.W+self.H)

        else:
            # Cross section area
            self.csA = math.pi*(self.Di/2)**2

            # Hydraulic diameter
            self.HD = self.Di

        if roughness==0:
            self.eD = fds.core.relative_roughness(self.HD)
        else:
            self.eD = fds.core.relative_roughness(self.HD,roughness)


        self.fluid = fluid
        self.T_fluid_m = T_fluid_m
        self.P_m = P_m

        self.mdot = mdot # kg/s
        self.T_fluid_m = T_fluid_m # K
        self.P_m = P_m # bar

        if self.mdot!=0: 
            self.compute_flow()

    def change_flow_rate(self,mdot,T_fluid_m=None,P_m=None):
        self.mdot = mdot
        if T_fluid_m!=None:
            self.T_fluid_m = T_fluid_m
        if P_m!=None:
            self.P_m = P_m

        self.compute_flow()
    
    def compute_flow(self):

        rho = PropsSI('D', 'P', self.P_m*1E5, 'T', self.T_fluid_m, self.fluid) # kg/m3
        mu = PropsSI('V', 'P', self.P_m*1E5, 'T', self.T_fluid_m, self.fluid) # Pa.s
        nu = mu/rho
        Pr = PropsSI('PRANDTL', 'P', self.P_m*1E5, 'T', self.T_fluid_m, self.fluid)
        k = PropsSI('CONDUCTIVITY', 'P', self.P_m*1E5, 'T', self.T_fluid_m, self.fluid) # conductivité thermique de l'eau

        # flow rate
        self.fr = self.mdot/rho
        # linear velocity
        self.lv = self.fr/self.csA
        self.Re = fds.Reynolds(self.lv,self.HD,rho,mu)
        self.Nu_int = ht.conv_internal.Nu_conv_internal(self.Re,Pr,self.eD,self.HD)
        self.h_int = (k*self.Nu_int)/self.HD


    def set_thermal(self,k_tube,e_ins,k_ins):
        self.k_tube = k_tube
        self.e_ins = e_ins
        self.k_ins = k_ins

    def compute_heat_transfer(self,T_fluid_m,T_amb, h_ext):

        if self.type == "rectangle":
            # Calculate the inner perimeter
            p_int = 2 * (self.W + self.H)
        elif self.type == "tube":
            # Calculate the inner perimeter
            p_int = 2 * math.pi * self.Di/2

        self.D1 = self.HD
        self.D2 = self.D1 + 2*self.t 
        self.D3 = self.D2 + 2*self.e_ins

        R_int = 1 / (self.h_int * p_int * self.L)
        R_tube = math.log(self.D2 / self.D1) / (2 * math.pi * self.k_tube * self.L)
        R_ins = math.log(self.D3 / self.D2) / (2 * math.pi * self.k_ins * self.L)
        R_ext = 1 / (h_ext * 2 * math.pi * (self.D3/2) * self.L)

        # Calculate the heat transfer coefficient K with insulation adaptation
        R_tot = R_int + R_tube + R_ins + R_ext
        self.R_tot = R_tot

        K = 1 / R_tot
        
        # Calculate the heat transfer Q
        Q = - K * (T_fluid_m - T_amb)
        
        # Calculate the wall temperature at location 'a', Twa
        # Twa = (self.h_int * T_fluid_m * math.pi * self.D1 - Q) / (math.pi * self.D1 * self.h_int)
        
        # Calculate the wall temperature at the outer surface of insulation, Twb
        Twb = (R_int+R_tube+R_ins)*Q + T_fluid_m
        # Twb = (h_ext * T_amb * math.pi * (self.D2 + 2 * self.e_ins) + Q) / (math.pi * (self.D2 + 2 * self.e_ins) * h_ext)

        # self.T_tube_int = Twa
        self.T_ins = Twb
        
        return Q
    
    def compute_heat_transfer_cylinder_free_convection(self,T_fluid_m,T_amb):
        Q0 = 0.
        Q1 = 1000.
        self.T_ins = T_fluid_m

        while abs(Q1-Q0)>0.001:
            Q0 = Q1
            self.h_ext = bht.back_h_cylinder(self.T_ins,T_amb,self.HD + 2*self.t + 2*self.e_ins) + 1*1E-3
            Q1 = self.compute_heat_transfer(T_fluid_m,T_amb,self.h_ext)

        return Q1