import math
import fluids as fds
import ht as ht

class duct:
    def __init__(self, type, diameter=0, width=0, height=0,length=0,thickness=0,eD=0,m_dot=0,fluid="water"):
        self.type = type
        self.Di = diameter
        self.W = width
        self.H = height
        self.L = length
        self.t = thickness
        self.eD = eD # roughness

        if self.type == "tube":
            # Cross section area
            self.csA = math.pi*(self.Di/2)**2

            # Hydraulic diameter
            self.HD = self.Di

        elif self.type == "rectangle":
            self.csA = self.W*self.H

            self.HD = (2*self.W*self.H)/(self.W+self.H)

        self.fluid = fluid

        self.m_dot = m_dot # kg/s

        if self.m_dot!=0: 

            rho = 997 # kg/m3
            nu = 0.896*1e-6
            eta = rho*nu
            Pr = 7.01
            k = 0.6 # conductivité thermique de l'eau

            # flow rate
            self.fr = self.m_dot/rho

            # linear velocity
            self.lv = self.fr/self.csA

            self.Re = fds.Reynolds(self.lv,self.HD,rho,eta)

            self.Nu_int = ht.conv_internal.Nu_conv_internal(self.Re,Pr,self.eD,self.HD)

            self.h_int = (k*self.Nu_int)/self.HD
        



