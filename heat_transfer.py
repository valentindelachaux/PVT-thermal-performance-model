import math
import openpyxl as opxl

import numpy as np
import scipy.constants as scc
import scipy.integrate as integrate
import ht
from CoolProp.CoolProp import PropsSI

from matplotlib import pyplot as plt

coeff_downward_cool_surface = 1.2

import scipy.constants as scc

def h_rad(eps,T_s,T_rad):
    """Calculates the radiative heat transfer coefficient between two surfaces
    
    Args:
        eps (float): Emissivity of the surface
        T_s (float): Temperature of the surface
        T_rad (float): Temperature of the radiating surface
        
    Returns:
        float: Radiative heat transfer coefficient"""
    
    return eps*scc.sigma*(T_s+T_rad)*(T_s**2+T_rad**2)

def air_rho(T):
    return -0.00439881*T+2.500535714

def air_c_p():
    return 1006

def air_mu(T):
    return (0.004791667*T+0.4065)*1e-5

def air_nu(T):
    return (0.008894048*T-1.097178571)*1e-5

def air_k(T):
    return 7.2607*1e-5*T+0.004365714

def air_Pr():
    return 0.711

def grashof(beta, DT, L_c, rho, mu):
    g = scc.g
    return (g*beta*DT*L_c**3*rho**2)/(mu**2)

def back_h_cylinder(T_tube,T_amb,D):
    T_film = (T_tube+T_amb)/2
    Gr = grashof(1/T_film, T_tube-T_amb, D, air_rho(T_film), air_mu(T_film))
    Nu = ht.conv_free_immersed.Nu_horizontal_cylinder(air_Pr(),abs(Gr))
    h = (air_k(T_film)/D)*Nu
    return h

def back_h_mixed(T_abs,T_amb,u_back,theta,longueur):
    h_back_free = back_h_simple(T_abs,T_amb,theta,longueur)
    h_back_forced = h_top_forced(T_abs,T_amb,u_back,longueur)
    return (h_back_free**3 + h_back_forced**3)**(1/3)

def back_h(T_abs,T_amb,theta,longueur,largeur,N_ailettes,a):
    DT = T_abs - T_amb
    T_mean = (T_abs+T_amb)/2

    g = scc.g

    rho = air_rho(T_mean)
    Cp = air_c_p()
    mu = air_mu(T_mean)
    nu = air_nu(T_mean)
    lambd = air_k(T_mean)
    alpha = (lambd)/(rho*Cp)
    Pr = air_Pr()
    beta = 1/T_mean
    
    """
    air = Mixture('air',T=T_mean,P=1e5)
    rho = air.rho
    Cp = air.Cp
    mu = air.mu
    nu = air.nu
    lambd = air.k
    alpha = (lambd)/(rho*Cp)
    Pr = (mu*Cp)/lambd
    beta = 1/T_mean
    """

    D = (largeur-N_ailettes*a)/(N_ailettes-1)
    
    Ra = ((rho**2)*g*math.cos(math.radians(theta))*beta*Cp*(D**4)*DT)/(mu*lambd*longueur)
    
    if DT > 0:
        Nu = (Ra/24)*(1-np.exp(-35/Ra))**(0.75)
    elif DT < 0:
        Nu = (-Ra/24)*(1-np.exp((coeff_downward_cool_surface*-35)/(-Ra)))**(0.75)

    h=(lambd/D)*Nu
    
    return h

def back_h_fins(T_abs,T_amb,theta,longueur,D,L_a):
    DT = T_abs - T_amb
    T_mean = (T_abs+T_amb)/2

    g = scc.g

    rho = air_rho(T_mean)
    Cp = air_c_p()
    mu = air_mu(T_mean)
    nu = air_nu(T_mean)
    lambd = air_k(T_mean)
    alpha = (lambd)/(rho*Cp)
    Pr = air_Pr()
    beta = 1/T_mean
    
    Ra = ((rho**2)*g*math.cos(math.radians(theta))*beta*Cp*(D**4)*DT)/(mu*lambd*longueur)

    Gr2=(g*beta*abs(DT)*D**3)*((L_a/longueur)**(1/2))*((D/L_a)**0.38)*(1/nu)**2
    Gr1=(g*beta*abs(DT)*D**4)/(math.sqrt(longueur*L_a)*nu**2)

    # if D>=0.050:
    #     return back_h_simple(T_abs,T_amb,theta, longueur)

    if DT<0 and theta<=30:
        crit=Gr2*Pr*math.sin((math.pi/2) - math.radians(theta))
        if crit<=20000:
            Nu = 0.0915*crit**0.436
            h = (lambd/D)*Nu
            return h
        else:
            return 2.
    else:
        crit=Gr1*Pr*math.cos((math.pi/2)-math.radians(theta))
        if crit<=250:
            Nu = 0.0929*crit**0.5
            h = (lambd/D)*Nu
            return h
        elif crit<=1000000:
            Nu = 0.2413*crit**(1/3)
            h = (lambd/D)*Nu
            return h
        else:
            return 2.
        
def ht_fins_forced_wiki(L_fin,D,u):
    return 2.11*((u/3.6)**0.71)*((D*1000)**0.44)*(L_fin*1000)**-0.14

def ht_fins_forced_wiki2(L_fin,D,u):
    return (2.47 - 2.55/(D*1000)**0.4)*((u/3.6)**0.9)*0.0872*(D*1000) + 4.31
  
def speed_natural_convection(T_abs,T_amb,theta,L):

    DT = T_abs - T_amb
    T_mean = (T_abs+T_amb)/2

    g = scc.g

    beta = 1/T_mean     
    Ra = Ra_L(T_abs,T_amb,theta,L)

    speed = ((scc.g*beta*abs(DT)*L)/Ra)**(1/2)

    return speed

def Ra_L(T_abs,T_amb,theta,L):

    DT = T_abs - T_amb

    T_mean = (T_abs+T_amb)/2

    g = scc.g

    rho = air_rho(T_mean)
    Cp = air_c_p()
    mu = air_mu(T_mean)
    nu = air_nu(T_mean)
    lambd = air_k(T_mean)
    alpha = (lambd)/(rho*Cp)
    Pr = air_Pr()
    beta = 1/T_mean

    Ra_L= (g*beta*math.cos(math.pi/2-math.radians(theta))*abs(DT)*(L**4))/(nu*alpha)

    return Ra_L

def back_h_simple(T_abs,T_amb,theta,longueur): # dans 'Inputs', theta est l'angle par rapport à l'horizontale donc c'est theta_p
    # températures en K
    # angle theta en °
    # longueur en m
    
    h_back_zero = 0.1

    DT = T_abs - T_amb

    T_mean = (T_abs+T_amb)/2

    g = scc.g

    rho = air_rho(T_mean)
    Cp = air_c_p()
    mu = air_mu(T_mean)
    nu = air_nu(T_mean)
    lambd = air_k(T_mean)
    alpha = (lambd)/(rho*Cp)
    Pr = air_Pr()
    beta = 1/T_mean

    # On vire ça ?
    # if abs(DT)<=0.05:
    #     return 0.5

    if DT<0:
        if theta>45:
            Ra_L=(g*beta*math.cos(math.pi/2-math.radians(theta))*abs(DT)*(longueur**4))/(nu*alpha)
            if Ra_L >= 1e4 and Ra_L <= 1e9:
                Nu_L = 0.68+0.67*Ra_L**(1/4)*(1+(0.492/Pr)**(9/16))**(-4/9)
                h = (lambd/longueur)*Nu_L
                return h
            elif Ra_L >= 1e9:
                Nu_L = 0.10*Ra_L**(1/3)
                h = (lambd/longueur)*Nu_L
                return h  
            else:
                print('back_h_simple DT<0 Ra_L < 1e4',Ra_L)
                return h_back_zero
        
        elif theta<=45 and theta>=2:
            Ra_L=(g*beta*math.sin(math.pi/2-math.radians(theta))*abs(DT)*(longueur**4))/(nu*alpha)
            if Ra_L>=1e5 and Ra_L<=1e11:
                Nu_L = 0.14*Ra_L**(1/3)*((1+0.0107*Pr)/(1+0.01*Pr))
                h = (lambd/longueur)*Nu_L
                return h
            else:
                print('back_h_simple DT<0 Ra_L<1e5',Ra_L)
                return h_back_zero

        else:
            print('theta',theta)
            print("theta should be greater than 2°")
            return "erreur"

    elif DT>0:
        if theta>=2:
            Ra_L=(g*beta*math.cos(math.pi/2-math.radians(theta))*abs(DT)*(longueur**4))/(nu*alpha)

            if Ra_L >= 1e5 and Ra_L <= 1e11:
                Nu_L = 0.68+0.67*Ra_L**(1/4)*(1+(0.492/Pr)**(9/16))**(-4/9)
                h = (lambd/longueur)*Nu_L
                return h
            else:
                print('back_h_simple DT>0 Ra_L<1e5',Ra_L)
                return h_back_zero
        else:
            print('back_h_simple DT>0 Ra_L',Ra_L)
            return h_back_zero

    else:
        print('DT',DT)
        return h_back_zero

#  local Nusselt number relations for a flat plate with a constant heat flux
# https://courses.ansys.com/wp-content/uploads/2021/02/LT4_C2_L3-Handout-v2.pdf

def Nu_forced_flat_plate_isoflux_lam(x,k,speed,nu,Pr): # 0.6 < Pr
    Re_x = (speed*x)/nu
    return (k/x)*0.453*Re_x**(1/2)*Pr**(1/3)

def Nu_forced_flat_plate_isoflux_turb(x, k,speed,nu,Pr): # 0.6 < Pr < 60
    Re_x = (speed*x)/nu
    return (k/x)*0.0308*Re_x**(4/5)*Pr**(1/3)

def h_top_custom(T_s,T_amb,speed,longueur):
    return 3+3*speed

def h_top_forced(T_s,T_amb,speed,longueur):

    T_mean = (T_s+T_amb)/2

    g = scc.g

    rho = air_rho(T_mean)
    Cp = air_c_p()
    mu = air_mu(T_mean)
    nu = air_nu(T_mean)
    lambd = air_k(T_mean)
    alpha = (lambd)/(rho*Cp)
    Pr = air_Pr()
    beta = 1/T_mean

    Re_c = 3.5 * 10**5

    x_c = (nu*Re_c)/speed

    if x_c < longueur:
        lam = integrate.quad(Nu_forced_flat_plate_isoflux_lam,0,x_c,args=(lambd,speed,nu,Pr))[0]
        turb = integrate.quad(Nu_forced_flat_plate_isoflux_turb,x_c,longueur,args=(lambd,speed,nu,Pr))[0]
    else:
        lam = integrate.quad(Nu_forced_flat_plate_isoflux_lam,0,longueur,args=(lambd,speed,nu,Pr))[0]
        turb = 0.

    return (1/longueur)*(lam+turb)

def h_top_forced_turbulent(T_s,T_amb,speed,longueur):

    T_mean = (T_s+T_amb)/2

    g = scc.g

    rho = air_rho(T_mean)
    Cp = air_c_p()
    mu = air_mu(T_mean)
    nu = air_nu(T_mean)
    lambd = air_k(T_mean)
    alpha = (lambd)/(rho*Cp)
    Pr = air_Pr()
    beta = 1/T_mean

    Re_c = 3.5 * 10**5

    x_c = (nu*Re_c)/speed

    if speed > 0.4:
        h_int = integrate.quad(Nu_forced_flat_plate_isoflux_turb,0,longueur,args=(lambd,speed,nu,Pr))[0]
    else:
        h_int = integrate.quad(Nu_forced_flat_plate_isoflux_lam,0.,longueur,args=(lambd,speed,nu,Pr))[0]

    return (1/longueur)*h_int
    

def top_h_simple(T_s,T_amb,theta,longueur):
    
    h_back_mean = 2.

    DT = T_s - T_amb

    T_mean = (T_s+T_amb)/2

    g = scc.g

    rho = air_rho(T_mean)
    Cp = air_c_p()
    mu = air_mu(T_mean)
    nu = air_nu(T_mean)
    lambd = air_k(T_mean)
    alpha = (lambd)/(rho*Cp)
    Pr = air_Pr()
    beta = 1/T_mean

    # if abs(DT)<=0.05:
    #     return 0.5

    if DT==0.:
        return 0.

    if DT>0:

        # Churchill and Chu for theta < 45°

        if theta>45:
            Ra_L=(g*beta*math.cos(math.pi/2-math.radians(theta))*abs(DT)*(longueur**4))/(nu*alpha)
            if Ra_L >= 1e4 and Ra_L <= 1e9:
                Nu_L = 0.68+0.67*Ra_L**(1/4)*(1+(0.492/Pr)**(9/16))**(-4/9)
                h = (lambd/longueur)*Nu_L
                return h
            elif Ra_L >= 1e9:
                Nu_L = 0.10*Ra_L**(1/3)
                h = (lambd/longueur)*Nu_L
                return h           
            else:
                print('top h simple DT>0 Ra_L',Ra_L)
                return h_back_mean
        
        # Raithby and Hollands

        elif theta<=45:
            Ra_L=(g*beta*math.sin(math.pi/2-math.radians(theta))*abs(DT)*(longueur**4))/(nu*alpha)
            if Ra_L>=1e6 and Ra_L<=1e12: # normalement c'est 1e6 to 1e11
                Nu_L = 0.14*Ra_L**(1/3)*((1+0.0107*Pr)/(1+0.01*Pr))
                h = (lambd/longueur)*Nu_L
                return h
            else:
                print('top h simple theta<=45 Ra_L',Ra_L)
                return h_back_mean

        else:
            print('theta',theta)
            return h_back_mean

    # Fujii and Imura

    if DT<0:
        if theta>=2:
            Ra_L=(g*beta*math.cos(math.pi/2-math.radians(theta))*abs(DT)*(longueur**4))/(nu*alpha)
            if Ra_L >= 1e5 and Ra_L <= 1e11:
                Nu_L = 0.68+0.67*Ra_L**(1/4)*(1+(0.492/Pr)**(9/16))**(-4/9)
                h = (lambd/longueur)*Nu_L
                return h
            else:
                print('top h simple DT<0 theta>=2 Ra_L',Ra_L)
                return h_back_mean
        else:
            print('Ra_L',Ra_L)
            return h_back_mean

    print('DT',DT)
    return h_back_mean