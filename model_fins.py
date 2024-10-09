import sys
import math
import copy
import pandas as pd
import numpy as np
import heat_transfer as bht
import ht

import model_ht as mht
import model_temperatures as mtemp

from CoolProp.CoolProp import PropsSI

import model_transfers as mtr

import sympy as sp

# Define symbolic variables
x, L_fin, lambd, Ac, k, Bi, T_ext, T_0 = sp.symbols('x L_fin \lambda Ac k Bi T_ext T_0')

def fin_analytical(bc):

    # Define the temperature profile solution
    A = sp.Symbol('A')
    B = sp.Symbol('B')

    # Temperature function T(x)
    T = T_ext + A * sp.cosh(sp.sqrt(2 * Bi) * x/lambd) + B * sp.sinh(sp.sqrt(2 * Bi) * x/lambd)

    # Boundary condition at x = 0 (T(0) = T_0)
    boundary_condition_0 = sp.Eq(T.subs(x, 0), T_0)

    # Solve for A
    A_solution = sp.solve(boundary_condition_0, A)[0]

    # Substitute the value of A in T(x)
    T = T.subs(A, A_solution)

    dTdx = sp.diff(T, x)

    if bc == 'free_end':
        boundary_condition_L_fin = sp.Eq( dTdx.subs(x, L_fin), -Bi/lambd * (T.subs(x, L_fin) - T_ext) )
        B_solution = sp.solve(boundary_condition_L_fin, B)[0]
    elif bc == 'adiabatic':
        boundary_condition_L_fin = sp.Eq( dTdx.subs(x, L_fin), 0)
        B_solution = sp.solve(boundary_condition_L_fin, B)[0]
    else:
        raise ValueError('Boundary condition not recognized')

    T = T.subs(B, B_solution)
    dTdx = dTdx.subs(B, B_solution)

    expression = (T-T_ext) / (T_0-T_ext)
    expression = sp.simplify(expression)

    if bc == "adiabatic":
        expression = sp.simplify(expression * sp.cosh(sp.sqrt(2 * Bi) * L_fin/lambd)) / sp.cosh(sp.sqrt(2 * Bi) * L_fin/lambd)

    gamma = sp.diff(expression,x).subs(x, 0)
    gamma = - k * Ac * sp.simplify(gamma)

    Qdot =  gamma * (T_0 - T_ext)

    return T, dTdx, expression, gamma, Qdot

# gamma_free_end = fin_analytical('free_end')[3]
# gamma_adia = fin_analytical('adiabatic')[3]

def gamma_fin(bc, L_fin_, lambd_, Ac_, k_, Bi_):

    if bc == "free_end":
        # return float(gamma_free_end.subs({L_fin: L_fin_, lambd: lambd_, Ac: Ac_, k: k_, Bi: Bi_}).evalf())
    
        return (Ac_ * k_ * (np.sqrt(2) * Bi_**(3/2) * np.cosh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) + 
               2 * Bi_ * np.sinh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_))) / (
               lambd_ * (np.sqrt(2) * np.sqrt(Bi_) * np.cosh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) + 
               Bi_ * np.sinh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_)))
    
    elif bc == "adiabatic":
        # return float(gamma_adia.subs({L_fin: L_fin_, lambd: lambd_, Ac: Ac_, k:k_, Bi: Bi_}).evalf())

        return (np.sqrt(2) * Ac_ * np.sqrt(Bi_) * k_ * np.tanh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_)) / lambd_
    else:
        raise ValueError('Boundary condition not recognized')

# T_free_end = fin_analytical('free_end')[0]
# T_adia = fin_analytical('adiabatic')[0]

def T_fin_mean(bc, L_fin_, lambd_, Ac_, k_, Bi_, T_ext_, T_0_):

    # if bc == "free_end":
    #     T_to_integrate = T_free_end
    # elif bc == "adiabatic":
    #     T_to_integrate = T_adia
    # else:
    #     raise ValueError('Boundary condition not recognized')
    
    # T_mean = (sp.integrate(T_to_integrate, (x, 0, L_fin)) / L_fin).subs({L_fin: L_fin_, lambd: lambd_, Ac: Ac_, k: k_, Bi: Bi_, T_ext: T_ext_, T_0: T_0_})
    
    if bc == "free_end":
        T_mean = (L_fin_ * T_ext_ + (
                                    (np.sqrt(2) * lambd_ * (T_0_ - T_ext_) * np.sinh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) / (2 * np.sqrt(Bi_))) -
                                    np.sqrt(2) * lambd_ * (
                                        np.sqrt(2) * np.sqrt(Bi_) * T_0_ * np.sinh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) -
                                        np.sqrt(2) * np.sqrt(Bi_) * T_ext_ * np.sinh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) +
                                        Bi_ * T_0_ * np.cosh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) -
                                        Bi_ * T_ext_ * np.cosh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_)
                                    ) * np.cosh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) / (
                                        2 * np.sqrt(Bi_) * (np.sqrt(2) * np.sqrt(Bi_) * np.cosh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) +
                                        Bi_ * np.sinh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_))
                                    ) + 
                                    np.sqrt(2) * lambd_ * (
                                        np.sqrt(2) * np.sqrt(Bi_) * T_0_ * np.sinh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) -
                                        np.sqrt(2) * np.sqrt(Bi_) * T_ext_ * np.sinh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) +
                                        Bi_ * T_0_ * np.cosh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) -
                                        Bi_ * T_ext_ * np.cosh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_)
                                    ) / (
                                        2 * np.sqrt(Bi_) * (np.sqrt(2) * np.sqrt(Bi_) * np.cosh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) +
                                        Bi_ * np.sinh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_))
                                    )
                                )
                ) / L_fin_   
    
    elif bc == "adiabatic":
        T_mean = (L_fin_ * T_ext_ + (
                                    (np.sqrt(2) * lambd_ * (T_0_ - T_ext_) * np.sinh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) / (2 * np.sqrt(Bi_))) -
                                    (np.sqrt(2) * lambd_ * (T_0_ - T_ext_) * np.cosh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) *
                                    np.tanh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) / (2 * np.sqrt(Bi_))) +
                                    (np.sqrt(2) * lambd_ * (T_0_ - T_ext_) * np.tanh(np.sqrt(2) * np.sqrt(Bi_) * L_fin_ / lambd_) / (2 * np.sqrt(Bi_)))
                                    )
                ) / L_fin_
    else:
        raise ValueError('Boundary condition not recognized')
    
    return T_mean