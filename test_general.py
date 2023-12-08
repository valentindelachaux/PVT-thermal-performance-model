
import numpy as np
import scipy.optimize as sco
import networkx as nx

class Node:
    def __init__(self, id, thermostat, temperature=0., volume=None):
        self.id = id
        self.temperature = temperature
        self.volume = volume
        self.thermostat = thermostat
        self.S = lambda T : 0.

    def source(self, S): # S is a function of T vector
        self.S = S

    def equation(self, f): # f is a function of T vector - equation = 0
        self.f = f

class Edge:
    def __init__(self, N1, N2, f):
        """
        Args:
            N1: Node 1
            N2: Node 2
            f: function which takes two temperatures and a graph and returns the corresponding heat transfer
            
        Returns:
            None"""
        
        self.N1 = N1
        self.N2 = N2
        self.N1.id = N1.id
        self.N2.id = N2.id
        self.f = f

    def conduction(self,R,S):
        self.f = lambda T : S*((T[self.N1.id] - T[self.N2.id])/R)

    def convection_or_radiation(self,h,S):
        """
        Args:
            h: heat transfer function which takes temperatures and returns the corresponding heat transfer coefficient
        """
        self.f = lambda T : h(T)*S*(T[self.N1.id] - T[self.N2.id])

    def convection_and_radiation(self,h_conv,S_conv,h_rad,S_rad):
        """
        Args:
            h: heat transfer function which takes two temperatures and returns the corresponding heat transfer coefficient
        """
        self.f = lambda T : h_conv(T)*S_conv*(T[self.N1.id] - T[self.N2.id]) + h_rad(T)*S_rad*(T[self.N1.id] - T[self.N2.id])

    def heat_transfer(self,T):
        self.ht = self.f(T)

class Graph:
    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges
    
    def energy_balance(self,T):
        balance = np.zeros_like(T)

        for node in self.nodes:
            if node.thermostat:
                balance[node.id] = T[node.id] - node.temperature

        for node in self.nodes:
            if not node.thermostat:
                q_net = 0
                flag=0
                for edge in self.edges:
                    if edge.N1 == node:
                        q_net -= edge.f(T)  # Heat transfer from N1 to N2
                        flag=1
                    elif edge.N2 == node:
                        q_net += edge.f(T)  # Heat transfer from N2 to N1
                        flag=1
                if flag==1:
                    balance[node.id] = q_net + node.S(T)  # Energy balance equation for the node
                else:
                    balance[node.id] = node.f(T)  # Energy balance equation for the node
        return balance

    def objective(self,T):
        balance = self.energy_balance(T)
        return np.sum(balance)  # We want to minimize the sum of the squares of the balances
    
    def solve(self,T_guess):

        # Solve the system of equations using fsolve
        T_final = sco.fsolve(self.energy_balance, T_guess)

        # Temperature bounds (example)
        # T_bounds = [(x-10, x+10) for x in T_guess]  # Adjust your bounds accordingly

        # T = T_guess

        # eps = 1e-3
        # residuals = self.energy_balance(T)

        # compt=0
        # while max(residuals) > eps:
        #     print("compt",compt,"max(residuals)",max(residuals))
        #     # randomly choose a node number
        #     node_id = np.random.randint(0, len(self.nodes))

        #     res = sco.minimize(lambda T : self.energy_balance(T)[node_id], T, method='L-BFGS-B', bounds=T_bounds)

        #     if res.success:
        #         T = res.x
        #         residuals = self.energy_balance(T)
        #     else:
        #         pass
        #     compt+=1

    
        # Update the temperatures of the nodes
        for i, node in enumerate(self.nodes):
            if not node.thermostat:
                node.temperature = T_final[i]



        # # Minimize the objective function
        # res = sco.minimize(self.objective, T_guess, method='L-BFGS-B', bounds=T_bounds)

        # # Check if the solver succeeded
        # if res.success:
        #     T_final = res.x
        #     # Update the temperatures of the nodes
        #     for i, node in enumerate(self.nodes):
        #         if not node.thermostat:
        #             node.temperature = T_final[i]
        #     for i, edge in enumerate(self.edges):
        #         edge.heat_transfer(np.array([node.temperature for node in self.nodes]))

        # else:
        #     print('The solver did not converge:', res.message)

    def show_temp(self):
        for node in self.nodes:
            print("Node", node.id, "temperature:", node.temperature)

    def show_ht(self):
        for edge in self.edges:
            print("Edge", edge.N1.id, edge.N2.id, "heat transfer:", edge.ht)
