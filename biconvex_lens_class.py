# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 13:49:53 2022

@author: ncwil
"""

import numpy as np
import matplotlib.pyplot as plt
from optical_element_class import OpticalElement

class BiconvexLens(OpticalElement):
    '''
    This class creates a biconvex lens. The curvatures must be positive, and its plane surface is at z_joint, and its height is defined.
    '''
    def __init__(self, curvature_l = 0.02, curvature_r = 0.03, z_joint = 0, refractive_indices = np.array([1,1.5]), height = 1):
        self.curvature_l = curvature_l
        self.curvature_r = curvature_r
        self.z_joint = z_joint
        self.refractive_indices = refractive_indices
        self.height = height
        if curvature_l < 0:
            raise Exception("Curvature must be positive")
        if curvature_r < 0:
            raise Exception("Curvature must be positive")
            
    def plot(self):
        z_values = np.linspace(self.z_joint, self.z_joint, 1000)
        x_values = np.linspace(-self.height / 2, self.height / 2, 1000)
        plt.plot(z_values, x_values, color = "midnightblue")
        thetal_values = np.linspace(np.pi - np.arcsin(self.height * self.curvature_l / 2), np.pi + np.arcsin(self.height * self.curvature_l / 2), 1000)
        radiusl_values = np.linspace(1 / self.curvature_l, 1 / self.curvature_l, 1000)
        thetar_values = np.linspace(-np.arcsin(self.height * self.curvature_r / 2), np.arcsin(self.height * self.curvature_r / 2), 1000)
        radiusr_values = np.linspace(1 / self.curvature_r, 1 / self.curvature_r, 1000)
        plt.plot(radiusl_values * np.cos(thetal_values) + self.z_joint + np.sqrt((1/self.curvature_l)**2 - (self.height**2/4)), radiusl_values * np.sin(thetal_values), color = "midnightblue")
        plt.plot(radiusr_values * np.cos(thetar_values) + self.z_joint - np.sqrt((1/self.curvature_r)**2 - (self.height**2/4)), radiusr_values * np.sin(thetar_values), color = "midnightblue")       

    def intercept_l(self, ray):
        '''
        Finds the intercept with the left curve.
        '''
        r = ray.position - np.array([0, 0, self.z_joint + np.sqrt((1/self.curvature_l)**2 - (self.height**2/4))])
        k_hat = ray.direction
        if (np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature_l) **2 ) < 0:
            raise Exception("There is no valid intercept") 
        if (np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature_l) **2 ) == 0:
            l = - np.dot(r, k_hat) + np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature_l) **2 ))
            q = ray.position + l * k_hat
            return q
        else:                
            l_min = - np.dot(r, k_hat) - np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature_l) **2 ))
            q = ray.position + l_min * k_hat
            return q
    
    def intercept_r(self, ray):
        '''
        Finds the intercept with the right curve.
        '''
        r = ray.position - np.array([0, 0, self.z_joint - np.sqrt((1/self.curvature_r)**2 - (self.height**2/4))])
        k_hat = ray.direction
        if (np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature_r) **2 ) < 0:
            raise Exception("There is no valid intercept") 
        if (np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature_r) **2 ) == 0:
            l = - np.dot(r, k_hat) + np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature_r) **2 ))
            q = ray.position + l * k_hat
            return q
        else:                
            l_max = - np.dot(r, k_hat) + np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature_r) **2 ))
            l_min = - np.dot(r, k_hat) - np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature_r) **2 ))
            q_max = ray.position + l_max * k_hat
            q_min = ray.position + l_min * k_hat
            if q_max[2] > 0:
                return q_max
            if q_min[2] > 0:
                return q_min
            
    def propagate_ray(self, incident_ray):
        intercept_l = self.intercept_l(incident_ray)
        surface_normal = (intercept_l - np.array([0, 0, self.z_joint + np.sqrt((1/self.curvature_l)**2 - (self.height**2/4))])) / np.linalg.norm(intercept_l - np.array([0, 0, self.z_joint + np.sqrt((1/self.curvature_l)**2 - (self.height**2/4))]))
        new_direction = self.refraction(incident_ray.direction, surface_normal , self.refractive_indices[0], self.refractive_indices[1])
        incident_ray.append(intercept_l, new_direction)
        
        intercept_r = self.intercept_r(incident_ray)
        surface_normal = (np.array([0, 0, self.z_joint - np.sqrt((1/self.curvature_r)**2 - (self.height**2/4))]) - intercept_r) / np.linalg.norm(np.array([0, 0, self.z_joint - np.sqrt((1/self.curvature_r)**2 - (self.height**2/4))]) - intercept_r)
        new_direction = self.refraction(incident_ray.direction, surface_normal , self.refractive_indices[1], self.refractive_indices[0])
        incident_ray.append(intercept_r, new_direction)