# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 16:32:02 2022

@author: ncwil
"""

import numpy as np
from optical_element_class import OpticalElement
import matplotlib.pyplot as plt

class SphericalRefraction(OpticalElement):
    '''
    This class creates a spherical refractor, defined by its z-intercept, its curvature, its aperture radius and the refractive indices.
    '''
    def __init__(self, z0 = 0, curvature = 1, refractive_indices = np.array([0,0]), aperture_radius = 1):
        self.z0 = z0
        self.curvature = curvature
        self.refractive_indices = refractive_indices
        self.aperture_radius = aperture_radius     
        if curvature == 0:
            raise Exception("This is a spherical refractor...")
            
    def intercept(self, ray):
        '''
        k_hat is the unit vector in the direction of the ray. Q is the intersection of the ray and the spherical refractor. 
        This method calculates the intercept of the ray with the refractor and returns it.
        '''
        r = ray.position - np.array([0, 0, self.z0 + 1 / self.curvature])
        k_hat = ray.direction
        if (np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ) < 0:
            raise Exception("There is no valid intercept") 
        if (np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ) == 0:
            l = - np.dot(r, k_hat) + np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ))
            q = ray.position + l * k_hat
            return q
        else:                
            l_min = - np.dot(r, k_hat) - np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ))
            q = ray.position + l_min * k_hat
            return q

    def propagate_ray(self, incident_ray):
        '''
        This method propagates the ray to the refractor, updating its position and direction.
        '''
        intercept = self.intercept(incident_ray)
        surface_normal = (intercept - np.array([0, 0, self.z0 + (1 / self.curvature)])) / np.linalg.norm(intercept - np.array([0, 0, self.z0 + (1 / self.curvature)]))
        new_direction = self.refraction(incident_ray.direction, surface_normal , self.refractive_indices[0], self.refractive_indices[1])
        incident_ray.append(intercept, new_direction)

                   
    def plot(self):
        r = 1 / self.curvature
        theta = np.linspace(np.pi / 2, 3 * np.pi / 2, 1000)
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        plt.plot(z + self.z0 + r, x, color = 'midnightblue')