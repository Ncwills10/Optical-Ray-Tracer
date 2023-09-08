# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 16:49:09 2022

@author: ncwil
"""

import numpy as np
import matplotlib.pyplot as plt
from optical_element_class import OpticalElement

class PlanoConvexLens(OpticalElement):
    '''
    This class creates a plano-convex lens, with a positive curvature resulting in the curved surface on the left, vice versa for a negative curvature. S is largest thickness of the lens.
    '''
    def __init__(self, z0 = 0, curvature = 1, refractive_indices = np.array([1, 1.5168]), s = 5):
        self.z0 = z0
        self.curvature = curvature
        self.refractive_indices = refractive_indices
        self.s = s
        height = 2 * np.sqrt(self.s * (2 * abs(1 / self.curvature) - self.s))
        self.height = height
        if curvature < 0:
            some_point_z = z0 - s
        if curvature > 0:
            some_point_z = z0 + s
        if curvature == 0:
            raise Exception("This curvature is invalid")
        some_point = np.array([0, 0, some_point_z])
        self.some_point = some_point

    def plot(self):
        radius = abs(1 / self.curvature)
        radius_array = np.linspace(radius, radius, 1000)
        x_values = np.linspace(-1 * self.height / 2, self.height / 2, 1000)
        if self.curvature < 0:
            z_values = np.linspace(self.z0 - self.s, self.z0 - self.s, 1000)
            theta_values = np.linspace(-1 * np.arccos((radius - self.s) / radius), np.arccos((radius - self.s) / radius), 1000)
            plt.plot(z_values, x_values, color = 'midnightblue')
            plt.plot(radius_array * np.cos(theta_values) + self.z0 - radius, radius_array * np.sin(theta_values), color = 'midnightblue')
        
        if self.curvature > 0:
            z_values = np.linspace(self.z0 + self.s, self.z0 + self.s, 1000)
            theta_values = np.linspace(np.pi - np.arccos((radius - self.s) / radius), np.pi + np.arccos((radius - self.s) / radius), 1000)
            plt.plot(z_values, x_values, color = 'midnightblue')
            plt.plot(radius_array * np.cos(theta_values) + self.z0 + radius, radius_array * np.sin(theta_values), color = 'midnightblue')
            
    def intercept_sphere(self, ray):
        '''
        This method finds the intercept with the spherical part of the lens.
        '''
        if self.curvature < 0:
            r = ray.position - np.array([0, 0, self.z0 + 1 / self.curvature])
            k_hat = ray.direction
            if (np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ) < 0:
                raise Exception("There is no valid intercept") 
            if (np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ) == 0:
                l = - np.dot(r, k_hat) + np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) ** 2 ))
                q = ray.position + l * k_hat
            else:                
                l_max = - np.dot(r, k_hat) + np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ))
                q = ray.position + l_max * k_hat
        if self.curvature > 0:
            r = ray.position - np.array([0, 0, self.z0 + 1 / self.curvature])
            k_hat = ray.direction
            if (np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ) < 0:
                raise Exception("There is no valid intercept") 
            if (np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ) == 0:
                l = - np.dot(r, k_hat) + np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ))
                q = ray.position + l * k_hat
            else:                
                l_min = - np.dot(r, k_hat) - np.sqrt((np.dot(r, k_hat)) ** 2 - ((np.linalg.norm(r, 2)) ** 2 - (1 / self.curvature) **2 ))
                q = ray.position + l_min * k_hat
        if np.sqrt(q[0] ** 2 + q[1] ** 2) > self.height:
            raise Exception("There is no valid intercept")
        else: 
            return q
        
    def intercept_plane(self, ray):
        '''
        This method finds the intercept with the plane part of the lens.
        '''
        d = -1 * np.dot(np.array([0, 0, 1]), self.some_point)
        if np.dot(np.array([0, 0, 1]), ray.direction) == 0:
            raise Exception("No valid intercept")
        else:
            scale_parameter = (-1 * (d + np.dot(np.array([0, 0, 1]), ray.position))) / np.dot(np.array([0, 0, 1]), ray.direction)
        q = ray.position + scale_parameter * ray.direction
        
        if np.sqrt(q[0] ** 2 + q[1] ** 2) > self.height:
            raise Exception("There is no valid intercept")
        else: 
            return q        
        
    def propagate_ray(self, incident_ray):
        '''
        This method propagates the ray through both refractive surfaces.
        '''
        if self.curvature < 0:
            intercept_sphere = self.intercept_sphere(incident_ray)
            surface_normal = (-intercept_sphere + np.array([0, 0, self.z0 + (1 / self.curvature)])) / np.linalg.norm(-intercept_sphere + np.array([0, 0, self.z0 + (1 / self.curvature)]))
            new_direction = self.refraction(incident_ray.direction, surface_normal, self.refractive_indices[1], self.refractive_indices[0])
            incident_ray.append(intercept_sphere, new_direction)
        
        if self.curvature > 0:
            intercept_sphere = self.intercept_sphere(incident_ray)
            surface_normal = (intercept_sphere - np.array([0, 0, self.z0 + (1 / self.curvature)])) / np.linalg.norm(intercept_sphere - np.array([0, 0, self.z0 + (1 / self.curvature)]))
            new_direction = self.refraction(incident_ray.direction, surface_normal, self.refractive_indices[0], self.refractive_indices[1])
            incident_ray.append(intercept_sphere, new_direction)
            intercept_plane = self.intercept_plane(incident_ray)
            surface_normal = np.array([0, 0, -1])
            new_direction = self.refraction(incident_ray.direction, surface_normal, self.refractive_indices[1], self.refractive_indices[0])

            incident_ray.append(intercept_plane, new_direction)