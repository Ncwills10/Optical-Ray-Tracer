# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 15:55:13 2022

@author: ncwil
"""

import numpy as np
from optical_element_class import OpticalElement
import matplotlib.pyplot as plt

class OutputPlane(OpticalElement):
    def __init__(self, surface_normal = np.array([0, 0, 0]), z_intercept = 0):
        
        '''
        This class creates an output plane, where the rays terminate.
        It is defined by its surface normal, and also its z-intercept. 
        d is part of the equation of a plane
        '''
        
        self.surface_normal = surface_normal
        some_point = np.array([0, 0, z_intercept])
        self.some_point = some_point
        
    def intercept(self, ray):
        '''
        scale_parameter is the lambda in the vector equation for the ray
        q is the intersection point. This method finds the intercept with the plane.
        '''
        d = -1 * np.dot(self.surface_normal, self.some_point)
        
        if np.dot(self.surface_normal, ray.direction) == 0:
            raise Exception("No valid intercept") 
        else:
            scale_parameter = (-1 * (d + np.dot(self.surface_normal, ray.position))) / np.dot(self.surface_normal, ray.direction)
        q = ray.position + scale_parameter * ray.direction
        return q
    
    def propagate_ray(self, incident_ray):
        '''
        This method propagates the ray to the output plane.
        '''
        
        intercept = self.intercept(incident_ray)
        incident_ray.append(intercept, np.array([0, 0, 0]))

    def plot(self):
        plt.vlines(self.some_point,-300,300, color = 'midnightblue')