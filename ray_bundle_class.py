# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 16:29:10 2022

@author: ncwil
"""

import numpy as np
from ray_class import Ray
import matplotlib.pyplot as plt

class RayBundle():
    def __init__(self, N, R, z_start = 0, x_start = 0, direction_bundle = np.array([0, 0, 1]), ray_array = np.array([])):
        '''
        This class is a bundle of evenly spaced rays, defined by the number of concentric rings in the cross section, the radius, and its starting point in the x and z dimensions.
        '''
        
        self.N = N
        self.R = R
        self.z_start = z_start
        self.ray_array = ray_array
        self.x_start = x_start
        
        mag = np.linalg.norm(direction_bundle,2)
        direction_bundle = direction_bundle / mag
        self.direction_bundle = direction_bundle
        
    def make_bundle(self):
        '''
        This method creates the start positions of the rays and adds them to the attribute ray_array, which contains all the rays in the bundle as objects.
        '''
        
        T = self.N * (3 * self.N + 3) + 1
        print("Number of points = " + str(T))      
        r_array_values = np.linspace(0, self.R, self.N + 1)
        r_array_values = r_array_values[1:]
        n_array = np.array([])
        for a in range(0, self.N):
            n = 6 * (a + 1)
            n_array = np.append(n_array, n)
        n_array = n_array.astype(int)
        theta_array_values = 2 * np.pi / n_array
        theta_array = np.array([])
        r_array = np.array([])
        for index, b in enumerate(n_array):
            for c in range(0, b):
                r_array = np.append(r_array, r_array_values[index])
                theta_array = np.append(theta_array, c * theta_array_values[index])
        theta_array = np.insert(theta_array, 0, 0)
        r_array = np.insert(r_array, 0, 0)
        ray_array_updated = np.array([])
        for d in range(0, len(r_array)):
            component_ray = Ray(position = np.array([r_array[d] * np.cos(theta_array[d]) + self.x_start, r_array[d] * np.sin(theta_array[d]), self.z_start]), direction = self.direction_bundle)
            ray_array_updated = np.append(ray_array_updated, component_ray) 
        self.ray_array = ray_array_updated
        
    def get_ray_array(self):
        return self.ray_array
    
    def plot(self, colour):
        for a in self.ray_array:
            a.plot(colour)
    
    def spot_diagram_plot(self, colour):
        '''
        This method plots the spot diagram, as y against x
        '''
        x_values = np.array([])
        y_values = np.array([])
        for b in self.get_ray_array():
            x_values = np.append(x_values, b.position[0])
            y_values = np.append(y_values, b.position[1])
        plt.scatter(x_values, y_values, color = colour)