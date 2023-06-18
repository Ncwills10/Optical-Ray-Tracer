# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 16:05:12 2022

@author: ncwil
"""

import numpy as np

class OpticalElement:    
    '''
    This is an abstract class, that contains methods that apply to all optical elements.
    '''    
    def refraction(self, incident_direction, surface_normal, n1, n2):
        '''
        This method implements Snell's Law of refraction in vector form. The surface normal is defined as facing the opposite way to the incident ray.
        '''
        if np.sin(np.arccos(np.dot(incident_direction, -1 * surface_normal))) > n2/n1:
            raise Exception("Internal reflection is occuring")
        else:
            return  (n1 / n2) * np.cross(surface_normal, np.cross(- surface_normal, incident_direction))  - surface_normal * np.sqrt(1 - np.dot(np.cross(surface_normal, incident_direction), np.cross(surface_normal, incident_direction)) * ((n1 / n2)**2))                                                                             
        
    def propagate_ray_bundle(self, ray_bundle):
        '''
        This method propagates a ray bundle through the optical element.
        '''
        ray_array_updated = np.array([])
        for a in ray_bundle.get_ray_array():
            self.propagate_ray(a)
            ray_array_updated = np.append(ray_array_updated, a)
        ray_bundle.ray_array = ray_array_updated