# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 15:06:07 2022

@author: ncwil
"""

import numpy as np
import matplotlib.pyplot as plt

class Ray():
    '''
    This is a class to create rays. Points is an array containing all previous 
    points on the line, trajectories is an array of all previous directions. 
    Both are numpy arrays of all the points as numpy arrays. 
    Any initial direction can be inputted, and it automatically normalises it.
    '''
    
    def __init__(self, position = np.array([0, 0, 0]), direction 
                 = np.array([1, 1, 1]), points 
                 = np.empty((1, 3)), directions 
                 = np.empty((1, 3))):
        mag = np.linalg.norm(direction, 2)
        direction = direction / mag
        
        self.position   = position
        self.directions = directions
        self.points     = points
        self.direction  = direction
        
        
    def __repr__(self):
        '''
        This method allows you to print the ray's current position and direction
        '''
        return "position: %r, direction: %r" % (self.position, self.direction)
    
    def k(self):  
        return self.direction
    
    def p(self):
        return self.position
    
    def append(self, p, k):
        '''
        This method adds a new current position and direction, as well as 
        appending the old ones to the storage lists.
        '''
        self.points     = np.vstack((self.points, self.position))
        self.directions = np.vstack((self.directions, self.direction))
        self.position   = p
        self.direction  = k
        
    def vertices(self):
        '''
        Returns all points of ray: therefore all previous points stored in 'points' plus the current position
        '''
        return np.vstack((self.points, self.position))
    
    def trajectories(self):
        '''
        Returns all previous directions of ray
        '''
        return np.vstack((self.directions, self.direction))
    
    def plot(self, colour):
        plt.plot(self.vertices()[1:,2], self.vertices()[1:, 0], color 
                 = str(colour))