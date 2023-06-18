# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 11:50:37 2022

@author: ncwil
"""

import numpy as np
import matplotlib.pyplot as plt
from ray_class import Ray
from ray_bundle_class import RayBundle
from output_plane_class import OutputPlane
from spherical_refraction_class import SphericalRefraction
from plano_convex_lens_class import PlanoConvexLens
from biconvex_lens_class import BiconvexLens


#Task 9


spherical_refractor = SphericalRefraction(z0 = 100, curvature = 0.03, refractive_indices = np.array([1, 1.5]), aperture_radius = 1 / 0.03)
outputplane = OutputPlane(surface_normal = np.array([0, 0, 1]), z_intercept = 250)

plt.figure(1)
plt.title("Task 9 - Trajectory of a few example rays through spherical refractive surface")
plt.xlim(0, 260)
plt.ylim(-50, 50)
plt.xlabel("z axis (mm)")
plt.ylabel("x axis (mm)")
spherical_refractor.plot()
outputplane.plot()

ray1 = Ray(position = np.array([0, 0, 0]), direction = np.array([0, 0, 1]))
ray2 = Ray(position = np.array([10, 0, 0]), direction = np.array([0, 0, 1]))
ray3 = Ray(position = np.array([30, 0, 0]), direction = np.array([-0.1, 0, 1]))
ray4 = Ray(position = np.array([-15, 0, 0]), direction = np.array([0, 0, 1]))
ray5 = Ray(position = np.array([0, 0, 0]), direction = np.array([0.2, 0, 1]))


listOfRays = [ray1, ray2, ray3, ray4, ray5]


map(lambda x: spherical_refractor.propagate_ray(x), listOfRays)

map(lambda x: outputplane.propagate_ray(x), listOfRays)

map(lambda x: ray1.plot(x), listOfRays)

'''
spherical_refractor.propagate_ray(ray1)
spherical_refractor.propagate_ray(ray2)
spherical_refractor.propagate_ray(ray3)
spherical_refractor.propagate_ray(ray4)
spherical_refractor.propagate_ray(ray5)

outputplane.propagate_ray(ray1)
outputplane.propagate_ray(ray2)
outputplane.propagate_ray(ray3)
outputplane.propagate_ray(ray4)
outputplane.propagate_ray(ray5)

ray1.plot("cyan")
ray2.plot("cyan")
ray3.plot("cyan")
ray4.plot("cyan")
ray5.plot("cyan")
'''

#Task 10

near_axis_ray = Ray(position = np.array([0.1, 0, 0]), direction = np.array([0, 0, 1]))
spherical_refractor.propagate_ray(near_axis_ray)
scale_parameter = -1 * near_axis_ray.position[0] / near_axis_ray.direction[0]
paraxial_focus_point = near_axis_ray.position[2] + scale_parameter * near_axis_ray.direction[2]
print(paraxial_focus_point)
#Task 11

spherical_refractor = SphericalRefraction(z0 = 100, curvature = 0.03, refractive_indices = np.array([1, 1.5]), aperture_radius = 1 / 0.03)
outputplane = OutputPlane(surface_normal = np.array([0, 0, 1]), z_intercept = paraxial_focus_point)

plt.figure(2)
plt.title("Task 11")
plt.xlim(0, 260)
plt.ylim(-50, 50)
plt.xlabel("z axis (mm)")
plt.ylabel("x axis (mm)")
spherical_refractor.plot()
outputplane.plot()

#Task 12

plt.figure(3)
plt.title("Task 12 - Trajectories of two ray bundles of radius 10 mm")
plt.xlim(0, 260)
plt.ylim(-40, 40)
plt.xlabel("z axis (mm)")
plt.ylabel("x axis (mm)")
spherical_refractor.plot()
outputplane.plot()

bundle1 = RayBundle(3, 10)
bundle1.make_bundle()
bundle2 = RayBundle(3, 10, x_start = -51, direction_bundle = np.array([0.5, 0, 1]))
bundle2.make_bundle()

spherical_refractor.propagate_ray_bundle(bundle1)
outputplane.propagate_ray_bundle(bundle1)
spherical_refractor.propagate_ray_bundle(bundle2)
outputplane.propagate_ray_bundle(bundle2)

bundle1.plot("cyan")
bundle2.plot("dodgerblue")

#Task 13

bundle1 = RayBundle(5, 10)
bundle1.make_bundle()
bundle2 = RayBundle(5, 10, x_start = -35, direction_bundle = np.array([0.3, 0, 1]))
bundle2.make_bundle()

plt.figure(4)
plt.title("Task 13 - spot diagram at z = 0 for cyan ray bundle")
plt.xlabel("x axis (mm)")
plt.ylabel("y axis (mm)")
plt.xlim(-15,15)
bundle1.spot_diagram_plot("cyan")

plt.figure(5)
plt.title("Task 13 - spot diagram at z = 0 for blue ray bundle")
plt.xlabel("x axis (mm)")
plt.ylabel("y axis (mm)")
plt.xlim(-20,-50)
bundle2.spot_diagram_plot("dodgerblue")

spherical_refractor.propagate_ray_bundle(bundle1)
outputplane.propagate_ray_bundle(bundle1)
spherical_refractor.propagate_ray_bundle(bundle2)
outputplane.propagate_ray_bundle(bundle2)

plt.figure(6)
plt.title("Task 13 - spot diagram at paraxial focal plane for cyan ray bundle")
plt.xlabel("x axis (mm)")
plt.ylabel("y axis (mm)")
plt.xlim(-0.3, 0.3)
bundle1.spot_diagram_plot("cyan")

plt.figure(7)
plt.title("Task 13 - spot diagram paraxial focal plane for blue ray bundle")
plt.xlabel("x axis (mm)")
plt.ylabel("y axis (mm)")
plt.ylim(-1.1, 1.1)
bundle2.spot_diagram_plot("dodgerblue")

#x and y values for cyan at paraxial focus point
x_values_cyan = np.array([])
y_values_cyan = np.array([])
for a in bundle1.get_ray_array():
    x_values_cyan = np.append(x_values_cyan, a.position[0])
for b in bundle1.get_ray_array():
    y_values_cyan = np.append(y_values_cyan, b.position[1])
    
x_values_blue = np.array([])
y_values_blue = np.array([])
for c in bundle2.get_ray_array():
    x_values_blue = np.append(x_values_blue, c.position[0])
for d in bundle2.get_ray_array():
    y_values_blue = np.append(y_values_blue, d.position[1])
    
#calculating the distance from the optical axis
distance_array_cyan = np.sqrt(x_values_cyan**2 + y_values_cyan**2)
distance_array_blue = np.sqrt(x_values_blue**2 + y_values_blue**2)

#calculating the rms
rms_cyan = np.sqrt(np.sum(distance_array_cyan**2) / len(distance_array_cyan))
rms_blue = np.sqrt(np.sum(distance_array_blue**2) / len(distance_array_blue))

print(rms_cyan)
print(rms_blue)

#Task 15 - testing the plano-convex lens class

bundle3 = RayBundle(5, 10)
bundle3.make_bundle()

lens = PlanoConvexLens(curvature = 0.02, z0 = 100)
outputplane = OutputPlane(surface_normal = np.array([0, 0, 1]), z_intercept = 200)

plt.figure(8)
plt.title("Task 15 - Trajectory of a bundle of radius 10 mm through a convex-plano lens with positive curvature")
plt.xlabel("z axis (mm)")
plt.ylabel("x axis (mm)")
plt.ylim(-35, 35)
plt.xlim(50, 210)

lens.plot()
outputplane.plot()

lens.propagate_ray_bundle(bundle3)
outputplane.propagate_ray_bundle(bundle3)

bundle3.plot("cyan")

plt.figure(9)
plt.title("Task 15 - Trajectory of a single ray - refraction occurs at both surfaces")
plt.xlabel("z axis (mm)")
plt.ylabel("x axis (mm)")
plt.ylim(-40, 40)
plt.xlim(90, 110)

ray6 = Ray(position = np.array([0, 0, 90]), direction = np.array([1.1, 0, 1]))

lens.plot()
outputplane.plot()

lens.propagate_ray(ray6)
outputplane.propagate_ray(ray6)

ray6.plot("cyan")

bundle5 = RayBundle(5, 10)
bundle5.make_bundle()

lens = PlanoConvexLens(curvature = -0.02, z0 = 100)
outputplane = OutputPlane(surface_normal = np.array([0, 0, 1]), z_intercept = 200)

plt.figure(10)
plt.title("Task 15 - Trajectory of a bundle of radius 10 mm through a plano-convex lens with negative curvature")
plt.xlabel("z axis (mm)")
plt.ylabel("x axis (mm)")
plt.ylim(-35, 35)
plt.xlim(50, 210)

lens.plot()
outputplane.plot()

lens.propagate_ray_bundle(bundle5)
outputplane.propagate_ray_bundle(bundle5)

bundle5.plot("cyan")

#Task 15 - optimising the rms for positive curvature

lens = PlanoConvexLens(curvature = 0.02, z0 = 100)
ray7 = Ray(position = np.array([0.1, 0, 0]), direction = np.array([0, 0, 1]))

lens.propagate_ray(ray7)

scale_parameter = -1 * ray7.position[0] / ray7.direction[0]
paraxial_focus_point = ray7.position[2] + scale_parameter * ray7.direction[2]

outputplane = OutputPlane(surface_normal = np.array([0, 0, 1]), z_intercept = paraxial_focus_point)

radius_array = np.arange(1, 11, 1)
rms_array = np.array([])
for e in range(1, 11):
    bundle6 = RayBundle(5, e)
    bundle6.make_bundle()

    lens.propagate_ray_bundle(bundle6)
    outputplane.propagate_ray_bundle(bundle6)

    x_values = np.array([])
    y_values = np.array([])
    for f in bundle6.get_ray_array():
        x_values = np.append(x_values, f.position[0])
    for g in bundle6.get_ray_array():
        y_values = np.append(y_values, g.position[1])

    distance_array = np.sqrt(x_values**2 + y_values**2)

    rms = np.sqrt(np.sum(distance_array**2) / len(distance_array))
    rms_array = np.append(rms_array, rms)

plt.figure(11)
plt.title("Task 15: RMS radius and diffraction scale against bundle radius for positive curvature")
plt.xlabel("Bundle radius (mm)")
plt.ylabel("RMS radius (mm)")
plt.plot(radius_array, rms_array, color = "turquoise")

'''
Diffraction scale is lambda * focal length / aperture diameter  - therefore when 
optimising by making the beam radius smaller,
we cannot go below the diffraction limit
'''

plt.xlim(1, 11)
plt.ylim(0, 0.08)
x_arbitrary = np.linspace(0, 15, 1000)
diffraction_scale = 588e-6 * paraxial_focus_point / (2 * x_arbitrary)
plt.plot(x_arbitrary, diffraction_scale, color = "darkblue")
plt.legend(["RMS radius", "Diffraction scale"])



#Task 15 - optimising the rms for negative curvature

lens = PlanoConvexLens(curvature = -0.02, z0 = 100)
ray8 = Ray(position = np.array([0.1, 0, 0]), direction = np.array([0, 0, 1]))

lens.propagate_ray(ray8)

scale_parameter = -1 * ray8.position[0] / ray8.direction[0]
paraxial_focus_point = ray8.position[2] + scale_parameter * ray8.direction[2]

outputplane = OutputPlane(surface_normal = np.array([0, 0, 1]), z_intercept = paraxial_focus_point)

radius_array = np.arange(1, 11, 1)
rms_array = np.array([])
for h in range(1, 11):
    bundle6 = RayBundle(5, h)
    bundle6.make_bundle()

    lens.propagate_ray_bundle(bundle6)
    outputplane.propagate_ray_bundle(bundle6)

    x_values = np.array([])
    y_values = np.array([])
    for i in bundle6.get_ray_array():
        x_values = np.append(x_values, i.position[0])
    for j in bundle6.get_ray_array():
        y_values = np.append(y_values, j.position[1])

    distance_array = np.sqrt(x_values**2 + y_values**2)

    rms = np.sqrt(np.sum(distance_array**2) / len(distance_array))
    rms_array = np.append(rms_array, rms)

plt.figure(12)
plt.title("Task 15: RMS radius and diffraction scale against bundle radius for negative curvature")
plt.xlabel("Bundle radius (mm)")
plt.ylabel("RMS radius (mm)")
plt.plot(radius_array, rms_array, color = "turquoise")

plt.xlim(1, 11)
plt.ylim(0, 0.33)
x_arbitrary = np.linspace(0, 15, 1000)
diffraction_scale = 588e-6 * paraxial_focus_point / (2 * x_arbitrary)
plt.plot(x_arbitrary, diffraction_scale, color = "darkblue")
plt.legend(["RMS radius", "Diffraction scale"])

#%%

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
from ray_class import Ray
from ray_bundle_class import RayBundle
from output_plane_class import OutputPlane
from spherical_refraction_class import SphericalRefraction
from plano_convex_lens_class import PlanoConvexLens
from biconvex_lens_class import BiconvexLens

'''
Lens optimisation - we know that the diffraction scale is lambda * focal length / aperture diameter  - therefore when 
optimising the lens by making the beam radius smaller,
we cannot go below the diffraction limit

Could not figure out how to optimise this
'''

def rms_func(parameters_array):
    curv1 = parameters_array[0]
    curv2 = parameters_array[1]
    lens = BiconvexLens(curvature_l = curv1, curvature_r = curv2, z_joint = 50, refractive_indices = np.array([1, 1.5168]), height = 43.6)
    near_axis_ray = Ray(position = np.array([0.1, 0, 0]), direction = np.array([0, 0, 1]))
    lens.propagate_ray(near_axis_ray)
    scale_parameter = -1 * near_axis_ray.position[0] / near_axis_ray.direction[0]
    paraxial_focus_point = near_axis_ray.position[2] + scale_parameter * near_axis_ray.direction[2]
    outputplane = OutputPlane(surface_normal = np.array([0, 0, 1]), z_intercept = paraxial_focus_point)
    bundle = RayBundle(5, 10)
    bundle.make_bundle()
    lens.propagate_ray_bundle(bundle)
    outputplane.propagate_ray_bundle(bundle)
    #plt.figure(13)
    #plt.title("Trajectory of a ray bundle through a biconcave lens")
    #lens.plot()
    #bundle.plot("cyan")
    x_values = np.array([])
    y_values = np.array([])
    for a in bundle.get_ray_array():
        x_values = np.append(x_values, a.position[0])
    for b in bundle.get_ray_array():
        y_values = np.append(y_values, b.position[1])
    distance_array = np.sqrt(x_values**2 + y_values**2)
    rms = np.sqrt(np.sum(distance_array**2) / len(distance_array))
    return rms

initial_guesses = [0.01, 0.01]

rms_func(np.array([0.02,0.03]))