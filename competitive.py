# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 15:33:32 2017

@author: Stephan
"""

from sympy import *
from sympy.geometry import *
import sympy.plotting as splt
#from sympy.plotting import plot

from random import randint, random
import matplotlib.pyplot as plt
import math
import numpy as np

# global vars
max_x = 100 # xrange
max_y = 100 # yrange
number_of_costumers = 10
point_X = [1,1] # make sure it is set properly
costumer_points = []

################## TODO imported
from math import acos
from math import sqrt
from math import pi

def length(v):
    return sqrt(v[0]**2+v[1]**2)
def dot_product(v,w):
    return v[0]*w[0]+v[1]*w[1]
def determinant(v,w):
    return v[0]*w[1]-v[1]*w[0]
def inner_angle(v,w):    
    try:
        cosx=dot_product(v,w)/(length(v)*length(w))
    except:
        cosx=0
    rad=acos(cosx) # in radians
    return rad*180/pi # returns degrees
def angle_clockwise(A, B):
    inner=inner_angle(A,B)
    det = determinant(A,B)
    if det<0: #this is a property of the det. If the det < 0 then B is clockwise of A
        return inner
    else: # if the det > 0 then A is immediately clockwise of B
        return 360-inner
##################
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    # !TODO: Problem wenn ein Vektor = Nullvektor
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    
    angle = math.atan2(v2[1]-v1[1], v2[0]-v1[0])
    # DEBUG
#    print "angle between ", v1, " and ", v2, " is ", angle
    return angle
def angle_between_old(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    dot_product =np.dot(v1_u, v2_u)
#    result = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) #TODO: check if correct
    result = np.arccos(dot_product) #TODO: check if correct
    print "angle between ", v1, " and ", v2, " is ", result
    return result
################## imported
    
    
import matplotlib.lines as mlines

def newline(p1, p2):
    ax = plt.gca()
    xmin, xmax = ax.get_xbound()

    if(p2[0] == p1[0]):
        xmin = xmax = p1[0]
        ymin, ymax = ax.get_ybound()
    else:
        ymax = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmax-p1[0])
        ymin = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmin-p1[0])

    l = mlines.Line2D([xmin,xmax], [ymin,ymax])
    ax.add_line(l)
    return l
##################

        

###### TODO: PROBLEM momentan werden sowohl x als auch y koordinaten geplottet anstatt als (x,y) Paar


### TODO: change to sympy point datatype or create point struct
def create_random_points(number_of_points):
    points = []
### TODO find out how to add point lavel/value
    for x in range(0,number_of_points): #
        new_point = [randint(0, max_x), randint(0, max_y)]
        points.append(new_point)
    return points

def plot_points(points_to_plot):
    # set scale of plot to max_x, max_y
    plt.axis([0, max_x, 0, max_y])
    plt.plot(*zip(*points_to_plot), marker='o', ls='', color='b') # damit die punkte wirklich als (x,y) koordinaten geplottet werden

### plot a line through point_X in the vector direction
def plot_halfplane(point_X, vector):
    foobar = [0,0]
    foobar[0] = point_X[0] + vector[0]
    foobar[1] = vector[1] + point_X[1]
    newline(point_X, foobar)
#def plot_line(point_X):
    #TODO: change to line
    #TODO: or plot vector and normal line to show half plane
#    plt.quiver(point_X[0], point_X[1], 1, 1)
#    plt.quiver(point_X[0], point_X[1], vector[0], vector[1])

###Problem_1##############################################
## TODO: find best location for Y
# TODO: find the direction of each demand point from X and sort these directions in an
# increasing order between 0 and 2pi.
def get_direction_from_X_to_point(point):
# TODO: REMOVE dummy code     
    # TODO: Replace with direction calc
#    angle = random()*2*math.pi
#    
#    vector[0] = point[0]-point_X[0]
#    vector[1] = point[1]-point_X[1]   
#    angle = angle_clockwise(point_X, point) ########################### TODO: returns only 270
    global point_X
    angle = angle_between(point, point_X)
    # !TODO check if correct
    #####################
    return angle
#print get_direction_from_X_to_point([0,0])

def sort_all_point_directions(points):
    dir_list = []
    for point in points:
        dir_list.append(get_direction_from_X_to_point(point))
    return sorted(dir_list)
    

#def point_is_in_half_plane(line, point):
def point_is_in_half_plane(point):
    # TODO: implement right of line check
    return True

def get_points_in_halfplane(direction, costumer_points):
    result = 0
    for point in costumer_points:
        if point_is_in_half_plane(point):
            result += 1
    # TODO: calculate normal line to direction
    #
    # TODO: plot line
    # TODO: check for all points if they're in the halfplane or not
    return result # TODO ? return points already?


    

### PROBLEM 1 ###
# Find optimal location for Y, given X and costumer points
def find_optimal_location_for_Y():
    global costumer_points
    costumer_points = create_random_points(number_of_costumers) # create random distribution of buying power in a certain area
    plt.clf()# refresh the plot
    
    plot_points(costumer_points) 
    global point_X
    point_X = [randint(0,max_x),randint(0,max_y)]
    plt.plot(point_X[0], point_X[1] ,marker='o', color='r')
    
    plt.show    
    ### pre computation test
    print "TEST", angle_between([1,0],[0,1])
    print "TEST", angle_between([-1,0],[1,0])
    
    # TODO: Calculate in turn the buying power of all possible Y-half-planes 
    point_directions = sort_all_point_directions(costumer_points)
    
    max_buying_power = 0
    best_half_plane = 0
    
    for direction in point_directions:
#        print "direction ", direction
        # TODO: get all points on half plane
        buying_power = get_points_in_halfplane(direction, costumer_points)
#        print "Buying_power: ", buying_power, "in dir:", direction
        if (buying_power>max_buying_power):
            max_buying_power=buying_power
            best_half_plane = direction
#    plot_line(point_X)
####################################
## TODO just dummy, replace through
#    normal_vector = (1,1) # TODO calculate from best_half_plane
#    plot_halfplane(point_X, normal_vector)    
####################################    
#TODO    plot_line(point_X, direction)
    # Debugprints
#    print "point direction angles", point_directions
    print "highest point direction ", np.amax(point_directions)
#    print "max buying power: ", max_buying_power
#    print "direction: ", direction
#    print "------------"
    

find_optimal_location_for_Y()


############################## sympy
## basic example
# x = symbols('x')
# p1 = plot(x*x)

##3d example
# from sympy.plotting import plot3d
# var('x y z')
# plot3d(x*y**3-y*x**3)




# p = Plot(axes='label_axes=True')
#
# c = Circle(Point(0,0), 1)
#
# t = c.tangent_line(c.random_point())
#
# p[0] = c
# p[1] = t
#
#p = Plot(axes='label_axes=True')
#
#In [7]: t = RegularPolygon(Point(0,0), 1, 5)
#
#In [8]: for i in range(len(t.sides)):
#   ....:    p[i] = t.sides[i]