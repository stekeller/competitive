# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 11:31:41 2017

@author: Stephan
"""

from random import randint
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from math import pi, atan2

# global vars
max_x = 100 # xrange
max_y = 100 # yrange
number_of_costumers = 100
point_X = [1,1] # make sure it is set properly
costumer_points = []

################## TODO imported


##################
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    # !TODO: Problem wenn ein Vektor = Nullvektor
    return vector / np.linalg.norm(vector)

# returns angle in form of 0 to 2*pi
def angle_between(v1, v2):
    angle = atan2(v2[1]-v1[1], v2[0]-v1[0])
    if(angle<0):
        angle = 2*pi+angle
# DEBUG
#    print "angle between ", v1, " and ", v2, " is ", angle
    return angle


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


### TODO: change to sympy point datatype or create point struct
def create_random_points(number_of_points):
    points = []
    global max_x, max_y
### TODO find out how to add point lavel/value
    for x in range(0,number_of_points): #
        new_point = get_random_point(0, max_x, 0, max_y)
        points.append(new_point)
    return points

# returns random x,y point between min and max
def get_random_point(min_x, max_x, min_y, max_y):
    point = [randint(min_x, max_x), randint(min_y, max_y)]    
    return point

def plot_points(points_to_plot):
    # set scale of plot to max_x, max_y
    plt.axis([0, max_x, 0, max_y])
    plt.plot(*zip(*points_to_plot), marker='o', ls='', color='b') # damit die punkte wirklich als (x,y) koordinaten geplottet werden

### plot a line through point_X in the radiant direction
def plot_halfplane(point_X, vector):
    #TODO calc norm vector
    n_vector = [-vector[1],vector[0]]
    tmp_point = [point_X[0]+n_vector[0], point_X[1]+n_vector[1]]
#    print tmp_point
    newline(point_X, tmp_point)
# OLD
#### plot a line through point_X in the radiant direction
#def plot_halfplane(point_X, vector):
#    tmp = [0,0]
#    tmp[0] = point_X[0] + vector[0]
#    tmp[1] = vector[1] + point_X[1]
#    newline(point_X, tmp)

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
    vector = [point[0]-point_X[0],point[1]-point_X[1]] # TODO ? normalvector
    
    # !TODO check if correct
    #####################
    return angle, vector
#print get_direction_from_X_to_point([0,0])


# calculate all angles from point_X
def calculate_and_sort_all_point_directions():
    global costumer_points    
    dir_list = []
    for point in costumer_points:
        # return angle and corresponding point
        direction, vector = get_direction_from_X_to_point(point)
        dir_list.append([direction, vector])
    return sorted(dir_list)
    

#def point_is_in_half_plane(line, point):
def point_is_in_half_plane(vector, point):
    def formula(x,y):
        d = vector[0]*point_X[0] + vector[1]*point_X[1]
        return (vector[0]*x + vector[1]*y -d)
    # TODO: implement line check
    # accept on line
    tmp = formula(point[0],point[1])
#    print tmp
    if (tmp>=0):
        return True
    else:
        return False

def get_points_in_halfplane(vector):
    result = 0
    global costumer_points
    for point in costumer_points:
        if point_is_in_half_plane(vector, point):
            result += 1 # eventually add buying power to points
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
#    plt.clf()# refresh the plot
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_points(costumer_points) 
    global point_X
#    point_X = [randint(0,max_x),randint(0,max_y)]
    # a bit more central point X
    point_X = get_random_point(max_x/4, max_x*3/4, max_y/4, max_y*3/4)
#    point_X = [max_x/2, max_y/2] # TODO Debug remove
#    plt.plot(point_X[0], point_X[1] ,marker='o', color='r', label="X")
    ax.plot(point_X[0], point_X[1] ,marker='o', color='r', label="X")
    ax.annotate('X', xy=(point_X[0], point_X[1]))
    plt.show    
    ### pre computation test
#    print "TEST 0=", angle_between([-1,0],[1,0]) # 0    
#    print "TEST pi/2=", angle_between([1,0],[1,1]) # pi/2
#    print "TEST pi=", angle_between([1,0],[-1,0]) # pi
#    print "TEST -pi/2=", angle_between([1,0],[1,-1]) # pi
#    print "TEST -pi/4=", angle_between([0,0],[1,-1]) # pi
#    print "TEST ~2pi=", angle_between([0,0],[1,-0.000001]) # pi
    point_directions = calculate_and_sort_all_point_directions()
#    print point_directions
    max_buying_power = 0
    best_half_plane_vector = [0,0]
    # TODO: Calculate in turn the buying power of all possible Y-half-planes
    for angle, vector in point_directions:
#        print "angle ", angle
#        print "vector ", vector
#         TODO: get all points on half plane
        buying_power = get_points_in_halfplane(vector)
#        print "Buying_power: ", buying_power, "in dir:", direction
        if (buying_power>max_buying_power):
            max_buying_power=buying_power
            best_half_plane_vector = vector
#            print "new max bp: ", max_buying_power
#            print "new best vector: ", vector
#    plot_line(point_X)
####################################
## TODO just dummy, replace through
#    best_half_plane_vector = [1,5] # TODO ! remove, testing only
    #??????????????
#    normal_dir = best_half_plane-pi/2 # TODO ? grenzfälle über 2pi und unter 0
    #??????????????
    
#    if (normal_dir<0):
#        normal_dir += 
    plot_halfplane(point_X, best_half_plane_vector)    
####################################    
#TODO    plot_line(point_X, direction)
    # Debugprints
#    print "point direction angles", point_directions
#    print "highest point direction ", np.amax(point_directions)
#    print "max buying power: ", max_buying_power
#    print "direction: ", direction
#    print "------------"
    
### START ###
find_optimal_location_for_Y()