# -*- coding: utf-8 -*-
"""
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
number_of_customers = 50
point_X = [1,1] # make sure it is set properly
customer_points = []

##################
### Functions
##################
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

# returns angle in form of 0 to 2*pi
def angle_between(v1, v2):
    angle = atan2(v2[1]-v1[1], v2[0]-v1[0])
    if(angle<0):
        angle = 2*pi+angle
    return angle


def newline(p1, p2):
    ax = plt.gca()
    l = mlines.Line2D([p1[0], p2[0]], [p1[1], p2[1]])
    ax.add_line(l)
    return l

def newline_infinite(p1, p2): # line until plot boundaries
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


###
def create_random_points(number_of_points):
    points = []
    global max_x, max_y
### TODO add point label/value
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
    plt.plot(*zip(*points_to_plot), marker='o', ls='', color='b') #plot as (x,y) coordinates

### plot a line through point_X in the radiant direction
def plot_halfplane(point_X, vector):
    # calc norm vector for plotting halfplane
    n_vector = get_norm_vector(vector)
    u_vector = unit_vector(n_vector)
    tmp_point = [point_X[0]+u_vector[0], point_X[1]+u_vector[1]]
    newline_infinite(point_X, tmp_point)
    print "vector: ", vector
    plt.quiver(point_X[0], point_X[1], vector[0], vector[1], color=['r'], angles='xy', scale_units='xy', scale=1)
    
def get_norm_vector(vector):
    return [-vector[1],vector[0]]

def get_vector(p1,p2):
    return [p2[0]-p1[0], p2[1]-p1[1]]

###Problem_1##############################################
## find best location for Y
# find the direction of each demand point from X and sort these directions in an
# increasing order between 0 and 2pi.
def get_direction_from_X_to_point(point):
    global point_X
    angle = angle_between(point, point_X)
    vector = [point[0]-point_X[0],point[1]-point_X[1]]
    
    #####################
    return angle, vector

# calculate all angles (for sorting) & vectors from point_X
def calculate_and_sort_all_point_directions():
    global customer_points    
    dir_list = []
    for point in customer_points:
        # return angle and corresponding point
        direction, vector = get_direction_from_X_to_point(point)
        dir_list.append([direction, vector])
    return sorted(dir_list)
    

# vector pointing into halfplane
def point_is_in_half_plane(vector, point):
    def formula(x,y):
        d = vector[0]*point_X[0] + vector[1]*point_X[1]
        return (vector[0]*x + vector[1]*y -d)
    # accept on line
    tmp = formula(point[0],point[1])
    if (tmp>=0):
        return True
    else:
        return False

def get_points_in_halfplane(vector):
    result = 0
    global customer_points
    for point in customer_points:
        if point_is_in_half_plane(vector, point):
            result += 1 # MOD add buying power to points
    return result


    

### PROBLEM 1
# Find optimal location for Y, given X and customer points
def find_optimal_location_for_Y():
    global customer_points
#    plt.clf()# refresh the plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_points(customer_points) 
    global point_X
    # a bit more central point X
    point_X = get_random_point(max_x/4, max_x*3/4, max_y/4, max_y*3/4)
    ax.plot(point_X[0], point_X[1] ,marker='o', color='r', label="X")
    ax.annotate('X', xy=(point_X[0], point_X[1]))
    plt.show    
    
    point_directions = calculate_and_sort_all_point_directions()
#    print point_directions
    max_buying_power = 0
    best_half_plane_vector = [0,0]
    # TODO: Calculate in turn the buying power of all possible Y-half-planes
    for angle, vector in point_directions:
        buying_power = get_points_in_halfplane(vector)

        if (buying_power>max_buying_power):
            max_buying_power=buying_power
            best_half_plane_vector = vector
            # "animated" plot
            plot_halfplane(point_X, best_half_plane_vector) 
            print "Buying_power: ", buying_power
            plt.pause(1.5) 


#################################################################################################
### PROBLEM 2
def minimize_buying_power_of_Y():
    print "solving problem 2"
    global customer_points
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_points(customer_points) 
### Values needed for calculation
    P_min = 0
    P_max = 0
    P_0 = 0
    
### Step 1. Calculate all lines through pairs of points and calculate all Pi for
#each half-plane defined by the lines.    
    bp_list = calculate_all_lines_and_buying_power()
    if(not bp_list):
        return 0
    
### Step 2. Sort Pi in decreasing order. Set Pmirl and P,,, to the lowest and
#highest Pi I espectively.
    sorted_bps = sort_buying_powers(bp_list)
    P_min = sorted_bps[0][0] # lowest BP
#    print ("Pmin", P_min)
    P_max = sorted_bps[len(sorted_bps)-1][0] # highest BP
#    print ("Pmax", P_max)
    
    ### TODO start loop until "last P_max"(?) here
   ### if sorted_bps is empty go to step 7
    while(sorted_bps):
        sorted_bps = update_list_P_min_max(sorted_bps, P_min, P_max)
### Step 3. Set Fe to the median value in the Pi vector for all rn,in < Pi < P,,,. If
#there is IIO Pi fulfilling Pmin<Pi < I’,,,,, go to Step ‘7.
    
        P_0 = find_median_value(sorted_bps)
        
    ### Step 4. Find if there is a feasible point to all half-planes for which Vim P,.
    #This can be done by linear programming.
        if(find_feasible_point(sorted_bps, P_0)):  
    ### Step 5. if there is a feasible solution point to the problem in Step 4 then
    #minxiS( <P,. Set P,,, to P,, and go to Step 3.
            P_max = P_0
    ### Step 6. Otherwise minx(f(X)) 2 PO. Set Pmin to P,, and go to Step 3.
        else:
            P_min = P_0
    
### Step 7. A feasible point for the last P,,, is an optimal solution. The value of
#the objective function is Ptnin.
    ### TODO while schleife um step 3 bis 6
    
### remove all tuples where BP is <P_min and >P_max
def update_list_P_min_max(sorted_bps, P_min, P_max):
    updated_list = list()    
    for bp_tuple in sorted_bps:
        p_i = bp_tuple[0]
        if (p_i> P_min and p_i<P_max):
            updated_list.append(bp_tuple)
    return updated_list

            

def calculate_all_lines_and_buying_power():
    global customer_points
    
    # create list of all point pairs (non redundant - to minimize iterations)
    point_pairs = []    
    for point1 in range(len(customer_points)):
        for point2 in range(point1+1, len(customer_points)):
                point_pairs.append((customer_points[point1], customer_points[point2]))

    buying_powers_list = list() # list of tuples (bp, p1, p2)
    
    for p1, p2 in point_pairs:
### calculate buying power for both sides of the line and create two entries in the buying power list

        newline(p1, p2)
#        newline_infinite(point_pair)
###PLOTTING#############
#        newline(point_pair)
########################        
        # calculate normvector in both directions
        norm_vector_left = get_norm_vector(get_vector(p1,p2))
        norm_vector_right = get_norm_vector(get_vector(p2,p1))
        # calc bp for both sides
        bp_1 = get_points_in_halfplane(norm_vector_left)
        bp_2 = get_points_in_halfplane(norm_vector_right)

        buying_powers_list.append((bp_1,p1,p2)) 
        buying_powers_list.append((bp_2,p2,p1))

    return buying_powers_list
    
    

def sort_buying_powers(bp_list):
    ### sort in decreasing order
    return sorted(bp_list)

def find_median_value(sorted_bp_list):
    
    if(not sorted_bp_list):
        return
    # get median of the buying powers
    median = np.median([(i[0]) for i in sorted_bp_list]) # only look at the buying power in the list of tuples
    print "MEDIAN ", median
    return median
    
def find_feasible_point(sorted_bps, P_0):
    ### TODO find if there is a feasible point for all half planes, for which  P_i >= P_0

### create half plane list for which P_i >= P_0
    reduced_list = [i for i in sorted_bps if (i[0]>=P_0)]
    print "find feasible point"
    print "P0: ", P_0
    print "Sorted_bps: ", sorted_bps
    print "Reduced: ",reduced_list
    if(True):
        return 1
    else:
        return 0
    
####################################

### Main ###
def solve_problem(number):
    global customer_points
    customer_points = create_random_points(number_of_customers) # create random distribution of buying power in a certain area

    if(number==1):
        find_optimal_location_for_Y()
    elif(number==2):
        minimize_buying_power_of_Y()
solve_problem(1)
#solve_problem(2)
