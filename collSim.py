#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:09:58 2024

@author: kieran
"""

import os
os.chdir('/Users/kieran/Documents/collisionSim/')
import random
import pandas as pd
import statistics

# Function to read in the parameter
def FindParameters(filename):
    
    dictionary = {}
    
    # Open the file
    file = open(filename)
    # Parse the contents
    for ln in file.readlines():

        ln = ln.split()
        
        if ln[0] == 'RADII_RANGE':
            new_list = list()
            
            # Remove the rubbish
            ln[1] = ln[1].strip('][')
            ln[1] = ln[1].split(',')
            # Convert to numeric
            for n in ln[1]:
                new_list.append(int(n))
            ln[1] = new_list
        else:
             
            # Add parameters to dictionary
            ln[1] = int(ln[1])
        dictionary[ln[0]] = ln[1]
        
    return(dictionary)


parameters = FindParameters('parameters.txt')

# Function to randomly generate a certain number of coordinates
def randCoord():

    # Isolate the variables from the parameters that is needed
    grid_size = int(parameters['GRID_SIZE'])
    num_part = int(parameters['NUM_PARTICLES'])
    all_coord = list()
    
    # Loop through every x value
    for x in range(grid_size):
        # Loop though every y value
        for y in range(grid_size):
                # Add the coordinate to the list
                all_coord.append([x,y])
    
    # Select the required number of coordinates (without replacement)
    list_of_coord = random.sample(all_coord, num_part)

    return(list_of_coord)
    
# Function to randomly choose a velocity
def randVelocity():
     
     # Identity the absobulte value of the maximum velocity
     max_vel = int(parameters['MAX_VELOCITY'])
     # Randomly generate a flot
     n = random.random()
     # Generate the velocity from the float
     vel_range = max_vel*2
     v = n*vel_range
     v = v-max_vel

     return(v)

def genAcceleration():
     
     return(0)

def randRadius():
     lower = parameters['RADII_RANGE'][0]
     upper = parameters['RADII_RANGE'][1]
     r = random.randint(lower,upper)
     #print(r)
     return(r)

# Function to generate a population of particles
def GenParticles():
    
    dictionary = {}
    num_part = int(parameters['NUM_PARTICLES'])

    # Generate a list of unique coordinates
    selected_coords = randCoord()

    # Loop through each particle
    for i in range(num_part):
        
        dictionary[i] = {
             # Define the coordinates of the particles
             'px':selected_coords[i][0],'py':selected_coords[i][1],
             # Define the velocities of the particles
             'vx':randVelocity(),'vy':randVelocity(),
             # Define the acceleration of the particle
             'a':genAcceleration(),
             # Define the radius
             'r':randRadius()
             }
    
    return(dictionary)

particle_dict = GenParticles()


def updateVelocity(vel,acc,delta_t):
     
     vel1 = vel + (acc*delta_t)

     return(vel1)

def updatePosition(pos,vel,delta_t):
     
     pos1 = pos + (vel * delta_t)

     return(pos1)

def borderCollisionDetection(x,y,vel_x,vel_y,r):

     xmin = 0
     xmax = parameters['GRID_SIZE']
     ymin = 0
     ymax = parameters['GRID_SIZE']

     # West border collision
     if x < (xmin + r):
          print("Collided with west border wall")

          # Correct the position
          distance_over = r-x
          x = r+distance_over
          # Flip the velocity
          vel_x = vel_x*-1



     # East border collision
     if x > (xmax-r):
          print("Collided with east border wall")
          
          # Correct the position
          distance_over = x-(xmax-r)
          x = xmax-r-distance_over
          # Flip the velocity
          vel_x = vel_x*-1

     # South border collision
     if y < (ymin+r):
          print("Collided with south border wall")

          # Correct the position
          distance_over = r-y
          y = r+distance_over
          # Flip the velocity
          vel_y = vel_y*-1
     
     # North border collision
     if y > (ymax-r):
          print("Collided with north border wall")

          # Correct the position
          distance_over = y-(ymax-r)
          y = ymax-r-distance_over
          # Flip the velocity
          vel_y = vel_y*-1
     

     return(x,y,vel_x,vel_y)

# Funciton to detect if there is a collision between two particles
def particleCollisionDetection(particle_ID1,particle_ID2):
    
    global particle_dict
    collision_bool = False
    
    # Identify first particles coordinates
    px1 = particle_dict[particle_ID1]['px']
    py1 = particle_dict[particle_ID1]['py']
    
    # Identify second particles coordinates
    px2 = particle_dict[particle_ID2]['px']
    py2 = particle_dict[particle_ID2]['py']
    
    # Identify particles' radii
    r1 = particle_dict[particle_ID1]['r']
    r2 = particle_dict[particle_ID2]['r']
    
    # Calculate their distance squares
    particle_dist_sq = ((px1-px2)**2) + ((py1-py2)**2)
    
    # Calculate the square of the minimum distance
    min_dist_sq = (r1+r2)**2
    
    # Determine collision
    if particle_dist_sq < min_dist_sq:
        collision_bool = True
        
    return(collision_bool)


# Function to extract one feature of a particle from the global dictionary
def extractParticleFeature(feature):
    
    # Initialise rubbish
    global particle_dict
    feature_dict = {}
    
    # Sequentially add the feature to a 1D dictionary
    for i in particle_dict.keys():
        
        feature_dict[i] = particle_dict[i][feature]
    
    return(feature_dict)

# Function that determines if two ranges overlap in any way
def detectOverlap(rangeA,rangeB):
    
    # I hate this
    if ((rangeA[0] <= rangeB[1]) and (rangeB[0] <= rangeA[1])):
        
        overlap_bool = True
    else:
        overlap_bool = False
        
    return(overlap_bool)
    


#dict(sorted(people.items(), key=lambda item: item[1]))
def makeRanges(df):
    
    # Import global data
    global particle_dict
    dictionary = {}
    
    # Loop through particles
    for i in df.keys():
        
        # Calculate the upper and lower bounds of their positions
        upper = df[i] + particle_dict[i]['r']
        lower = df[i] - particle_dict[i]['r']
        
        # Collate into a dictionary
        dictionary[i] = [lower,upper]
    
    return(dictionary)
        

# Function to detect candidates for true collision detection
def detectCandidates():
    
    # Import the particle metadara into the function
    global particle_dict
    dictionary = {}
    
    # Get coordinates across each axis
    x_coords = extractParticleFeature(feature = 'px')
    y_coords = extractParticleFeature(feature = 'py')
    
    # Calculate variance across each axis
    x_var = statistics.variance(list(x_coords.values()))
    y_var = statistics.variance(list(y_coords.values()))
    
    # Decide which axis to use
    if x_var >= y_var:
        data = x_coords
        axis = 'px'
    else:
        data = y_coords
        axis = 'py'
    
    # Order the particles in the dictionary by their axis coordinate
    data = dict(sorted(data.items(), key=lambda item: item[1]))
    
    # Add in their intervals along the axis
    range_data = makeRanges(df = data)
    
    # Loop through the particles
    for i in data.keys():
        
        # Save identify candiate collisons
        
        # Isolate particle of interest's full range
        range_i = range_data[i]
        
        # Loop through potential candidates
        for j in range_data.keys():
            
            # Skip if it is the same particle
            if i == j:
                continue
            
            # Isolate the potential candidate'f full range
            range_j = range_data[j]
            
            # If there is overlap of the ranges
            if detectOverlap(range_i,range_j):
                
                # Add the identity candidaite to a list
                # Declare list if necessary
                if i not in dictionary.keys():
                    dictionary[i] = [j]
                else:
                    dictionary[i].append(j)
    return(dictionary)
            


def updateParticles(particle_ID):
     
     global particle_dict
     global potential_collisions
     global resolved_this_frame

     dt = 1/int(parameters['FPS'])
     a = particle_dict[particle_ID]['a']


     # Isolate the current velocities
     vx0 = particle_dict[particle_ID]['vx']
     vy0 = particle_dict[particle_ID]['vy']
     # Isolate the current positions
     px0 = particle_dict[particle_ID]['px']
     py0 = particle_dict[particle_ID]['py']

     # Update the velocities
     vx1 = updateVelocity(vel = vx0,
                          acc = a,
                          delta_t = dt)
     vy1 = updateVelocity(vel = vy0,
                          acc = a,
                          delta_t = dt)
     # Update the positions
     px1 = updatePosition(pos = px0,
                          vel = vx1,
                          delta_t = dt)
     py1 = updatePosition(pos = py0,
                          vel = vy1,
                          delta_t = dt)
     

     # Insert collision detection and correction here
     px1,py1,vx1,vy1 = borderCollisionDetection(x = px1,y = py1,
                                                vel_x = vx1, vel_y = vy1,
                                                r = particle_dict[particle_ID]['r'])
         
     
     # Update the global dictionary
     particle_dict[particle_ID]['vx'] = vx1
     particle_dict[particle_ID]['vy'] = vy1
     particle_dict[particle_ID]['px'] = px1
     particle_dict[particle_ID]['py'] = py1
     
     # Only do particle-particle collision if there are candidate particles
     if particle_ID in potential_collisions.keys():
         
         # Add this particle to the list of resolved particles
         # These will be skipped in future ture particle collision detection tests
         resolved_this_frame.append(particle_ID)
         
         # Isolate the list of potentials
         potentials = potential_collisions[particle_ID]
         
         # Loop through the candidate particles
         for candidate in potentials:
             # Skip those who have previously been resolved
             if candidate in resolved_this_frame:
                 continue
             else:
                 
                 # If the particles are truly colliding, correct them
                 if particleCollisionDetection(particle_ID,candidate):
                     
                     print("Particle-Particle collision!")
                     
                         
             
             

# Function to export a dictionary as a '.csv'
def getArray(dictionary):
    data = pd.DataFrame.from_dict(dictionary)
    return(data.transpose())



num_frames = parameters['FPS']*parameters['TIME']


for i in range(num_frames):
    
     print(i)
     
     # Particle collision detection
     potential_collisions = detectCandidates()
     resolved_this_frame = list()
     print(potential_collisions)
     
     for p in particle_dict.keys():
         
         updateParticles(p)
     
     filename = (str(i)+'_frame_data.csv')

     getArray(particle_dict).to_csv('output_data/' + filename)
















