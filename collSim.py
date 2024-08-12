#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:09:58 2024

@author: kieran
"""

import os
os.chdir('/Users/kieran/Documents/collisionSim/')
import random

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


def ccdBorder(particleID):
     
     # Identify the borders
     borders = {'xmin':0,'xmax':parameters['GRID_SIZE'],
                'ymin':0,'ymax':parameters['GRID_SIZE']}
     # Identiy coordiante relevant to the border
     relCoord = {'xmin':'px','xmax':'px',
                 'ymin':'py','ymax':'py'}
     
     # Loop through each border
     for b in borders.keys():
          
          # Identify if a collision could happen in next frame
          diff = abs(borders[b] - particle_dict[particleID][relCoord[b]])
          
          if diff > parameters['MAX_VELOCITY']:
               continue
          
          else:
               print(diff)





def updateParticles():
     
     global particle_dict
     dt = 1/int(parameters['FPS'])

     print(particle_dict)
     # Loop through each particle
     for i in particle_dict.keys():

          # Update the x velocity
          particle_dict[i]['vx'] = updateVelocity(vel = particle_dict[i]['vx'],
                                                  acc = particle_dict[i]['a'],
                                                  delta_t = dt)
          # Update the y velocity
          particle_dict[i]['vy'] = updateVelocity(vel = particle_dict[i]['vy'],
                                                  acc = particle_dict[i]['a'],
                                                  delta_t = dt)
          
          # Update the x position
          particle_dict[i]['px'] = updatePosition(pos = particle_dict[i]['px'],
                                                  vel = particle_dict[i]['vx'],
                                                  delta_t = dt)
          # Update the y position
          particle_dict[i]['py'] = updatePosition(pos = particle_dict[i]['py'],
                                                  vel = particle_dict[i]['vy'],
                                                  delta_t = dt)


print(particle_dict)
ccdBorder(0)
