#!/usr/bin/env python
# encoding: utf-8
"""
Emulate epidemic process using SIR model.
"""
from __future__ import print_function
from collections import defaultdict

import random#for using functionalities related to random
import math
import sys#for handling command line arguments
import operator
import yaml#for handling .yaml file
import numpy as np#for using lomax distribution

cfg = yaml.full_load(open(sys.argv[1], 'r'))#to parse a .yaml file

max_infections = cfg['maxtags']#maximum number of infections from one host
dbgprint = cfg['dbgprint']# for debugging; if false no debugging 
mincluster = cfg['mincluster']#minimum size of cluster
lambda_heal = cfg['lambda_heal']#recovery rate
minsize = cfg['minsize']#minimum size of cascade

edges = sys.argv[2] # file with edges
rate_param = sys.argv[3] # parameter to generate infections
MAXITER = int(sys.argv[4])#number of iterations

#Creating adjacency list
adjacency_list = defaultdict(set)

for item in open(edges):

    #first we remove white spaces from the beginning and end using strip() 
    #after that we split the string into a list of string after breaking the string using the \t character
    edge = item.strip().split("\t")

    #adding  second node i.e. edge[0] in adjacency list of first node i.e.edge[0] and vice versa
    #before adding since edge nodes are string so we convert them into int
    adjacency_list[int(edge[0])].add(int(edge[1]))
    adjacency_list[int(edge[1])].add(int(edge[0]))

_iter = 0

# Algorithm to infect a node:
# 1.A random node is chosen as the starting node and mark it as infected.
# 2.Then we visit every neighbour the infected node and check if the neighbour node can be infected or not.
# 3.  if the neighbour node is infected after the source node is recovered then the neighbour node cannot not be infected.
#     if the maximum number of infections from a node has reached the maximum_infections, the neighbour node cannot not infected.
# 4.If a neighbour node is infected , we check if it has already been infected or not 
# 5.the neighbour infected node is then added to infected list.
# 6.step 2 to step 5 are repeated until there is no more infected nodes. 

while True:

    start_node = random.choice(list(adjacency_list.keys()))#it randomly picks any node as starting node
    print(start_node)
    lambda_infect = np.random.pareto(float(rate_param))#sampling infection rate using lomax distribution to model a variety of cascades

    recovered = set()
    infected = list()#infected stores [time of infection,infected node,source of infection(node)] 
    infected.append( [0, start_node, -1] )#for starting node its source of infection is unknown so we choose its source of infection = -1
   
    seen_modules = set()
    cascade=[]

    new_seen_nodes = set()#this stores the list of  infected nodes for a given iteration
    new_seen_edges = set()#this stores the list of edges whose end nodes are infected
    new_nodes_counter = 0#count the number of infected nodes for a given iteration

    while len(infected)>0:

        current_time,current_node,from_node = infected[0]

        cascade.append("%d %d %f" %(from_node, current_node, current_time))#adding source of infection node,infected node,time of infection to cascade 
        
        new_nodes_counter += 1#increment the number of infected node 
        new_seen_nodes.add(current_node)#a newly infected node is added to this cascade
        
        if from_node != -1: 
            new_seen_edges.add((min(from_node,current_node),max(from_node,current_node)))

        heal_time = current_time + random.expovariate(lambda_heal)#healing time of starting node
        
        recovered.add(current_node)#infected node added to recovered node after it is healed
        
        infected = infected[1:]#the starting node is removed from the infected node

        # for degugging
        if dbgprint: 
            print(current_node,'by',from_node,'ill since',current_time,'till',heal_time)
        
        new_nodes = set(adjacency_list[current_node]) - recovered
        
        current_infected = 0#stores the number of infections from a node
    
        #visiting all neighbours of currently infected node
        for node in new_nodes:

            #calculating the infection time for every neighbour
            infection_time = current_time + random.expovariate(lambda_infect)
            
            #if infection time for a neighbour is greater than the healing time for the source, the neighbour could not be infected 
            if infection_time > heal_time:
                continue

            #if maximum number of infections is reached, the neighbour couldnt be infected    
            if current_infected > max_infections and max_infections>0: 
                continue

            current_infected += 1#counting number of infections from a node
            found = False
        
            for i,queued in enumerate(infected): #visiting all infected nodes 

                if queued[1]==node:#if the newly infected node has already been infected

                    if queued[0]>infection_time:

                        infected[i][0] = infection_time
                        infected[i][2] = current_node

                    found = True
                    break

            if not found:
                infected.append( [infection_time, node, current_node] )

        infected = sorted(infected, key=operator.itemgetter(0))#infected nodes are sorted in ascending order of their infection time
    
    if len(cascade)>= minsize:#checking if the cascade size >= minimum size of cascade

        out=["#",_iter]
        out.extend(cascade)
        
        print ("\t".join(map(str,out)))#adding space between the strings
        
        _iter += 1#incrementing the number of iteration

        if _iter >= MAXITER:#if maximum number of iterations has reached
            break
