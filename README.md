# The-traveling-Salesperson-Problem

## Part 1. Using an approximation algorithm to solve TSP 
Given an undirected graph G = (V, E) we want to find a minimum length cycle that visits each
vertex once and then returns to the start vertex. In general, this problem is NP-Hard. In this exercise, we
will find a cycle that may not be of minimal length but will be no more than twice the minimum. The
solution that we will build will run in polynomial time.
Your task is to implement the following algorithm from “Introduction to Algorithms” by Cormen,
Lieserson, Rivest and Stein:
Approx-TSP-Tour(G,c)
1. Select a vertex r e V[G] to be a root vertex
2. Compute a minimum spanning tree T for G from root r using MST-Prim(G,c,r)
3. Let L be the list of vertices visited in a preorder tree walk of T
4. Return the Hamiltonian cycle H that visits the vertices in the order L

The data file that you will use is found on the course schedule and is named:
http://www.andrew.cmu.edu/user/mm6/95-771/CrimeData/CrimeLatLonXY1990.csv. This file contains a
list of crime records from 1990 in Pittsburgh. It contains three types of addresses. The first, which we will
not use, is the street address of the crime. The second is the latitude and longitude of the crime (we will use
these data only in Part 3 of this assignment. The third is the (X,Y) - coordinates of the crime. These (X, Y) - coordinates are from the State Plane Coordinate System. These data may be used to compute distances
between vertices in South Western PA using the Pythagorean theorem. In Part 1, we will be using the (X,
Y) pairs to compute the distance between each crime. To convert from feet to miles, simply multiply the
computed distances by 0.00018939.
Your program will prompt the user for two dates. If the user enters 1/1/90 and 1/1/90 as input (the
same date), your program will find an approximate TSP tour that visits the four crime locations of crimes
that occurred on that date. There were four crimes on 1/1/90.
If the user enters 1/14/90 and 1/15/90 as input, your program will find an approximate TSP tour
that visits the crime locations of those crimes that occurred between 1/14/90 and 1/15/90 inclusive. The
root of the minimum spanning tree will always be the first index.
The output of your program will include a list of the crime records processed by your TSP
approximation algorithm. It will provide a list of vertices showing a tour generated from the Approx-TSPTour algorithm shown above. It will also show the length of the tour. When you read the records into
memory, you will assign an index of 0 to the first record that the user wishes to process. That is, the first
crime that occurred on the first date of interest, will have an index of 0. For example, if the user enters
1/14/90 and 1/15/90 as input, your solution will be a list of vertices numbered from 0 and ending at 8.

## Part 2. Finding an optimal solution to TSP 
One way to find an optimal tour is to simply list every possible tour and compute the length of
each. Then, simply select the tour of minimum length. In Part 2, your task is to find an optimal tour using
this brute force approach. Note that there are |V| ! permutations of the |V| vertices. Each of these is a
different Hamiltonian cycle. But half of these are the same cycle travelled in a different direction. 

## Part 3. Displaying the output to Google Earth
In Parts 1 and 2, we computed a tour of crime scenes in two ways. The first was a tour that may or
may not have been optimal. The second was a minimum length tour. In Part 3, your task is to generate a
KML file that contains both tours (the first tour will be shown in blue and the second will be shown in red).
This KML file, when loaded into Google Earth, will show both tours. In order for both tours to be visible,
you may add a very small bit of latitude and longitude to each vertex in one of the tours. In this way, the
lines will not run exactly the same and one will not completely cover the other. Below is a screen shot of
my solution over the vertices 1 through 5. The name of the output file will be PGHCrimes.kml. The user
interaction is the same as above. That is, your program will prompt the user for two dates (which may be
the same) and then the crime data will be processed and the KML will be generated.
