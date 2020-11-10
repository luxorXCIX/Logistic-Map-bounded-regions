# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 15:33:39 2020

@author: Louis Coffin

This code is intended to generate a rough sketch of the set of complex numbers for which iterating rx(1-x) for a given r
by performing a Monte Carlo simulation on complex numbers. If some iterate leaves a particular bounding circle, then it and
all of its iterates will be put in a list of points that tend to infinity, infBasin, and otherwise that point will be assumed
to remain bounded and put in a list of points that remain bounded, finBasin.

Whether or not a point x is in the bounding circle for r is determined by the magnitude of r/2 - x*(r^2 - 1)/2. If this
magnitude is greater than the maximum of 2 and the magnitude of r(r - 2)/4, then iterating x with the logistic map will give a sequence
which tends to infinity. The converse is not necessarily true. The bounding circle is the set of points where these two
magnitudes are equal.

The reason this works is that the iterates of x under rx(1-x) can be determined by iterating z^2 + c with c = r(r - 2)/4
and initial value r/2 - x(r^2 - 1)/2 and then performing a particular transformation. The benefit of doing this is that there
is a convenient bounding circle for z^2 + c. If |z| > max(|c|, 2), then the iterates of z will tend towards infinity.
"""

import matplotlib.pyplot as plt
import random

#Takes in complex coordinates x and r and returns the complex numbers xres and c which makes the iteration x^2 + c 
#equivalent to the iteration rx(1-x).
#
def toMandel(x, r):
    xres = r * (0.5 - x)
    c = r*(r-2)/4
    return (xres, c)

#Checks if a complex coordinate is in the bounding circle.

def checkInBound(x, r):
    mand = toMandel(x, r)
    return abs(mand[0]) < max(2, abs(mand[1] * (2 - mand[1]) / 4))

#This will be the r used in rx(1-x)
r = 3

k = max(2, abs(r * (2-r) / 4))
bound = k / abs(r) + 1/2


#infBasin is the set of points which will tend towards infinity, finBasin is the set which haven't yet
infBasin = []
finBasin = []

#Search depth and number of repetitions of the Monte Carlo simulation. One thing to note is that while a deeper search
#will theoretically lead to a more accurate map, the chaotic nature of the logistic map means that the small errors
#that build up over the iterations can lead to completely different values. Therefore, after a certain point, increasing
#depth is likely to decrease the accuracy of the sketch.
depth = 15
numTests = 20000


for i in range(numTests):
    checklist = []
    x, y = 0, 0
    
    #print(infBasin.count((x,y)), finBasin.count((x,y)), checkInBound(complex(x,y), r))
    
    #Will generate coordinates inside of the bounding circle which have not already been checked. The random numbers
    #are generated within a box containing the bounding circle, and then regenerated if they are not within the circle.
    while infBasin.count((x,y)) > 0 or finBasin.count((x,y)) > 0 or not checkInBound(complex(x,y), r):
        x = random.uniform(-bound, bound)
        y = random.uniform(-bound, bound)
    
    z = complex(x,y)
    checklist.append((x,y))
    
    #Performs the iteration, checking at every step if the value of z is in the bounding circle.
    for d in range(depth):
        z = r * z * (1 - z)
        checklist.append((z.real,z.imag))
        if not checkInBound(z, r):
            break
    
    #Adds the current list of checked points to the appropriate basin.
    if checkInBound(z, r):
        finBasin.append(checklist[0])
    else:
        infBasin.extend(checklist)
        
#Plots the points in the finite basin. Color and plotting style can be changed by changing
#the last input in plt.plot, see the section on format strings here: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.plot.html
for (x,y) in finBasin:
    plt.plot(x,y, 'g,')
plt.axis('equal') #Makes the axes equal, i.e. circles look like circles, not ovals
plt.title('r = ' + str(r), fontsize=10)
plt.show()
    
