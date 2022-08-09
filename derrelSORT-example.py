#!/bin/python3

import numpy as np
import derrelSORT as dS
import time as t
import random as rdm
import sys

try:
  array_len = int(sys.argv[1])
except IndexError:
  array_len = 100000000

# Create an array with array_len elements
print(50*'-')
print("Creating array of", array_len, "random integers.")
t0 = t.time()
x = np.asfortranarray(np.array([round(100000*rdm.random(),0) 
                      for i in range(array_len)]).astype(np.int32))
t1 = t.time()
print('Creation time:', round(t1-t0, 2), 'seconds')


# Sort the array using python's sorted function
print(50*'-')
print("Sorting the array with Python's sorted function.")
t0 = t.time()
y = sorted(x)
t1 = t.time()
print('Sorting time:', round(t1-t0, 2), 'seconds')

# Sort the array using derrelSORT
print(50*'-')
print("Sorting the array with derrelSORT.")
t0 = t.time()
y = dS.isort(x)
t1 = t.time()
print('Sorting time:', round(t1-t0, 2), 'seconds')
print(50*'-')
