#!/bin/python

# Calculate the velocity autocorrelation functions

import numpy as np
import math
import sys

""" Read steps of frame from Input File
     @natoms: number of atoms
     @steps: number of steps
     @fp: file descriptor of this file
     @frame: array to store frame information
"""
def readFrame(natoms, steps, fp, frames):
  for my_step in range(0,steps):
  #  for iline in range(0,9):
  #    line = fp.readline().split()
  
    for i in range(0, natoms):
      line = fp.readline().split()
      frames[my_step][i] = [float(x) for x in line[1:5]]


def correlation(steps, natoms, frames):
  shape = frames.shape
  length = shape[0]
  print(shape[0])
  if steps > shape[0]:
    sys.exit()

  cc = np.zeros(steps)
  cp = np.zeros(steps)
  cv = np.zeros(steps)
  for dt in range(0, steps, 5):
    for i in range(0, natoms):
      vv = frames[0:length-dt, i, 3:4]
      ff = frames[dt:length, i, 0:1]
      vp = frames[dt:length, i, 1:2]
      vr = frames[dt:length, i, 3:4]
      cc[dt] += sum(sum(np.multiply(vv, ff))) / float(vv.shape[0])
      cp[dt] += sum(sum(np.multiply(vv, vp))) / float(vv.shape[0])
      cv[dt] += sum(sum(np.multiply(vv, vr))) / float(vv.shape[0])
    cc[dt] = cc[dt] / float(natoms)
    cp[dt] = cp[dt] / float(natoms)
    cv[dt] = cv[dt] / float(natoms)
    print(dt)

  np.savetxt("cfv.dat", cc, delimiter="\t")
  np.savetxt("cvpv.dat", cp, delimiter="\t")
  np.savetxt("cvv.dat", cv, delimiter="\t")

"""Unit Test"""
if __name__ == "__main__":
  args = sys.argv
  input_file = args[1]
  steps = int(args[2])
  fp = open(input_file)
  natoms = 125
  nframes = 300000

  frames = np.zeros((nframes, natoms, 4))
  readFrame(natoms, nframes, fp, frames)
  correlation(steps, natoms, frames)
