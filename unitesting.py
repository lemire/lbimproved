#! /usr/bin/env python
import sys
sys.path.append("..")
sys.path.append(".")
import time
import dtw
from data.syntheticdata import * 

ok = True

print "Please wait a few minutes..."
n=10
querydataset = [randomwalk(n)  for i in range(100)]
trainingdataset = [randomwalk(n)  for i in range(1000)]

reducdim = 10
for reducdim in [10,5]:
 for c in [0,1,2]:
	rtree=dtw.TimeSeriesTree("mytmpfile.bin",c,reducdim)
	for x in trainingdataset:
	  rtree.add(x)
	for i in range(len(trainingdataset)):
	  readcopy = rtree.readTimeSeries(i)
	  for j in range(n):
	    if(abs(readcopy[j]-trainingdataset[i][j])>0.00001):
	      print "possible bug!"
	      break
	for x in querydataset:
	  treesolution = rtree.getNearestNeighborCost(x,rtree.NAIVE)
	  golden = dtw.LB_Keogh(x,c)
	  golden2 = dtw.LB_Improved(x,c)
	  golden3 = dtw.DimReducedLB_Keogh(x,c,reducdim)
	  for z in trainingdataset:
	    golden.test(z)
	    golden2.test(z)
	    golden3.test(z)
	  if((treesolution <> golden.getLowestCost()) or (treesolution <> golden2.getLowestCost()) or (treesolution <> golden3.getLowestCost())) : 
	    print treesolution, golden.getLowestCost(),golden2.getLowestCost(),golden3.getLowestCost()
	    print "BUG!!!"
	    ok = False
	    break
	rtree.close()
rtree=dtw.TimeSeriesTree("mytmpfile.bin")
for x in trainingdataset:
	  rtree.add(x)
for i in range(len(trainingdataset)):
	  readcopy = rtree.readTimeSeries(i)
	  for j in range(n):
	    if(abs(readcopy[j]-trainingdataset[i][j])>0.00001):
	      print "possible bug!"
	      break
rtree.close()

for c in [0,1,2]:
  print c
  lb1  = [dtw.LB_Keogh(x,c)  for x in querydataset]
  lb2  = [dtw.LB_Improved(x,c)  for x in querydataset]
  for x in trainingdataset:
    for i in range(len(lb1)):
     lb1[i].test(x)
     lb2[i].test(x)
     if(lb1[i].getLowestCost()<>lb2[i].getLowestCost()):
      print "Possible BUG = ", lb1[i].getLowestCost(),lb2[i].getLowestCost()
      ok = False
      break
if(ok):
  print "your code checks out"
else:
  print "you appear to have a bug"