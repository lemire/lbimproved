

def readFlat(filename, delimiter):
  f = open(filename)
  ans = []
  for line in f:
   ans.append(map(lambda x:float(x), filter(lambda x:len(x)>0,line.split(delimiter))))
  return ans

