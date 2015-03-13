import re
import time

def readyahoofinance(filename):
  """ as the name suggests"""
  regex = "([0-9]+-[A-Za-z]+-[0-9]+),([0-9]*\.[0-9]*),([0-9]*\.[0-9]*),([0-9]*\.[0-9]*),([0-9]*\.[0-9]*),([0-9]+)"
  # date format seems to be month/day/year
  data = []
  in_file =open(filename,'r')
  for line in in_file:
    #print line
    matchedcontent = re.match(regex,line)
    #print matchedcontent
    if(matchedcontent<> None):
      #print matchedcontent.group(1)
      try:
        date =round(time.mktime( time.strptime(matchedcontent.group(1),"%d-%b-%y"))/(3600*24)) #days since epoch
      except OverflowError:
        continue
      #print date
      openval = float(matchedcontent.group(2))
      #print open
      high = float(matchedcontent.group(3))
      #print high
      low = float(matchedcontent.group(4))
      #print low
      close = float(matchedcontent.group(5))
      #print close
      volume = int(matchedcontent.group(6))
      #print volume
      data.append((date,low))#close))
  print "read ", len(data), " data points!"
  in_file.close()
  assert len(data)>0
  data.sort()
  #data.reverse()
  return data

def readmoneycentral(filename):
  """ This reads data from Microsoft Money Central web site"""
  regex = "([0-9]+/[0-9]+/[0-9]+),([0-9]*\.[0-9]*),([0-9]*\.[0-9]*),([0-9]*\.[0-9]*),([0-9]*\.[0-9]*),([0-9]+)"
  # date format seems to be month/day/year
  data = []
  in_file =open(filename)
  for line in in_file:
    #print line
    matchedcontent = re.match(regex,line)
    #print matchedcontent
    if(matchedcontent<> None):
      date =round(time.mktime( time.strptime(matchedcontent.group(1),"%m/%d/%Y"))/(3600*24)) #days since epoch
      #print date
      openval = float(matchedcontent.group(2))
      #print open
      high = float(matchedcontent.group(3))
      #print high
      low = float(matchedcontent.group(4))
      #print low
      close = float(matchedcontent.group(5))
      #print close
      volume = int(matchedcontent.group(6))
      #print volume
      data.append((date,low))#close))
  print "read ", len(data), " data points!"
  in_file.close()
  assert len(data)>0
  data.sort()
  #data.reverse()
  return data

