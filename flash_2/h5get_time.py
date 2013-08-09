#!/usr/bin/env python

import subprocess as sub
import sys

h5file = sys.argv[1]

'''
The following set of piped subprocesses performs the following shell command 

  h5ls -l -d -S ${h5file}/simulation\ parameters | grep time | cut -d = -f 3 | cut -d ' ' -f 1

This extracts the current time of the plotfile being passed to the script as input.
All of this would be easier, of course, if bash allowed for floating point arithmetic.
'''

# Pipe the commands
p1 = sub.Popen(["h5ls","-l","-d","-S",h5file+"/simulation parameters"], stdout=sub.PIPE)
p2 = sub.Popen(["grep", "time"], stdin=p1.stdout, stdout=sub.PIPE)
p3 = sub.Popen(["cut","-d","=","-f","3"], stdin=p2.stdout, stdout=sub.PIPE)
p4 = sub.Popen(["cut","-d"," ","-f","1"], stdin=p3.stdout, stdout=sub.PIPE, stderr=sub.PIPE)

# Allow earlier open subprocesses to receive a SIGPIPE of p4 exits
p1.stdout.close()
p2.stdout.close()
p3.stdout.close()

# Get the time (output) and standard error from p4
time, errors = p4.communicate()

# Convert time from string output to floating point number
time = float(time)

# Convert from seconds to years
secyr = 3.15569e7
time = time/secyr

print time
