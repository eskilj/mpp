#!/usr/bin/env python

import sys

pixels="   1   1   1   1   1   1   1   1   1   1\n"

x=int(sys.argv[1])
y=int(sys.argv[2])

f=open("synthetic"+str(x)+"x"+str(y)+".pgm",'w')

f.write("P2\n")
f.write("# Synthetic image\n")
f.write(str(x)+" "+str(y)+"\n")
f.write("255\n")
for i in range(int(x*y/10)):
    f.write(pixels)

f.close()
