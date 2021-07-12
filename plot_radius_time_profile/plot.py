from locale import atof
import matplotlib.pyplot as plt
import numpy as np

r=[]
t=[]

def read_file(filename,data):
    lines = []
    with open(filename) as f:
        lines = f.readlines()

    count = 0
    for line in lines:
        count += 1
        data.append(atof(line))

def read_group(filenames,r_g,t_g):
    for i in range(len(filenames)):
        r = []
        read_file(filenames[i],r)
        t =  np.linspace(0,len(r)*0.47/4.,len(r))
        r_g.append(r)
        t_g.append(t)

filenames = []
read_group(filenames,r,t)


plt.plot(t,r)
plt.xlabel('t (Gyr)')
plt.ylabel('r (kpc)')
plt.legend(loc = 'lower left')

plt.show()