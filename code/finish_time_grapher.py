import numpy as np
from matplotlib import pyplot as plt

ALPHA = 1
BETA = 1
Q = 1
M = 0.5
c_0 = 1
a_0 = 1
s_0 =0.5
def asymps_s(x):
    return -ALPHA*Q*M*(c_0**ALPHA)*(a_0**BETA)*x + s_0


# sets time inverval between files
dt = 0.1

#creates start of volume arrays
t_list_gph = np.array([0.5])

#variable for tracking number of files
total_points = 0
while True:
    # finds file name and exits loop once finished
    N = (total_points+1)*dt 
    string = "%.1f" % float(N)
    file_name = "finish_time_data/Q="+string+".dat"
    try:
        t_list_full = np.loadtxt(file_name, delimiter=' ', usecols=(0), unpack=True)
    except FileNotFoundError:
        print("file not found"+file_name)
        break
    try:
        t_list = t_list_full[-1]
    except TypeError:
        t_list=0
    total_points+=1
    # sorts data and calculates boundary conditions
    t = t_list

    #print(boundary_con_c,boundary_con_a)
    
    # intergrates and finds total volume of cleanser and agent
    t_list_gph = np.append(t_list_gph,t)

t_list_gph = np.delete(t_list_gph, 0)
# plots and saves graphs
plt.xlabel("Q value")
plt.ylabel("time until decontamination")
x = np.arange(total_points)*dt + dt
plt.plot(x,t_list_gph,linestyle='-',linewidth=2, color='r', label='Model')
plt.legend()
plt.savefig('time until decontamination over Q.png')
plt.clf() 

