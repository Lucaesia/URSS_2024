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
dt = 0.002

#creates start of volume arrays
s_list_gph = np.array([0.5])

#variable for tracking number of files
total_points = 0
while True:
    # finds file name and exits loop once finished
    N = (total_points)*dt 
    string = "%.3f" % float(N)
    file_name = "data_output/time="+string+".dat"
    try:
        x, c, a,s_list,s_t_list = np.loadtxt(file_name, delimiter=' ', usecols=(0, 1,2,3,4), unpack=True)
    except FileNotFoundError:
        print("file not found"+file_name)
        break
    total_points+=1

    # sorts data and calculates boundary conditions
    s = s_list[0]

    #print(boundary_con_c,boundary_con_a)
    
    # intergrates and finds total volume of cleanser and agent
    s_list_gph = np.append(s_list_gph,s)

s_list_gph = np.delete(s_list_gph, 0)
# plots and saves graphs
plt.xlabel("time")
plt.ylabel("distance")
plt.ylim([0,1])
x = np.arange(total_points)*dt
plt.plot(x,s_list_gph,linestyle='-',linewidth=2, color='r', label='Model')
y = asymps_s(x)
plt.plot(x,y, linestyle='--', color='b', label='Asymptotics')
plt.legend()
plt.savefig('interface position over time.png')
plt.clf() 

