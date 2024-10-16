import numpy as np
from matplotlib import pyplot as plt

ALPHA = 1
BETA = 1
Q = 1
M = 0.5
c_0 = 1
a_0 = 1
s_0 =0.5



# sets time inverval between files
dt = 0.1

#creates start of volume arrays
t_list_gph = np.array([0.5])

#variable for tracking number of files
### QQQQQQQQQQQQQQQQQQQQQQQ #####
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

Q_list = np.delete(t_list_gph, 0)

### DDDDDDDDDDDDDDDDDDDDDD ###
t_list_gph = np.array([0.5])
total_points = 0
while True:
    # finds file name and exits loop once finished
    N = (total_points+1)*dt 
    string = "%.1f" % float(N)
    file_name = "finish_time_data/D_a="+string+".dat"
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

D_list = np.delete(t_list_gph, 0)

### MMMMMMMMMMMMMMMMMMMMMMMMM ###
t_list_gph = np.array([0.5])
total_points = 0
while True:
    # finds file name and exits loop once finished
    N = (total_points+1)*dt 
    string = "%.1f" % float(N)
    file_name = "finish_time_data/M="+string+".dat"
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

M_list = np.delete(t_list_gph, 0)

### CCCCCCCCCCCCCCCCCCCCC ###
t_list_gph = np.array([0.5])
total_points = 0
list_of_c_vals = [1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.5,2.6,2.8,
3,3.2,3.4,3.5,3.6,3.8,
4,4.2,4.4,4.5,4.6,4.8,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20]
for i in list_of_c_vals:
    # finds file name and exits loop once finished
    N = (total_points+10)*dt 
    string = "%.1f" % float(i)
    file_name = "finish_time_data/c_0="+string+".dat"
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

C_list = np.delete(t_list_gph, 0)


# plots and saves graphs
plt.xlabel("Q value")
plt.ylabel("time until decontamination")
x = np.arange(20)*dt + dt
x_2 = np.arange(10)*dt + dt

# plots and saves graphs
'''
plt.figure(figsize=(8,2))
ax1 = plt.subplot2grid((1, 4), (0, 0), rowspan=1)
ax2 = plt.subplot2grid((1, 4), (0, 1), colspan=1)
ax3 = plt.subplot2grid((1, 4), (0, 2), colspan=1)
ax4 = plt.subplot2grid((1, 4), (0, 3), colspan=1)

# Main plot
ax1.set_xlabel("c_0 value")
ax1.set_ylabel("Time until decontamination")
ax1.plot(list_of_c_vals,C_list,linestyle='-',linewidth=2, color='#648fff', label='Mass of Cleanser')

ax2.set_xlabel("D_a value" )
#ax2.set_ylabel("Time until decontamination")
ax2.plot(x,D_list,linestyle='-',linewidth=2, color='#fe6100', label='Mass of Cleanser')

ax3.set_xlabel("Q value")
#ax3.set_ylabel("Time until decontamination")
ax3.plot(x,Q_list,linestyle='-',linewidth=2, color='#785ef0', label='Mass of Cleanser')

ax4.set_xlabel("M value")
#ax4.set_ylabel("Time until decontamination")
ax4.plot(x_2,M_list,linestyle='-',linewidth=2, color='#dc267f', label='Mass of Cleanser')



plt.tight_layout()
plt.savefig('time until decontamination ALL.png')
plt.clf() 
'''

plt.plot(x,D_list,linestyle='-',linewidth=2, color='#fe6100', label='D_a')
plt.plot(x,Q_list,linestyle='-',linewidth=2, color='#785ef0', label='Q')
plt.plot(list_of_c_vals[:11],C_list[:11],linestyle='--',linewidth=2, color='#648fff', label='c_0')
plt.plot(x_2,M_list,linestyle='--',linewidth=2, color='#dc267f', label='M')

plt.legend()
plt.savefig('time until decontamination together.png')
plt.clf() 
