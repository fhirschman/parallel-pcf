"""
THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

"""
A short script to generate pair correlation functions from simulation 
trajectories accelerated by using the python multiprocessing package 
(https://docs.python.org/2/library/multiprocessing.html).  

The sample trajectory "traj.gsd" consists of 161 frames of a Brownian Dynamics 
simulation of 343 coarse-grained molecules, where each of them consists of six 
equally sized spheres. Coordinates are in units of bead-diameters. There are steric 
(repulsive) interactions between molecules and attractive interactions between 
beads within the same molecule, whose centers-of-geometry are considered for the 
spatial correlation. The simulation was conducted with periodic boundary conditions,
implying that all distances x bigger than BOX_LENGTH/2 need to be recalculated to 
(BOX_LENGTH - x).  

MDAnalysis (https://www.mdanalysis.org/) is used for processing the trajectory.

Parallelization is implemented by distributing the pairwise distance calculations
with (N_MOLEC*(N_MOLEC-1))/2 terms for each frame to all available CPU threads. 
"""

import MDAnalysis as mda
import numpy as np
import multiprocessing as mp
import ctypes as c
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
from pcf_parallel_analysis import handle_mol

# Bin-width in units of bead-diameters
D_R=0.01

# Fraction of frames in the file to use 
FRACTION=0.5

# Trajectory path
T_PATH = "traj.gsd"

# Generate universe, extract/calculate several quantities
u = mda.Universe(T_PATH)
BOX_LENGTH=u.dimensions[0]
CG_SITES = u.select_atoms('all')
N_BEADS=len(CG_SITES)
N_BEAD_TYPES=len(set(CG_SITES.types))
N_MOLEC=N_BEADS//N_BEAD_TYPES
N_FRAMES=len(u.trajectory)
N_FRAMES_REDUCED=int(N_FRAMES*FRACTION)
START_FRAME=N_FRAMES-N_FRAMES_REDUCED
N_BINS=int((BOX_LENGTH*2.)/(D_R)+1)
R_VEC=np.linspace(0.,BOX_LENGTH*2.,num=N_BINS+1)
DIM=3
 
# Declare result array 
pcf_c=mp.Array(c.c_double,N_BINS+1)

# Use all available CPU cores
N_PROC = mp.cpu_count()
print("Using " + str(N_PROC) + " threads.\n")

# Number of molecules for each process
N_STRIDES=N_MOLEC//N_PROC
# Number of molecules for each remaining process
N_MOLEC_REMAIN=N_MOLEC%N_PROC


u = mda.Universe(T_PATH)
for ts in u.trajectory:
	if(ts.frame<START_FRAME):
		continue
	print("Current frame: ",ts.frame)

	# Do N_PROC*N_STRIDES molecules
	CG_SITES = u.select_atoms('all')
	processes=[]
	for proc_ind in range(N_PROC):
		if(len(processes)< N_PROC):
			current_proc = mp.Process(target=handle_mol,args=(proc_ind,N_PROC,N_STRIDES,pcf_c,N_BINS,N_BEAD_TYPES,CG_SITES,BOX_LENGTH,N_MOLEC,N_FRAMES_REDUCED,DIM,D_R))
			current_proc.start()
			processes.append(current_proc)
		else:
			for proc_ind in range(len(processes)):
				if(not(processes[proc_ind].is_alive())):
					processes.pop(proc_ind) 
					break
			time.sleep(1)

	# Do all leftover molecules and wait until all are finished
	n_procs_remain=0
	while( (len(processes) != 0) or (n_procs_remain != N_MOLEC_REMAIN)):
		time.sleep(1)
		for proc_ind in range(len(processes)):
			if(not(processes[proc_ind].is_alive())):
				processes.pop(proc_ind) 
				if(n_procs_remain != N_MOLEC_REMAIN):
					current_proc = mp.Process(target=handle_mol,args=(N_MOLEC-(N_MOLEC_REMAIN-n_procs_remain),0,1,pcf_c,N_BINS,N_BEAD_TYPES,CG_SITES,BOX_LENGTH,N_MOLEC,N_FRAMES_REDUCED,DIM,D_R))
					current_proc.start()
					processes.append(current_proc)
					n_procs_remain += 1 
				break
pcf = np.frombuffer(pcf_c.get_obj()).reshape((N_BINS+1))

# Normalize
for bin_ind in range(N_BINS+1):
	pcf[bin_ind] = (pcf[bin_ind])*(1.0/(( (R_VEC[bin_ind]+D_R)**3.0  - R_VEC[bin_ind]**3.0 )))

# Save result
with open("pcf.dat", "w") as f:
	f.write("# Distace [sigma] \t PCF \n")
	np.savetxt(f, np.column_stack([R_VEC,pcf]), delimiter="\t")


# Plot result
PKF_FIGSIZE=[5,4]
FONTSIZE = 12
LINEWIDTH = 1 
plt.rc('font', family='stixgeneral', serif='Times New Roman')
mpl.rc('xtick', labelsize=FONTSIZE)
mpl.rc('ytick', labelsize=FONTSIZE)
plt.rcParams['mathtext.fontset'] = 'stix'
fig,ax = plt.subplots(figsize=[5,4])
ax.grid()
ax.set_xlabel('$\sigma$',fontsize=FONTSIZE)
ax.set_ylabel('PCF',fontsize=FONTSIZE)
ax.ticklabel_format(style='sci', axis='x',scilimits=(0,0), useMathText=True)
pkf_data = np.loadtxt('pcf.dat',delimiter ='\t',dtype=float,comments=[" ","#","@"])
ax.plot(np.hsplit(pkf_data,2)[0],np.hsplit(pkf_data,2)[1],c="red",linewidth=LINEWIDTH,label="PCF COG")
plt.xlim([-0.05,5])
plt.ylim([-0.1,2.0])
ax.legend(fontsize=FONTSIZE)
plt.tight_layout()
fig.savefig("pcf.png",bbox_inches = 'tight')
