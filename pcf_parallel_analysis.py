"""
THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import numpy as np

def handle_mol(mol_ind,stride_len,n_strides,pcf_c,n_bins,n_bead_types,cg_sites,box_length,n_molec,n_frames_reduced,dim,d_r,lock):

	"""
	Calculates differences of centers-of-geometry between molecules
	within the same simulation frame and sorts into the pcf 
	histogram vector "pcf_c".

	Input
	----------
	mol_ind : :obj:`int`
		Index of the first molecule to calculate 
		distances to all not previously visited 
		other molecules.
	stride_len : :obj:`int`
		Increment of mol_ind after a single sweep.
	n_strides : :obj:`int`
		Number of sweeps. 
	pcf_c : (n_bins,) multiprocessing.Array of :obj:`ctypes.c_double`
		Aray to hold the pcf histogram.
	n_bins : :obj:`int`
		Number of pcf histogram bins.
	n_bead_types : :obj:`int`
		Number of bead types of the goarse-grained molecule.
	cg_sites: :obj:`MDAnalysis.core.groups.AtomGroup`
		Container of positions of all molecules.
	box_length : :obj:`float`
		Simulation box length.
	n_molec : :obj:`int`
		Number of molecules in the trajectory file.
	n_frames_reduces : :obj:`int`
		Number of frames of the trajectory to be used
		to calculate the pcf.
	dim : :obj:`int`
		Dimension.
	d_r : :obj:`floar`
		Bin-width of the pcf histogram.
	lock : :obj:`multiprocessing.Lock()
		Lock object to prevent loss of sumands through simultanteous trials to write into
		shared memory`

	Output
	----------
	None


	"""


	pcf = np.frombuffer(pcf_c.get_obj()).reshape(n_bins+1)

	# Each of these loops correlates one molecule (mol_ind)with all other molecules.
	for stride_ind in range(n_strides):
		
		# Index of current molecule to correlate all other molecules
		if(stride_ind!=0):
			mol_ind = mol_ind + stride_len
		else:
			continue

		# Extract coordinates of current molecule
		coords_pbc=np.zeros((n_bead_types,dim))
		for type_ind in range(n_bead_types):
			bead_ind = mol_ind*n_bead_types + type_ind
			coords_pbc[type_ind]=cg_sites[bead_ind].position
		
		# Define first bead as reference and enforce PBC such that all sites are in the same periodic box
		for type_ind in range(n_bead_types-1):
			for dim_ind in range(dim):
				if(coords_pbc[0][dim_ind]-coords_pbc[type_ind+1][dim_ind]<box_length/2.0):
					coords_pbc[type_ind+1][dim_ind]-=box_length
				if(coords_pbc[0][dim_ind]-coords_pbc[type_ind+1][dim_ind]>=box_length/2.0):
					coords_pbc[type_ind+1][dim_ind]+=box_length

		# Center-of-geometry
		cog=np.zeros((dim))
		for dim_ind in range(dim):
			for type_ind in range(n_bead_types):
				cog[dim_ind] += (1.0/n_bead_types)*coords_pbc[type_ind][dim_ind]
		
		# Loop over all other molecules in the box but only those who have not been visited before
		for mol_ind_other_h in range(n_molec-(mol_ind+1)):
			
			coords_pbc=np.zeros((n_bead_types,dim))
			mol_ind_other=mol_ind_other_h+mol_ind+1 

			# Extract coordinates of other molecule
			for type_ind in range(n_bead_types):
				bead_ind = mol_ind_other*n_bead_types + type_ind
				coords_pbc[type_ind]=cg_sites[bead_ind].position
			
			
			# Define first bead as reference and enforce PBC such that all sites are in the same periodic box
			for type_ind in range(n_bead_types-1):
				for dim_ind in range(dim):
				
					if(coords_pbc[0][dim_ind]-coords_pbc[type_ind+1][dim_ind]<box_length/2.0):
						coords_pbc[type_ind+1][dim_ind]-=box_length
					if(coords_pbc[0][dim_ind]-coords_pbc[type_ind+1][dim_ind]>=box_length/2.0):
						coords_pbc[type_ind+1][dim_ind]+=box_length
			
			# Center-of-geometry
			cog_other=np.zeros((dim))
			for dim_ind in range(dim):
				for type_ind in range(n_bead_types):
					cog_other[dim_ind] += (1.0/n_bead_types)*coords_pbc[type_ind][dim_ind]

			dx_cog =  abs(cog_other[0] - cog[0])
			dy_cog =  abs(cog_other[1] - cog[1])
			dz_cog =  abs(cog_other[2] - cog[2])
			
			# Enforce PBC
			if(dx_cog > box_length/2.0):
				dx_cog = box_length - dx_cog
			if(dy_cog > box_length/2.0):
				dy_cog = box_length - dy_cog 
			if(dz_cog > box_length/2.0):
				dz_cog = box_length - dz_cog 

			dist_cog = np.sqrt(dx_cog**2.0 + dy_cog**2.0 + dz_cog**2.0)
			bin_ind = int( dist_cog / d_r  + 0.5 )
			
			# Increment corresponding bin
			lock.acquire()
			pcf[bin_ind] += (2.0/(n_frames_reduced))*(box_length**3/(n_molec))*(1.0/(n_molec*4.0*np.pi/3.0))
			lock.release()
			

	return
