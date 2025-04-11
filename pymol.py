#pymol.py

## importing pymol packages
from pymol import cmd
import requests
import  numpy as np
import sys
import os
## Pymol Functions

#1. welcome
def Welcome():
	print("###################")
	print("Welcome to the APP")
	print("###################")
	return "Welcome"

#. load pdb

lst_pdb = []

def load(pdb):
	'''
	input -> pdb code
	action -> it hides the previous pdb files before loading new pdb files
	'''
	pdb_lst_active = cmd.get_object_list()

	while pdb_lst_active:
		pdb_remove = pdb_lst_active.pop(0)
		hide(pdb_remove)
		print(f'PDB {pdb_remove} Hidden')
	try:
		cmd.fetch(code = pdb, name = pdb)
	except:
		pass
	else:
		lst_pdb.append(pdb)
		pdb_speak = ' '.join(pdb)
		chains = cmd.get_chains(pdb)
		pdb_name = _RCSB_API(pdb)
		# speak_text(f'protein structure file {pdb_speak} {pdb_name}loaded. Which has {len(chains)} chains, {",".join(chains)}') ## to the speak function
	main() ## calling the main function for reseting everthing
	
	return "pdb_name"

## RCSB API
def _RCSB_API(pdb):
	pdb_code = "4cha"
	api_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_code}"
	headers = {"Accept": "application/json"}
	response = requests.get(api_url, headers=headers)

	# Check for successful response
	if response.status_code == 200:
		structure_data = response.json()
		pdb_title = structure_data['struct']['title']
	else:
		print(f"Error: {response.status_code}")
	return pdb_title

## Delete
def delete(pdb = 'all'):
	cmd.delete(pdb)

## zoom in and out work 
## zoom in
buffer_in_out = 0 # this is resposible for a smooth zoomin and zoomout
def zoomin():
	global buffer_in_out
	cmd.zoom(buffer = buffer_in_out - 1)
	buffer_in_out -= 1

## zoom out
def zoomout():
	global buffer_in_out
	cmd.zoom(buffer = buffer_in_out + 1)
	buffer_in_out += 1

## Rotate
#rotate x
def rotate_x(deg):
	cmd.rotate('x', deg)
#rotate y
def rotate_y(deg):
	cmd.rotate('y', deg)
#rotate z
def rotate_z(deg):
	cmd.rotate('z', deg)

def rotate_c(deg):
	cmd.rotate('z', deg)

def rotate_a(deg):
	cmd.rotate('z', -deg)

# DepthCue
def Depth_Cue():
	cmd.set("depth_cue", 1)
	cmd.set("fog_start", 0.03)

## remove water
def clean():
	cmd.remove('sol')

## set bg color
	
def bgcolor(color):
	cmd.bg_color(color)

## color
	
def chcolor(color):
	cmd.color(color, 'all')

## shadow
def shadowoff():
	cmd.set('ray_shadows', 'off')

## Reflection off

def reflection():
	cmd.set('specular', 'off')

## light
	
def light():
	cmd.set('light', [0, 0, 0]) # position of the light

## orthoscopic
	
def orthoscopic():
	cmd.set('orthoscopic', 'on')

def showonlychains(sel):   ## > sel can be a string
	'''
	sel > string
	need to do : after hiding chains, need to change the position of other chains, based on COM
	'''
	if sel.lower() == 'all':  ## if all the in the input > show all chains
		cmd.show('cartoon', 'all')

	else:
		chain_lst = cmd.get_chains()
		for chain in chain_lst:
			if chain in sel:
				cmd.show('cartoon', f'chain {chain}')
			else:
				cmd.hide('cartoon', f'chain {chain}')

# hiding chains
def hidechains(sel): 
	'''
	sel > string
	it will hide all the chains in the string
	'''
	chain_lst = cmd.get_chains()
	for chain in chain_lst:
		if chain in sel:
			cmd.hide('everything', f'chain {chain}')

# Hide 
def hide(pdb = 'all'):
	'''
	it hides the pdb file if name is given
	otherwise, it hides everything
	'''
	cmd.hide('everything', pdb)

# N-terminal and C-terminal 
def sphere(n_term_lst, c_term_lst, chain_lst):
	n_term_lst = n_term_lst.split('_')
	c_term_lst = c_term_lst.split('_')
	chain_lst = chain_lst.split('_')
	length = len(chain_lst)

	for term in range(length):
	## N-term
		cmd.select(f'sel{term}', f'chain {chain_lst[term]} and resi {n_term_lst[term]}') # and name CA
		cmd.set('sphere_scale', 0.8, f'sel{term}')
		cmd.show('spheres', f'sel{term}')

	## C-term
		cmd.select(f'sel{term}', f'chain {chain_lst[term]} and resi {c_term_lst[term]}') # and name CA
		cmd.set('sphere_scale', 0.8, f'sel{term}')
		cmd.show('spheres', f'sel{term}')
		cmd.deselect()

	## image 1
	cmd.png('3-D-printing/image1.png')
	cmd.rotate('y',90)
	cmd.png('3-D-printing/image2.png')

	## 3d printing file wrl
	cmd.save('3-D-printing/3d-image.stl')


# cartoon representation
def show(pdb):
	cmd.show('cartoon', pdb)

# Reseting to defaults
def reset():
	global buffer_in_out
	buffer_in_out = 0
	cmd.reset()
	cmd.hide('spheres', 'all')

# orient the screen
def orient(pdb = ""): 
	cmd.orient('all', state = 0)
	global buffer_in_out 
	buffer_in_out = 0

# Harshit


## RG
def calculate_rg(pdb_id):
   # Fetch the protein structure from the PDB
   cmd.fetch(pdb_id, "protein")

   # Calculate the center of mass of all atoms
   center_of_mass = cmd.centerofmass("protein")

   # Get the coordinates of all atoms
   coords = cmd.get_coords("protein")

   # Calculate the distances of each atom from the center of mass
   distances_sq = np.sum((coords - center_of_mass) ** 2, axis=1)

   # Compute the radius of gyration
   rg = np.sqrt(np.mean(distances_sq))
   print("Radius of gyration (Rg) for all atoms in the protein:", rg)

## Contact
def contacts(pdb_id):
   try:

       tmp_pdb_file = f"temp/{pdb_id}.pdb"

       # Load the protein structure from the temporary PDB file
       pdb = md.load(tmp_pdb_file)

       # Calculate contacts between all pairs of residues
       contacts = md.compute_contacts(pdb)

       # Get residue indices corresponding to contacts
       residue_pairs = contacts[1]

       # Get unique residue pairs and their counts
       unique_pairs, counts = np.unique(residue_pairs, axis=0, return_counts=True)

       # Get indices for highest, least, median, and closest contacts
       highest_index = counts.argmax()
      # least_index = counts.argmin()
       median_index = np.argsort(counts)[len(counts) // 2]
       closest_index = contacts[0].argmin()

       # Get residue pairs for highest, least, median, and closest contacts
       highest_pair = unique_pairs[highest_index]
       #least_pair = unique_pairs[least_index]
       median_pair = unique_pairs[median_index]
       closest_pair = residue_pairs[closest_index]

       # Print highest, least, median, and closest contacts
       print("Highest contact:")
       print_residue_pair_info(pdb, highest_pair)
       #print("Least contact:")
      #print_residue_pair_info(pdb, least_pair)
       print("Median contact:")
       print_residue_pair_info(pdb, median_pair)
       print("Closest contact:")
       print_residue_pair_info(pdb, closest_pair)

   except Exception as e:
       print("An error occurred:", e)

def print_residue_pair_info(pdb, pair):
   res1_name = pdb.topology.residue(pair[0]).name
   res1_number = pdb.topology.residue(pair[0]).resSeq
   res2_name = pdb.topology.residue(pair[1]).name
   res2_number = pdb.topology.residue(pair[1]).resSeq
   print(f"{res1_name} (residue {res1_number}) - {res2_name} (residue {res2_number})")


# Axis
def calculate_principal_axes(coords):
    # Perform principal component analysis
    pca = PCA(n_components=3)
    pca.fit(coords)
    # Get principal axes
    principal_axes = pca.components_
    return principal_axes

def construct_rotation_matrix(principal_axes):
    # Construct rotation matrix to align principal axes with standard axes
    rotation_matrix = np.transpose(principal_axes)
    return rotation_matrix

def align_principal_axes_with_standard_axes(pdb_id):
    # Fetch PDB file
    # cmd.fetch(pdb_id, 'protein', async_=0)
	
    # Extract protein coordinates
    coords = []
    for atom in cmd.get_model('protein').atom:
        coords.append([atom.coord[0], atom.coord[1], atom.coord[2]])
    coords = np.array(coords)

    # Calculate principal axes
    principal_axes = calculate_principal_axes(coords)

    # Construct rotation matrix
    rotation_matrix = construct_rotation_matrix(principal_axes)

    # Apply rotation to the loaded structure
    #pymol.cmd.rotate('x', 90)
    cmd.rotate('y',90)
    cmd.transform_object('protein', rotation_matrix.flatten().tolist())

    # Save aligned structure
    # cmd.save(f'aln_{pdb_id}.pdb', 'protein')
# Main function
def main():
	print("Main Function Implemented")
	Depth_Cue() # Depth cue
	clean() # remove water
	bgcolor('white') # change bg color
	chcolor('grey') # change color
	shadowoff() # off shadow
	reflection() # off reflection
	light()
	orthoscopic()



# Calling main function
main()
if __name__ == "__main__":
	main()
