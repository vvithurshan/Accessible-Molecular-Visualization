from flask import Flask, render_template
from flask import request
import requests
import os
import shutil
import MDAnalysis as mda
import mdtraj as md
from MDAnalysis.analysis import contacts
import numpy as np
from collections import defaultdict
import google.generativeai as genai
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

## making folder for 3d printing and rendering
# if not os.path.existsp('3-D-printing') and os.path.isdir('3-D-printing'):
#     os.mkdir('3-D-printing')

# making temp folder for working directory
if os.path.exists('temp') and os.path.isdir('temp'):
    shutil.rmtree('temp', ignore_errors=True)
    os.mkdir('temp')
else:
    os.mkdir('temp')

## making folder for contact map
if not os.path.exists('contactmap') and os.path.isdir('contactmap'):
    os.mkdir('contactmap')

## flask app
app = Flask(__name__)

## global variable for pdb file
pdb_code = None

# function for downloading pdb from RSCB
def download_pdb(pdb_id):
    # Define the URL for the PDB file
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    # Send a GET request to download the PDB file
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Save the PDB file
        with open(f"{pdb_id}.pdb", "wb") as pdb_file:
            pdb_file.write(response.content)
        #     ## move to the temp file
        original_path = os.path.join(os.getcwd(),f'{pdb_id}.pdb',)
        target_path = os.path.join(os.getcwd(), 'temp')
    
        if not os.path.exists(os.path.join(target_path, os.path.basename(original_path))):
            shutil.move(original_path, target_path)
        # If Download Success
        return True
    else:
        # if Download failed
        return False

# for a given pdb code it returns: number of alpha helics and number of beta sheets
def count_secondary_structure(pdb_id ):

    pdb_file = f"temp/{pdb_id}.pdb"
    # Load the fetched protein structure using MDTraj
    protein = md.load(pdb_file)
    # Compute secondary structure using DSSP
    dssp = md.compute_dssp(protein)
    # Count alpha helices and beta sheets
    num_alpha_helices = sum(True for ss in dssp[0] if ss in ['H', 'G', 'I'])
    num_beta_sheets = sum(True for ss in dssp[0] if ss in ['E', 'B'])
    
    return f"Number of alpha helices: {str(num_alpha_helices)}, Number of beta sheets: {str(num_beta_sheets)}"

## RCSB API : it returns the name of the protein.
def _RCSB_API(pdb = pdb_code):
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

## AA composition
def analyze_aa_composition(protein_sequence):
    """
    Analyze the amino acid composition of a protein sequence.

    Args:
    - protein_sequence (str): The protein sequence in single-letter amino acid code.

    Returns:
    - dict: A dictionary containing the count of each amino acid.
    """
    # Create a ProteinAnalysis object from the protein sequence
    protein_analysis = ProteinAnalysis(protein_sequence)

    # Get the amino acid composition
    aa_composition = protein_analysis.count_amino_acids()

    return aa_composition


## It returns the radius of gyration for a given pdb file
def Radius_of_Gy(pdb_id):
            # Step 1: Load the protein structure
        u = mda.Universe(f"temp/{pdb_id}.pdb")

        # Step 2: Select protein atoms
        protein_atoms = u.select_atoms("protein")

        # Step 3: Compute the center of mass (COM) of the protein
        com = protein_atoms.center_of_mass()

        # Step 4: Calculate the distance of each atom from the COM
        distances = protein_atoms.positions - com

        # Step 5: Compute the mean square distance from the COM
        msd = (distances ** 2).sum(axis=1).mean()

        # Step 6: Calculate the radius of gyration (Rg)
        rg = (msd ** 0.5)
        formatted_number = format(rg, '.2f')
        return f"Radius of gyration (Rg): {str(formatted_number)}"

## it calculates the contact map of a given pdb file
def generate_contact_map(pdb_id):
    # Load the protein structure using BioPython
    parser = PDBParser()
    structure = parser.get_structure('protein', f"temp/{pdb_id}.pdb")
    
    # Extract coordinates of alpha carbons
    alpha_carbons = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue):
                    try:
                        alpha_carbons.append(residue['CA'].get_coord())
                    except KeyError:
                        continue
    
    # Compute distances between alpha carbons
    contact_map = np.zeros((len(alpha_carbons), len(alpha_carbons)), dtype=np.float32)
    for i, coord_i in enumerate(alpha_carbons):
        for j, coord_j in enumerate(alpha_carbons):
            contact_map[i, j] = np.linalg.norm(coord_i - coord_j)
    
    # Plot the contact map
    base_name = os.path.basename(f"temp/{pdb_id}.pdb")
    plt.imshow(contact_map, cmap='Greys', origin='upper')
    plt.colorbar(label='Contacts')
    plt.xlabel('Residue Index',fontsize=20)
    plt.ylabel('Residue Index',fontsize=20)
    plt.title('Contact Map',fontsize=20)
    plt.savefig(f'contactmap/{base_name}_cmap.png', bbox_inches='tight', dpi=300)
    plt.show()
    ## returing the path of contact map
    path_contact_map = os.path.join(os.getcwd(),'contactmap')
    return path_contact_map


## It creates the input for a language model. So that LLM can create the summary. 
def seq_input_generator_llm(pdb_id):
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
    }
    # Load the protein structure file 
    u = mda.Universe(f"./temp/{pdb_id}.pdb")

    # Select the protein atoms
    protein_atoms = u.select_atoms("protein")

    # Get the sequence of the protein
    chain_ids = [i.segid for i in protein_atoms.segments] 

    chain_sequences = {}

    for chain_id in chain_ids:
        chain_residues = protein_atoms.select_atoms(f"segid {chain_id}").residues.resnames
        sequence_one_letter = ''.join([three_to_one[res] for res in chain_residues])
        chain_sequences[chain_id] = sequence_one_letter

    ## protein seuqence 
    protein_sequence = "".join([chain_sequences[key] for key in chain_ids])
    aa_composition = analyze_aa_composition(protein_sequence)

    ## get Radius of gyration
    radius_of_Gy = Radius_of_Gy(pdb_code)

    ## get Secondary structure details
    secondary_structure = count_secondary_structure(pdb_code)

    ## it gets the protein name
    protein_name = _RCSB_API(pdb_code)

    qns1 = "can you give the summary of the protein using the following features in a passge? Also, explain how one can interpret the amino acid composition, radius of gyratio, and secondary structure."
    feature_1 = f'name of the protein: {protein_name}'
    feature_2 = f'{secondary_structure}'
    feature_3 = f'Amino acid composition: {aa_composition}'
    feature_4 = f'radius of gyration: {radius_of_Gy}'

    input_to_LLM = f'{qns1}. {feature_1}. {feature_2}. {feature_3}. {feature_4}'

    return input_to_LLM # it will later be send to the LLM for writing the summary
     
#Root directry
@app.route('/')
def index():  
    return render_template('index.html', final_result = "Welcome To the Program. With this application, you can do the following things such as loading P D B file from online, compute the radius of gyration of the protein, calclate the number of alph helix and bete sheets of the protein, print the contact map, and finnaly get the summary of the protien using large language model. To begin with, you can load a  P D B structure using the command  'load' followed by the file name")


#Flask application
@app.route('/output', methods = ['GET','POST'])
def output():
    #loading PDB file in Pymol
    print("#####################")
    input_data = request.form.get('input_data')
    if 'load' in input_data:
        global pdb_code
        pdb_code = input_data.split()[1]
        pdb_response = download_pdb(pdb_code) ## downloading pdb file into temp folder
        respose = requests.get(f'http://localhost:8080/apply/load?pdb={pdb_code}')
        if pdb_response and respose.status_code == 200:
            u = mda.Universe(os.path.join(os.getcwd(),'temp',f'{pdb_code}.pdb'))
            u = u.select_atoms('protein')
            ## Chain IDs
            chain_ids = [i.segid for i in u.segments] 

            N_term = []
            C_term = []

            ## Extracting N and C terminal resids
            for chain in chain_ids:
                ch = u.select_atoms(f"segid {chain}")
                n_term = ch.residues.resids[0]
                c_term = ch.residues.resids[-1]
                N_term.append(n_term)
                C_term.append(c_term)
                # print(n_term, c_term)
            
            # int to string conversion
            n_term_lst = [str(item) for item in N_term]
            c_term_lst = [str(item) for item in C_term]

            # preparing arguments
            n_term = "_".join(n_term_lst)
            c_term = "_".join(c_term_lst)
            chain_ids_lst = "_".join(chain_ids)
            #http request
            respose = requests.get(f'http://localhost:8080/apply/sphere?n_term_lst={n_term}&c_term_lst={c_term}&chain_lst={chain_ids_lst}')
            return render_template('index.html', final_result = f"PDB file {_RCSB_API(pdb_code).replace('*', '')} with a code name {' '.join(list(pdb_code))} is loaded into pymol successfully. In this visualization, the N-terminal and C-terminal of the protein are denoted by spheres. If you'd like to compute the radius of gyration of the protein, type the command rg")
        else:
            return render_template('index.html', final_result = f"PDB file {pdb_code} is  not loaded into Pymol. Please, try again with a valid code by typing the command load followd by a P D B code")
            
    #Deliting PDB
    elif input_data == 'delete':
        respose = requests.post('http://localhost:8080/apply/delete')
        # global pdb_code
        temp = pdb_code
        pdb_code = None
        return render_template('index.html', final_result = f"Your loaded file {' '.join(list(temp))} is deleted")
    #Zoomin
    elif input_data == 'zoomin' or input_data == 'zoom in':
        if pdb_code:
            respose = requests.post('http://localhost:8080/apply/zoomin')
            return render_template('index.html', final_result = "Zooming in")
        else:
            return render_template('index.html', final_result = "Please load any file")
    #Zoomout
    elif input_data == 'zoomout' or input_data == 'zoom out':
        if pdb_code:
            respose = requests.post('http://localhost:8080/apply/zoomout')
            return render_template('index.html', final_result = "Zooming out")
        else:
            return render_template('index.html', final_result = "Please load any file")
    #Rotate x
    elif input_data == 'rotate_x' or input_data == 'rotatex':
        if pdb_code:
            respose = requests.get('http://localhost:8080/apply/rotate_x?deg=5')
            return render_template('index.html', final_result = "Rotating in the X axis")
        else:
            return render_template('index.html', final_result = "Please load any file")
    #Rotate y
    elif input_data == 'rotate_y' or input_data == 'rotatey':
        if pdb_code:
            respose = requests.get('http://localhost:8080/apply/rotate_y?deg=5')
            return render_template('index.html', final_result = "Rotating in the Y axis")
        else:
            return render_template('index.html', final_result = "Please load any file")
    #Rotate z
    elif input_data == 'rotate_z' or input_data == 'rotatez':
        if pdb_code:
            respose = requests.get('http://localhost:8080/apply/rotate_z?deg=5', final_result = "Rotating in the Z axis")
            return render_template('index.html')
        else:
            return render_template('index.html', final_result = "Please load any file")

    #Reset
    elif input_data.lower() == 'reset':
        respose = requests.post('http://localhost:8080/apply/reset')
        return render_template('index.html', final_result = "All Settings are set to the default")
    #Orient
    elif input_data.lower() == 'orient':
        respose = requests.post('http://localhost:8080/apply/orient')
        return render_template('index.html')

    ## N and C term Marker
    elif input_data.lower() == 'marker':
        u = mda.Universe(os.path.join(os.getcwd(),'temp',f'{pdb_code}.pdb'))
        u = u.select_atoms('protein')
        ## Chain IDs
        chain_ids = [i.segid for i in u.segments] 

        N_term = []
        C_term = []

        ## Extracting N and C terminal resids
        for chain in chain_ids:
            ch = u.select_atoms(f"segid {chain}")
            n_term = ch.residues.resids[0]
            c_term = ch.residues.resids[-1]
            N_term.append(n_term)
            C_term.append(c_term)
            # print(n_term, c_term)
        
        # int to string conversion
        n_term_lst = [str(item) for item in N_term]
        c_term_lst = [str(item) for item in C_term]

        # preparing arguments
        n_term = "_".join(n_term_lst)
        c_term = "_".join(c_term_lst)
        chain_ids_lst = "_".join(chain_ids)
        #http request
        respose = requests.get(f'http://localhost:8080/apply/sphere?n_term_lst={n_term}&c_term_lst={c_term}&chain_lst={chain_ids_lst}')
        return render_template('index.html')
    
    ## Returns Secondary structure to web application
    elif  'secondary' in input_data.lower():
        output =  f'{count_secondary_structure(pdb_code)}' # now it works
        return render_template('index.html', final_result = f'{output}. If you would like to print the contact map of the protein type of command contactmap')
    
    # Returns radius of gyration to the web application
    elif 'radius of gyration' in input_data.lower() or 'rg' in input_data.lower():
        output = f'{Radius_of_Gy(pdb_code)}'
        output = output.replace("Radius of gyration (Rg): ", "")
        prot_name = _RCSB_API(pdb_code)
        prot_name = prot_name.replace("-"," ").replace("*"," ")
        return render_template('index.html', final_result = f'Radius of gyration of protein {prot_name} is {output}. If you want to get the details of secondary structure of the protein like number of alpha helix, number of beta sheets, type the command secondary')
    
    ## it prints contact map
    elif input_data.lower() == 'contact map' or input_data.lower() == 'contactmap':
        # output = f'{Contact_Map(pdb_code)}'
        prot_name = _RCSB_API(pdb_code)
        prot_name = prot_name.replace("-"," ").replace("*"," ")
        path_cmap = generate_contact_map(pdb_code)
        return render_template('index.html', final_result = f'Contact map for the  {prot_name} is created and saved in {path_cmap}. Finally, If you want to get the overall summary of the protein with large language model, type the command LLM')

    
    ## it is related to Large language model
    elif input_data.lower().replace(' ','') == "llm":
        llm_question = seq_input_generator_llm(pdb_code)
        ### Language Model
        # return llm_question

        ## if the following is empty, type your API key
        GOOGLE_API_KEY = "AIzaSyCOxZ9dCoofIUg7tSVYGVA1JDSSfX7GD44"
        genai.configure(api_key= GOOGLE_API_KEY)

        ## function to load gemini pro model and get resposes
        model = genai.GenerativeModel("gemini-pro")

        ## LLM API## 
        def get_gemini_response(question):
            respose = model.generate_content(question)
            print(respose.text)
            return respose.text
        LLM_Response = get_gemini_response(llm_question)
        LLM_Response = LLM_Response.replace('*', '')
        LLM_Response += ". If you want to hear this again press zero,Thank you!"
        return render_template('index.html',final_result=LLM_Response)
    # return respose.text



## Flask Main Function
if __name__== '__main__':
    app.debug = True
    app.run(host = '0.0.0.0', port=5000)
