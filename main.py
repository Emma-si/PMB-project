import sys
import argparse
import pyGeno.bootstrap as B
from pyGeno.importation.SNPs import *
from pyGeno.Genome import *
import multiprocessing as mp   
import time
import tqdm 



# TO DO: ELIMINATE GLOBAL VARIABLE
pbar = None   

def protein_worker(protein_id, genome_name, snp_name):   
    '''
    Input: genome name, snp name, protein id

    Output: protein dictionary

    This function extracts the modified sequence of the single protein.
    '''
 
    # Import genome in the single process, required because it is not possible to access fro other processes
    myGeno = Genome(name = genome_name, SNPs = snp_name)

    # Selection of the single protein
    protein = myGeno.get(Protein, id = protein_id)[0]

    # Dictionary definition
    prot_dict = {}
    prot_dict["id"] = protein.id
    prot_dict["transcript_id"] = protein.transcript.id
    prot_dict["name"] = protein.name
    prot_dict["chromosome_number"] = protein.chromosome.number
    try:
        prot_dict["sequence"] = protein.sequence 
    except:
        prot_dict["sequence"] = "none"  

    return prot_dict, snp_name


def write_on_file(prot, snp_name): 
    '''
    Input: snp name, protein dictionary

    Callback function that write on file the sequences.
    '''

    table_name = snp_name + ".txt"
    
    with open(table_name, 'a') as file:
        file.write(prot["id"] + "\t" +
                prot["transcript_id"] + "\t" +
                prot["name"] + "\t" +
                prot["chromosome_number"] + "\t" +
                prot["sequence"] + "\t" +
                '\n')

        # Update the progress bar
        pbar.update(1) 






if __name__ == "__main__":

    # Command line parser
    argparser = argparse.ArgumentParser(prog = 'PMB_project', description = 'PMB_project')
    argparser.add_argument('--genome', type = str, help = 'set the reference genome', default = "", required = True)
    argparser.add_argument('--snp', type = str, help = 'set the snp path', default = "", required = False)
    argparser.add_argument('--snp_name', type = str, help = 'set the snp name', default = "", required = True)

    # Parameters
    args = argparser.parse_args()
    genome_name = args.genome
    snp_path = args.snp
    snp_name = args.snp_name

    # If the user specify the snp path to insert, the snp is imported
    # TO DO ADD TRY EXCEPT
    if snp_path != '':
        importSNPs(snp_path)

    # Import the reference genome and if not present download it
    try:
        myGeno = Genome(name = genome_name, SNPs = snp_name)
    except KeyError:
        B.importRemoteGenome(genome_name)
        myGeno = Genome(name = genome_name, SNPs = snp_name)

    # Get proteins from genome
    proteins = myGeno.get(Protein)
    number_of_proteins = len(proteins) 

    # Progress bar definition
    pbar = tqdm.tqdm(total=number_of_proteins)

    # Create a list of protein ids to pass to processes
    protein_ids = [p.id for p in proteins]

    # Define the processes based on the available cpu
    pool = mp.Pool(mp.cpu_count()-3)

    # Start processes in asyncronous way
    for protein_id in protein_ids:
        pool.apply_async(protein_worker, args = (protein_id, genome_name, snp_name, ), 
                                         callback = write_on_file) 

    # Close the current process and wait for the other to close
    pool.close() 
    pool.join() 
    