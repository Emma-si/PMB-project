import argparse
import glob
import utils
import multiprocessing as mp   
import tqdm 
import time
import os
import warnings

from pyGeno.SNPFiltering import SNPFilter
from pyGeno.SNPFiltering import SequenceSNP
import pyGeno.bootstrap as B
from pyGeno.importation.SNPs import importSNPs, deleteSNPs
from pyGeno.Genome import Genome, Protein



def protein_worker( args ):
    '''
    Input: args = [process id, reference genome name, snp name, list protein ids]

    This function extracts the modified sequence of a list of proteins
    and saves it into a temporary tab separated file.
    '''

    process_id = args[0]
    reference_genome = args[1]
    snp_name = args[2]
    protein_ids = args[3]

    progress_bar = tqdm.tqdm(total=len(protein_ids), desc = "Process %d" % process_id, position = process_id)

    table_name = "tmp%d_%s.txt" % (process_id, snp_name)

    with open(table_name, 'a+') as file :
        for row in genome_to_proteinlist_generator(protein_ids, reference_genome, snp_name) :
            file.write(row) 
            progress_bar.update(1) 




def genome_to_proteinlist_generator(protein_ids, reference_genome, snp_name) :
    '''
    Input: process id, reference genome name, snp name

    Generator that creates a row with the protein information extracted 
    from the genome
    '''
    
    myGeno = Genome(name = reference_genome, SNPs = snp_name, SNPFilter = MyFilter())

    for protein_id in protein_ids:
        protein = myGeno.get(Protein, id = protein_id)[0]


        transcript_id = protein.transcript.id
        name = protein.name
        chromosome_number = protein.chromosome.number
        # Assign none if non-existent protein sequence.
        try:
            sequence = protein.sequence 
        except:
            sequence = "none"

        row = protein_id + "\t" + \
              transcript_id + "\t" + \
              name + "\t" + \
              chromosome_number + "\t" + \
              sequence + "\t" + \
              '\n'
        yield row



class MyFilter(SNPFilter) :
    '''
    Custom filter that, given the snp alteration, subsitutes it into the protein sequence. 
    '''
    def init(self) :
        pass

    def filter(self, chromosome, **SNPs) :
    
        for snpSet, snp in SNPs.items() :
            return SequenceSNP(snp.alt)
        
        return None



if __name__ == "__main__":

    # Command line parser
    argparser = argparse.ArgumentParser(prog = 'PMB_project', description = 'PMB_project')
    argparser.add_argument('--genome', type = str, help = 'set the reference genome', default = "", required = True)
    argparser.add_argument('--vcf_file', type = str, help = 'set the vcf file name', default = "", required = True)
    argparser.add_argument('--num_processes', type = int, help = 'set the number of processes', default = 0, required = False)


    # Parameters
    args = argparser.parse_args()
    genome_name = args.genome
    vcf_file = args.vcf_file
    n_processes = args.num_processes

    if n_processes == 0: 
        warnings.warn("No number of processes specified, the program will use all available cpu. This may slow down your pc.")
        n_processes = mp.cpu_count()


    start_time = time.time()

    # Create snp file
    gzipped_vcf_file = utils.zip_vcf_file(vcf_file)
    snp_file, snp_name = utils.create_snp_file(gzipped_vcf_file)

    # Print message if already existing SNP
    try:
        importSNPs(snp_file)
    except KeyError:
        print("Existing SNP")

    # Import the reference genome and if not present download it
    try:
        myGeno = Genome(name = genome_name, SNPs = snp_name)
    except KeyError:
        B.importRemoteGenome(genome_name)
        myGeno = Genome(name = genome_name, SNPs = snp_name)

    
    proteins = myGeno.get(Protein)
    protein_ids = [p.id for p in proteins]

    # Create the sublists to give in input to the different processes
    chunks = utils.split_list(protein_ids, n_processes)

    # Start processes in asyncronous way
    with mp.Pool(processes=n_processes) as pool:
        for i in range(n_processes):
            pool.apply_async(protein_worker, args = ([i, genome_name, snp_name, chunks[i]],))
        # Close the current process and wait for the other to close
        pool.close() 
        pool.join() 
    
    # Delete SNP to avoid memory issues
    deleteSNPs(snp_name)

    table_name = snp_name + ".txt"
    tmp_tables_list = glob.glob("tmp*.txt")
    utils.merge_tmp_tables(table_name, tmp_tables_list)

    # Delete temporary files to avoid memory issues
    for tmp_table in tmp_tables_list:
        os.remove(tmp_table)
    
    os.remove("manifest.ini")
    os.remove(gzipped_vcf_file)
    os.remove(snp_file)

    # Print of the total time used for the processing 
    end_time = time.time()
    total_time = (end_time - start_time)/60
    print("Genome processed in %.2f minutes" % total_time)



