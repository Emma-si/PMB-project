import argparse
import glob
import utils
import multiprocessing as mp   
import tqdm 

from pyGeno.SNPFiltering import SNPFilter
from pyGeno.SNPFiltering import SequenceSNP
import pyGeno.bootstrap as B
from pyGeno.importation.SNPs import *
from pyGeno.Genome import *


def protein_worker( args ):
    '''
    Input: reference genome name, snp name, list protein ids

    This function extracts the modified sequence of a list of proteins
    and saves it into a temporary tab separated file.
    '''
    # print(args)
    i = args[0]
    reference_genome = args[1]
    snp_name = args[2]
    protein_ids = args[3]
    pbar = tqdm.tqdm(total=len(protein_ids), desc = "Process %d"%i, position = i)

    table_name = "tmp%d_%s.txt" % (i, snp_name)

    # Import genome in the single process, required because it is not possible to access fro other processes
    myGeno = Genome(name = reference_genome, SNPs = snp_name, SNPFilter = MyFilter())

    for protein_id in protein_ids:
        # Selection of the single protein
        protein = myGeno.get(Protein, id = protein_id)[0]

        # Dictionary definition
        transcript_id = protein.transcript.id
        name = protein.name
        chromosome_number = protein.chromosome.number
        try:
            sequence = protein.sequence 
        except:
            sequence = "none"

        with open(table_name, 'a+') as file:
            file.write(protein_id + "\t" +
                    transcript_id + "\t" +
                    name + "\t" +
                    chromosome_number + "\t" +
                    sequence + "\t" +
                    '\n')

            # Update the progress bar
            pbar.update(1)  






class MyFilter(SNPFilter) :
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

    # Parameters
    args = argparser.parse_args()
    genome_name = args.genome
    vcf_file = args.vcf_file

    # Create snp file
    gzipped_vcf_file = utils.zip_vcf_file(vcf_file)
    snp_file, snp_name = utils.create_snp_file(gzipped_vcf_file)

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

    # Get proteins from genome
    proteins = myGeno.get(Protein)
    

    # Create a list of protein ids to pass to processes
    protein_ids = [p.id for p in proteins]

    n_processes = mp.cpu_count()-3
    chunks = split_list(protein_ids, n_processes)

    # Define the processes based on the available cpu
    pool = mp.Pool(n_processes)

    # Start processes in asyncronous way
    for i in range(n_processes):
        pool.apply_async(protein_worker, args = ([i, genome_name, snp_name, chunks[i]],))

    # Close the current process and wait for the other to close
    pool.close() 
    pool.join() 
    
    deleteSNPs(snp_name)

    table_name = snp_name + ".txt"
    seg_tables_list = glob.glob("tmp*.txt")

    merge_tmp_tables(table_name, seg_tables_list)


