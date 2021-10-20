import sys
import argparse
import pyGeno.bootstrap as B
from pyGeno.importation.SNPs import *




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

    # Import the reference genome and if not present is downloaded
    try:
        myGeno = Genome(name = genome_name, SNPs = snp_name)
    except KeyError:
        B.importRemoteGenome(genome_name)

    