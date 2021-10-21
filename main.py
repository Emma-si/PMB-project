import sys
import argparse
import pyGeno.bootstrap as B
from pyGeno.importation.SNPs import *
from pyGeno.Genome import *
import multiprocessing as mp   
import time
import tqdm 




pbar = None   #variabile globale della barra di progresso

def protein_worker(protein_id):   #viene chiamata dal processo
    # Tiro fuori la sequenza dalle proteine 

   
    #importo il genoma nel processo (non posso accederci da altro processo)
    myGeno = Genome(name = 'GRCh37.75', SNPs="mySNPmerged02")
    #prendo la proteina con l'ID passato in argomento
    protein = myGeno.get(Protein, id = protein_id)[0]
    #definisco il dizionario
    prot_dict = {}
    prot_dict["id"] = protein.id
    prot_dict["transcript_id"] = protein.transcript.id
    prot_dict["name"] = protein.name
    prot_dict["chromosome_number"] = protein.chromosome.number
    try:
        prot_dict["sequence"] = protein.sequence 
    except:
        prot_dict["sequence"] = "none"  #stampa none se non è presente sequenza

    return prot_dict

def write_on_file(prot):    #gli passo il dizionario (output del worker)

    with open('tabella_merged02.txt', 'a') as file:
        file.write(prot["id"] + "\t" +
                prot["transcript_id"] + "\t" +
                prot["name"] + "\t" +
                prot["chromosome_number"] + "\t" +
                prot["sequence"] + "\t" +
                '\n')
        pbar.update(1)  #update della barra di progresso






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

     proteins = myGeno.get(Protein)
    number_of_proteins = len(proteins) 


    # Progress bar definition
    pbar = tqdm.tqdm(total=number_of_proteins)
    

    #creo una lista di id delle proteine (che posso passare)
    protein_ids = [p.id for p in proteins]

    #creo n processi
    pool = mp.Pool(mp.cpu_count()-3)

   
    for protein_id in protein_ids:
        pool.apply_async(protein_worker, args = (protein_id, ), callback = write_on_file) 
        #applica in modo asincrono le funzioni protein worker (con i suoi args) e callback

    pool.close() #chiude il processo corrente
    pool.join()  #aspetta che tutti siano chiusi
    