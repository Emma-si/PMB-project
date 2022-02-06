import gzip
import tarfile

def zip_vcf_file(vcf_file : str):
    """
    Function that zip the vcf file in .gz format
    Input:
        vcf_file : path of the vcf file
    Output:
        vcf_zipped_file : path of the gzipped file
    """
    vcf_zipped_file = vcf_file + ".gz"

    with open(vcf_file, 'rb') as vcf_input:
        with gzip.open( vcf_zipped_file, 'wb') as vcf_zipped_output:
            vcf_zipped_output.writelines(vcf_input)
    
    return vcf_zipped_file

def create_snp_file(vcf_zipped_file : str):
    """
    Function that create the snp tar file according to the Pygeno parser
    Input:
        vcf_zipped_file : path of the vcf zipped file
    Outputs:
        vcf_tar_file : path of the new compressed file in the format .tar.gz
        snp_name : name of the created snp
    """
    # Definition of the snp name to associate to each individual
    snp_filename = vcf_zipped_file[:-7]
    individual_name = snp_filename.split("/")[-1]
    snp_name = "SNP_" + individual_name

    vcf_filename = vcf_zipped_file.split("/")[-1]

    # Definition and creation of the manifest.ini file required by Pygeno
    with open("manifest.ini", 'w') as manifest_file:
        manifest_file.write("[package_infos]\n")
        manifest_file.write("description = SNPs for reasearch purposes\n")
        manifest_file.write("maintainer = Person\n")
        manifest_file.write("maintainer_contact = Person contact\n")
        manifest_file.write("version = 1\n\n")

        manifest_file.write("[set_infos]\n")
        manifest_file.write("species = human\n")
        manifest_file.write("name = " + snp_name + "\n")
        manifest_file.write("type = dbSNPSNP\n")
        manifest_file.write("source = Source\n\n")

        manifest_file.write("[snps]\n")
        manifest_file.write("filename = " + vcf_filename + "\n")

    # Creation of the archive according to Pygeno requirements
    vcf_tar_file = snp_filename + ".tar.gz"
    
    with tarfile.open( vcf_tar_file, "w:gz" ) as tar:
        tar.add("manifest.ini")
        tar.add(vcf_zipped_file, arcname = vcf_filename)

    return vcf_tar_file, snp_name

def split_list(complete_list : list, n : int):
    """
    Function that split a list in n chunks
    Inputs:
        complete_list : list to be divided
        n : number of chuncks to split the list in
    Outputs:
        chunks : list of sublists
    """
    k, m = divmod(len(complete_list), n)
    chunks = [complete_list[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n)]

    return chunks

def merge_tmp_tables(table_name : str, tmp_tables_list : list):
     """
    Function that merges together the temporary files created by the 
    different processes into a single final table
    Inputs:
        table_name : filename of the final table
        tmp_tables_list : list of temporary tables filenames
    """
    # Sort list to obtain ordered chromosomes in the final table
    tmp_tables_list.sort()

    with open(table_name,'wb+') as table:
        for f in tmp_tables_list:
            with open(f,'rb') as tmp:
                shutil.copyfileobj(tmp, table)