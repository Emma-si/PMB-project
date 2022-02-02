import gzip
import tarfile

def zip_vcf_file(vcf_file):

    vcf_zipped_file = vcf_file + ".gz"

    with open(vcf_file, 'rb') as vcf_input:
        with gzip.open( vcf_file + ".gz", 'wb') as zipped_vcf_file:
            zipped_vcf_file.writelines(vcf_input)
    
    return vcf_zipped_file

def create_snp_file(vcf_zipped_file):

    snp_filename = vcf_zipped_file[:-7]
    snp_name = "SNP_" + snp_filename

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
        manifest_file.write("filename = " + vcf_zipped_file.split("/")[-1] + "\n")

    vcf_tar_file = snp_filename + ".tar.gz"
    
    with tarfile.open( vcf_tar_file, "w:gz" ) as tar:
        tar.add("manifest.ini")
        tar.add(vcf_zipped_file, arcname = vcf_zipped_file.split("/")[-1])

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

def merge_tmp_tables(table_name, seg_tables_list):
    seg_tables_list.sort()

    with open(table_name,'wb+') as table:
        for f in seg_tables_list:
            with open(f,'rb') as seg:
                shutil.copyfileobj(seg, table)