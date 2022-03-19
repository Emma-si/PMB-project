from utils import split_list, merge_tmp_tables, create_snp_file, create_manifest_file, zip_vcf_file
import os
from main import genome_to_proteinlist_generator
from pyGeno.Genome import Genome, Protein
import pyGeno.bootstrap as B
import tarfile
import gzip



def test_chuncks_division_1():
    """
    Tests:
        - if, given a known even list, all sublists (chunks) are of equal length
        - if the sublists content is the expected one 
    """

    l = [1,2,3,4,5,6]
    n = 2
    chunks = split_list(l, n)
    assert chunks[0] == [1,2,3]
    assert chunks[1] == [4,5,6]

def test_chuncks_division_2():
    """
    Tests:
        - if, given a known odd list, all sublists (chunks) are of acceptable length
        - if the sublists content is the expected one 
    """

    l = [1,2,3,4,5,6,7]
    n = 2
    chunks = split_list(l, n)
    assert chunks[0] == [1,2,3,4]
    assert chunks[1] == [5,6,7]

def test_chuncks_division_3():
    """
    Tests:
        - if, given a known list, all sublists (chunks) are of acceptable length
        - if the sublists content is the expected one 
    """

    l = [1,2,3,4,5,6,7,8,9,10]
    n = 3
    chunks = split_list(l, n)
    assert chunks[0] == [1,2,3,4]
    assert chunks[1] == [5,6,7]
    assert chunks[2] == [8,9,10]

def test_chunks_division_file():
    """
    Test:
        Apply the split_list function to a known file in order to check if the content 
        and lenght of the created files matches the original one line by line.
    """

    # Get list of lines in the original file
    with open("test/test_split_list.txt", "r+") as original_file: 
        all_lines = original_file.readlines()
    
    # Divide the list of all lines in n sublists
    n = 5
    chunks = split_list(all_lines, n)

    # Save the lines in the sublists in n temporary files 
    for i, chunk in enumerate(chunks):
        chunk_filename = "test/test_split_list_chunk%d.txt" % i
        with open(chunk_filename, "a+") as chunk_file:
            for row in chunk: 
                chunk_file.write(row)

    # Open again the original file in read mode
    original_file = open("test/test_split_list.txt", "r+")

    # Compare line by line the original file and the temporary files
    for i in range(n):
       chunk_filename = "test/test_split_list_chunk%d.txt" % i 
       with open(chunk_filename, "r+") as chunk_file:
           for line in chunk_file:
               assert line == original_file.readline()
    
    # Delete temporary files
    original_file.close()
    for i in range(n):
        chunk_filename = "test/test_split_list_chunk%d.txt" % i 
        os.remove(chunk_filename)


def test_merge_tmp_tables():
    """
    Test:
        Apply the merge_tmp_tables function to a known splitted file in order to check if the content 
        and lenght of the created temporary tables matches the original one after the merge.
    """

    # Get list of lines in the original file
    table_name = "test/SNP_chr21_NA20502.txt"
    with open(table_name, "r+") as original_file: 
        all_lines = original_file.readlines()
    
    # Divide the list of all lines in n sublists
    n = 8
    chunks = split_list(all_lines, n)

    # Save the lines in the sublists in n temporary files and keep track of tmp filenames
    tmp_files = []
    for i, chunk in enumerate(chunks):
        chunk_filename = "test/SNP_chr21_NA20502_tmp%d.txt" % i
        tmp_files.append(chunk_filename)
        with open(chunk_filename, "a+") as chunk_file:
            for row in chunk: 
                chunk_file.write(row)

    # Merge temporary files
    merged_file = "test/merged_SNP_chr21_NA20502.txt"
    merge_tmp_tables(merged_file, tmp_files)

    # Open again the original file in read mode
    original_file = open("test/SNP_chr21_NA20502.txt", "r+")

    # Compare line by line the original file and the merged file
    with open(merged_file, "r+") as merged:
        for line in merged:
            assert line == original_file.readline()

    # Delete temporary and merged files       
    original_file.close()
    for i in range(n):
        chunk_filename = "test/SNP_chr21_NA20502_tmp%d.txt" % i 
        os.remove(chunk_filename)
    os.remove(merged_file)


def test_genome_to_proteinlist_generator():
    """
    Tests:
        If the output appears as expected. 
        If, from the reference genome, the rows created are 
        a list of five elements separated from tab spacing.
    """

    # Load genome
    try:
        myGeno = Genome(name = "GRCh37.75")
    except KeyError:
        B.importRemoteGenome("GRCh37.75")
        myGeno = Genome(name = "GRCh37.75")

    proteins = myGeno.get(Protein)
    protein_ids = [p.id for p in proteins]

    # Assert that the rows are composed of five elements
    for i in range(len(protein_ids)):
        row = genome_to_proteinlist_generator(protein_ids, "GRCh37.75", [])
        row = row.__next__().strip()
        assert len(row.split("\t")) == 5

        # Too much time required to analyze all the genome 
        if i == 1000:
            break


def test_zip_vcf_file():
    """
    Tests: 
        - if, given a known vcf file in input, the gzipped vcf file, 
        once extracted, is equal to the original line by line.
    """
    
    vcf_filename = "test/test_create_vcf.vcf"

    # Zip the input vcf file
    gzipped_vcf_file = zip_vcf_file(vcf_filename)

    # Open again the original file in read mode
    original_file = open(vcf_filename, "r+")

    # Extract the gzip and compare line by line with the original file 
    with gzip.open(gzipped_vcf_file, "rt") as zip_file:
        for line in zip_file:
            assert line == original_file.readline()

    original_file.close()
    os.remove(gzipped_vcf_file)



def test_create_manifest_file():
    """
    Tests: 
        - if, given a handwritten manifest.ini file in input, the created one
        is equal to it line by line.
    """

    test_manifest = "test/test_manifest.ini"

    # Definition of expected variables 
    snp_name = "SNP_test_create_vcf"
    vcf_filename = "test_create_vcf.vcf.gz"

    create_manifest_file(snp_name, vcf_filename)

    # Open the handwritten file in read mode
    test_file = open(test_manifest, "rt")

    # Compare line by line the handwritten file and the created one
    with open("manifest.ini", "rt") as f:
        for line in f:
            assert line == test_file.readline()

    test_file.close()
    os.remove("manifest.ini")


def test_create_snp_file():
    """
    Tests: 
        - if the created snp archive contains the two
        expected files (manifest.ini and the input gzipped vcf).
        - if the gzipped_vcf_file contained in the archive
        is equal to the original vcf file.
    """
    
    vcf_filename = "test/test_create_vcf.vcf"

    # Create snp file
    gzipped_vcf_file = zip_vcf_file(vcf_filename)
    snp_file, snp_name = create_snp_file(gzipped_vcf_file)

    # Delete gzipped file created by the zip_vcf_file function (already tested)
    os.remove(gzipped_vcf_file)

    # Open the archive 
    with tarfile.open(snp_file, "r:gz") as tar:
        contained_files = tar.getnames()

        # Assert that the archive contains exactly the two expected files
        assert len(contained_files) == 2
        assert "manifest.ini" in contained_files
        assert gzipped_vcf_file.split("/")[-1] in contained_files

        # Extract the contained files from the archive in the test folder
        tar.extractall("test/")

        # Assert that the extracted gzipped file is the same as 
        # the original vcf file line by line
        original_file = open(vcf_filename, "rt")
        extracted_gzipped_file = "test/" + contained_files[1]

        with gzip.open(extracted_gzipped_file, "rt") as f:
            for line in f:
                assert line == original_file.readline()

        original_file.close()
        # Remove extracted files from the test folder
        os.remove(extracted_gzipped_file)
        os.remove("test/manifest.ini")

    # Remove files created by the create_snp_file function
    os.remove(snp_file)
    os.remove("manifest.ini")

