from utils import split_list, merge_tmp_tables
import os
from main import genome_to_proteinlist_generator
from pyGeno.Genome import Genome, Protein
import pyGeno.bootstrap as B



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



  