from hypothesis import given, settings, strategies as st
from utils import split_list
import os


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

    with open("test/test_split_list.txt", "r+") as original_file: 
        all_lines = original_file.readlines()
    
    n = 5
    chunks = split_list(all_lines, n)

    for i, chunk in enumerate(chunks):
        chunk_filename = "test/test_split_list_chunk%d.txt" % i
        with open(chunk_filename, "a+") as chunk_file:
            for row in chunk: 
                chunk_file.write(row)

    original_file = open("test/test_split_list.txt", "r+")

    for i in range(n):
       chunk_filename = "test/test_split_list_chunk%d.txt" % i 
       with open(chunk_filename, "r+") as chunk_file:
           for line in chunk_file:
               assert line == original_file.readline()
           
    original_file.close()
    for i in range(n):
        chunk_filename = "test/test_split_list_chunk%d.txt" % i 
        os.remove(chunk_filename)
        

def test_merge_tmp_tables():
    """
    Tests:
        - if the sum of tmp tables has the same length of the entire table
    """

    table_name = "test/SNP_chr21_NA20502.txt"
    seg_tables_list = ['test/tmp0_SNP_chr21_NA20502.txt', 
                        'test/tmp1_SNP_chr21_NA20502.txt', 
                        'test/tmp2_SNP_chr21_NA20502.txt', 
                        'test/tmp3_SNP_chr21_NA20502.txt', 
                        'test/tmp4_SNP_chr21_NA20502.txt', 
                        'test/tmp5_SNP_chr21_NA20502.txt', 
                        'test/tmp6_SNP_chr21_NA20502.txt', 
                        'test/tmp7_SNP_chr21_NA20502.txt', 
                        'test/tmp8_SNP_chr21_NA20502.txt']

    with open(table_name, 'r') as ft:
        table_lines = len(ft.readlines())


    seg_lines = 0
    for seg_table in seg_tables_list:
        with open(seg_table, 'r') as st:
            l = len(st.readlines())
            seg_lines += l

    assert table_lines == seg_lines
