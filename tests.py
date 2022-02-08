from hypothesis import given, settings, strategies as st
from utils import split_list


@given(st.lists(st.integers(), min_size=500, max_size=2000), st.integers(min_value = 10, max_value=20))
@settings(max_examples = 50)
def test_chuncks_same_length(l, n):
    """
    Tests:
        - if all sublists (chunks) are of equal length
        - if the union of the sublists returns the original list
    """

    dim = len(l) // n
    good_dims = [dim - 1, dim, dim + 1]

    chunks = split_list(l, n)
    for chunk in chunks:
        assert len(chunk) in good_dims
    
    union = chunks[0]
    for chunk in chunks[1:]:
        union += chunk
    assert union == l


def test_merge_tmp_tables(table_name, seg_tables_list):
    """
    Tests:
        - if the sum of tmp tables has the same length of the entire table
    """
    with open(table_name, 'r') as ft:
        table_lines = len(ft.readlines())


    seg_lines = 0
    for seg_table in seg_tables_list:
        with open(seg_table, 'r') as st:
            l = len(st.readlines())
            seg_lines += l

    assert table_lines == seg_lines

if __name__ == "__main__":
    test_chuncks_same_length()

    test_merge_tmp_tables("test/SNP_chr21_NA20502.txt", ['test/tmp0_SNP_chr21_NA20502.txt', 
                                                          'test/tmp1_SNP_chr21_NA20502.txt', 
                                                          'test/tmp2_SNP_chr21_NA20502.txt', 
                                                          'test/tmp3_SNP_chr21_NA20502.txt', 
                                                          'test/tmp4_SNP_chr21_NA20502.txt', 
                                                          'test/tmp5_SNP_chr21_NA20502.txt', 
                                                          'test/tmp6_SNP_chr21_NA20502.txt', 
                                                          'test/tmp7_SNP_chr21_NA20502.txt', 
                                                          'test/tmp8_SNP_chr21_NA20502.txt'])
