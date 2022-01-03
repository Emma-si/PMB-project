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

if __name__ == "__main__":
    test_chuncks_same_length()