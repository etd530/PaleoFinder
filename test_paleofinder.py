#### LIBS ####
import paleofinder as pf

#### Test for alien_index ####
# Test that it computes correct value when given two lists of e-values
def test_alien_index_standard_case():
	belonging_list = [0.0001, 0.1, 0.1, 1]
	nonbelonging_list = [0.1, 1, 0.1, 1]
	assert pf.alien_index(belonging_list, nonbelonging_list) == 4

# Test it returns zero if both lists consist entirely of zeroes
def test_alien_index_both_smallest_zero():
	belonging_list = [0, 0, 0, 0]
	nonbelonging_list = [0, 0, 0, 0]
	assert pf.alien_index(belonging_list, nonbelonging_list) == 0

# Test it gives -np.inf when there are no belonging hits
def test_alien_index_zero_belonging_hits():
	belonging_list = []
	nonbelonging_list = [0, 2, 0.5, 14]
	assert pf.alien_index(belonging_list, nonbelonging_list) == -pf.np.inf

# Test it gives np.inf when there are no non-belonging hits
def test_alien_index_zero_nonbelonging_hits():
	belonging_list = [0, 2, 0.5, 14]
	nonbelonging_list = []
	assert pf.alien_index(belonging_list, nonbelonging_list) == pf.np.inf

# Test it computes ok whith some zero values in both lists
def test_alien_index_both_smallest_zero():
	belonging_list = [0, 0.0001, 0, 0.1]
	nonbelonging_list = [0, 0.1, 0.1, 1]
	assert pf.alien_index(belonging_list, nonbelonging_list) == 310.6575773191778

# Test it gives -np.inf when the smallest non-belonging hit is zero
def test_alien_index_best_nonbelonging_is_zero():
	belonging_list = [0.0001, 0.2, 23, 2]
	nonbelonging_list = [0, 2, 0.5, 14]
	assert pf.alien_index(belonging_list, nonbelonging_list) == -303.4752371108451

# Test it gives np.inf when the smallest belonging hit is zero
def test_alien_index_best_nonbelonging_is_zero():
	belonging_list = [0, 2, 0.5, 14]
	nonbelonging_list = [0.0001, 0.2, 23, 2]
	assert pf.alien_index(belonging_list, nonbelonging_list) == 303.4752371108451
