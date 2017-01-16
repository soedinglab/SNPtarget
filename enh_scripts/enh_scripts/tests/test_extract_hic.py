import unittest
import scipy.sparse as sp
import numpy.testing as nt
import numpy as np

from enh_scripts.extract_hic import extract_background
from enh_scripts.rao_hic_to_matrix import matrix_to_diagdict

class HiCExtractTestCase(unittest.TestCase):

    def test_matrix43_to_diagdict(self):
        m, n = 4, 3
        mat = np.arange(m * n).reshape(m, n)
        sp_mat = sp.lil_matrix(mat)

        exp_dict = {
            0: [0, 4, 8],
            1: [1, 5],
            2: [2],
        }
        diag_dict = matrix_to_diagdict(sp_mat, 2)
        self.assertEqual(exp_dict.keys(), diag_dict.keys())
        for d, sp_vec in diag_dict.items():
            nt.assert_allclose(exp_dict[d], sp_vec.toarray().ravel())

    def test_matrix34_to_diagdict(self):
        m, n = 3, 4
        mat = np.arange(m * n).reshape(m, n)
        sp_mat = sp.lil_matrix(mat)

        exp_dict = {
            0: [0, 5, 10],
            1: [1, 6, 11],
            2: [2, 7],
        }
        diag_dict = matrix_to_diagdict(sp_mat, 2)
        self.assertEqual(exp_dict.keys(), diag_dict.keys())
        for d, sp_vec in diag_dict.items():
            nt.assert_allclose(exp_dict[d], sp_vec.toarray().ravel())

    def test_diagdict_too_large_d(self):
        m, n = 3, 4
        mat = np.arange(m * n).reshape(m, n)
        sp_mat = sp.lil_matrix(mat)

        diag_dict = matrix_to_diagdict(sp_mat, 10)
        self.assertEqual(len(diag_dict), 4)

    def test_new_extract46(self):
        m, n = 4, 6
        mat = np.arange(m*n).reshape(m, n)
        sp_mat = sp.lil_matrix(mat)

        test_cases = [
            # i, j, bw, exp_output
            (1, 1, 1, [0, 7, 14]),
            (1, 1, 2, [0, 7, 14, 21]),
            (1, 3, 2, [2, 9, 16, 23]),
            (0, 5, 3, [5]),
            (3, 3, 3, [0, 7, 14, 21]),
            (0, 0, 3, [0, 7, 14, 21]),
            (0, 2, 2, [2, 9, 16]),
        ]

        for i, j, bw, exp_output in test_cases:
            diag_dict = matrix_to_diagdict(sp_mat, j-i)
            bg = extract_background(diag_dict, i, j, m, n, bw)
            nt.assert_allclose(exp_output, np.ravel(bg))

    def test_new_extract64(self):
        m, n = 6, 4
        mat = np.arange(m*n).reshape(m, n)
        sp_mat = sp.lil_matrix(mat)

        test_cases = [
            # i, j, bw, exp_output
            (1, 1, 1, [0, 5, 10]),
            (1, 1, 2, [0, 5, 10, 15]),
            (1, 3, 2, [2, 7]),
            (0, 3, 3, [3]),
            (3, 3, 3, [0, 5, 10, 15]),
            (0, 0, 3, [0, 5, 10, 15]),
            (0, 2, 2, [2, 7]),
        ]

        for i, j, bw, exp_output in test_cases:
            diag_dict = matrix_to_diagdict(sp_mat, j-i)
            bg = extract_background(diag_dict, i, j, m, n, bw)
            nt.assert_allclose(exp_output, np.ravel(bg))
