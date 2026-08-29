import os
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from app.tools import (
    calcule_Delta_G,
    calcule_probability,
    read_arguments,
    read_data,
    save_data,
)


class ToolsTestCase(unittest.TestCase):
    def test_read_data_concatenates_duplicate_states_without_merging(self):
        with tempfile.TemporaryDirectory() as directory:
            path_a = os.path.join(directory, 'a.tsv')
            path_b = os.path.join(directory, 'b.tsv')
            pd.DataFrame({'RG': [1.0, 1.0], 'RMSD': [2.0, 2.0]}).to_csv(path_a, sep='\t', index=False)
            pd.DataFrame({'RG': [1.0], 'RMSD': [2.0]}).to_csv(path_b, sep='\t', index=False)

            data = read_data(path_a, path_b)

            self.assertEqual(len(data), 3)
            self.assertEqual(data['model'].tolist(), ['a', 'a', 'b'])

    def test_probability_is_aligned_with_each_row(self):
        data = pd.DataFrame(
            {'RG': [1.0, 1.0, 2.0], 'RMSD': [2.0, 2.0, 3.0]}
        )

        self.assertEqual(calcule_probability(data), [2 / 3, 2 / 3, 1 / 3])

    def test_delta_g_rejects_invalid_probability(self):
        with self.assertRaises(ValueError):
            calcule_Delta_G([0.0], 300, ['a'])

    def test_save_data_writes_tsv(self):
        with tempfile.TemporaryDirectory() as directory:
            save_data([1.0], [2.0], [3.0], directory, ['a'], 'result')
            output = pd.read_csv(os.path.join(directory, 'result_out.dat'), sep='\t')

            self.assertEqual(output.columns.tolist(), ['RMSD', 'RG', 'DG', 'model'])
            self.assertEqual(output.iloc[0].to_dict(), {'RMSD': 1.0, 'RG': 2.0, 'DG': 3.0, 'model': 'a'})

    def test_big_traj_is_rejected(self):
        with patch('sys.argv', ['main.py', '-big_traj']):
            with self.assertRaises(SystemExit):
                read_arguments()


if __name__ == '__main__':
    unittest.main()
