import unittest
import os
import lifto
import shlex


class TestLifto(unittest.TestCase):
    # run tests from parent directory `$python -m unittest test.lifto_tests`
    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.test_path = "test/data_tests/"
        self.input = f"{self.test_path}hg19_data_test.csv"
        self.output = "out/hg19_to_hg38_from_test_data.csv"

    def test_basic(self):
        print(os.getcwd())
        lifto.lifto(self.input, output_type="csv")
        os.path.isfile(self.output)

    def test_only_chrompos_data(self):
        input = f"{self.test_path}hg19_onlychrompos_test.csv"
        result = lifto.lifto(input, output_type="dict")
        self.assertIsNotNone(result)
        self.assertEqual(len(list(result.keys())[1]), 2)

    def test_argparse(self):
        args = f'--in_file "{self.input}" --output_type csv --show_err True'
        lifto.main(shlex.split(args))


if __name__ == "__main__":
    unittest.main()
