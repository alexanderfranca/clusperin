import sys
import os
sys.path.insert(0,  os.getcwd() + '/../')
import unittest
from Clusperin import *
import re


# TODO: This class have to be in fact implemented.
class TestClusperin(unittest.TestCase):

    def setUp( self ):
        self.cl = Clusperin(log_file='/var/log/cluspering.log', 
                destination_directory='/var/kegg/clustering/clusters/',
                source_directory='/var/kegg/clustering/',
                cutoff=120,
                )


    def test_execute_analysis( self ):

        self.cl.execute_analysis()











if __name__ == "__main__":
    unittest.main()
