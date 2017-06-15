import sys
import os
sys.path.insert(0,  os.getcwd() + '/../')
import unittest
from Clusperin import *
import re


# TODO: This class have to be in fact implemented.
class TestClusperin( unittest.TestCase ):

    def setUp( self ):
        self.cl = Clusperin()


    def test_getConfiguration( self ):

        self.cl.setConfigurationFile( './fixtures/clusperin.conf' )
        expectedValue = '/var/kegg/clustering'

        result = self.cl.getConfiguration( 'clustering', 'ec_files' )

        self.assertEquals( expectedValue, result )


#    def test_createTrackingFile( self ):
#
#        self.cl.setConfigurationFile( './fixtures/anendb.conf' )
#
#        self.cl.createTrackingFile()
#
#
#    def test_executeAnalysis( self ):
#
#        self.cl.setConfigurationFile( './fixtures/anendb.conf' )
#
#        self.cl.executeAnalysis()
#










if __name__ == "__main__":
    unittest.main()
