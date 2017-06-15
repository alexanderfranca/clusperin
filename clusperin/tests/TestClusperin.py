import sys
import os
sys.path.insert(0,  os.getcwd() + '/../')
import unittest
from AnenpiClusteringMethod import *
import re


class TestAnenpiClusteringMethod( unittest.TestCase ):

    def setUp( self ):
        self.cl = AnenpiClusteringMethod()


#    def test_getConfiguration( self ):
#
#        self.cl.setConfigurationFile( './fixtures/anendb.conf' )
#        expectedPass = 'darkmatter'
#
#        result = self.cl.getConfiguration( 'kegg2017', 'pass' )
#
#        self.assertEquals( expectedPass, result )
#
#
#    def test_createTrackingFile( self ):
#
#        self.cl.setConfigurationFile( './fixtures/anendb.conf' )
#
#        self.cl.createTrackingFile()
#
#
    def test_executeAnalysis( self ):

        self.cl.setConfigurationFile( './fixtures/anendb.conf' )

        self.cl.executeAnalysis()











if __name__ == "__main__":
    unittest.main()
