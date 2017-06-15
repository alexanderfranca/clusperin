import re
import os
import sys
from collections import defaultdict
from sqlalchemy.sql import func
import time
import decimal
import datetime
import csv
import glob
import shutil
import subprocess
import logging
import logging.handlers
from Config import *
from AnendbFileSystem import *

result   = defaultdict(dict)
clusters = defaultdict(dict)
finalClusters = defaultdict(dict)
proteinStack = list()


class AnenpiClusteringMethod():

    # This class is about generating this dictionary and files from it.
    result = defaultdict(dict)

    def __init__( self ):

        # Configuration options system.
        self.config = Config()

        # File system helper.
        self.afs    = AnendbFileSystem()

        # Keep tracking of what EC files is being clustered.
        self.currentEcFile = None 

        # Log system
        self.log = None


    def createLogSystem(self):
        """
        Set all logger parameters (file log path, for example), output format and set the class property that stores the logging system.

        """

        logFile  = self.getConfiguration( 'log', 'log_file' )

        log = logging.getLogger('')
        log.setLevel(logging.DEBUG)
        format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(format)
        log.addHandler(ch)

        fh = logging.handlers.RotatingFileHandler( logFile , maxBytes=0, backupCount=0)
        fh.setFormatter(format)
        log.addHandler(fh)

        self.log = log


    def printConfigurationFileExample( self ):

        print( "\n" )
        print( "Here's an example of .anendb.conf file." )
        print( "You should put this file in your home directory. ")
        print( "If your home directory is /home/claypool , you create your file as: " )
        print( "/home/claypool/.anendb.conf" )
        print( "-------------------------------------------------------------------------" )
        print( "[clustering]" )
        print( "ec_files = /var/kegg/clustering/" )
        print( "cluster_files = /var/kegg/clustering/clusters/" )
        print( "" )
        print( "[log]" )
        print( "log_file = /var/kegg/clustering/clustering.log" )
        print( "-------------------------------------------------------------------------" )


    def isConfigurationsCorrect( self ):
        """
        Doesn't allow the clustering process without the correct configurations available.

        """

        errors = 0

        # Check configuration file.
        expectedClusteringOptions  = [ 'ec_files', 'cluster_files' ]
        expectedLogOptions = [ 'log_file' ]

        if not self.conf.has_section( 'clustering' ):
            print( "ERROR: The anendb.conf section 'clustering' wasn't found!" )
            self.printConfigurationFileExample()
            sys.exit()

        if not self.conf.has_section( 'log' ):
            print( "ERROR: The anendb.conf section 'log' wasn't found!" )
            self.printConfigurationFileExample()
            sys.exit()


        for option in expectedClusteringOptions:
            if not self.conf.has_option( 'clustering', option ):
                print( "ERROR: The anendb.conf option '" + str(option) + "' in the 'clustering' section wasn't found!" )
                errors += 1

        for option in expectedLogOptions:
            if not self.conf.has_option( 'log', option ):
                print( "ERROR: The anendb.conf option '" + str(option) + "' in the 'log' section wasn't found!" )
                errors += 1




        # Check blast exists.
        command = 'makeblastdb -help'

        result = subprocess.call( command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Zero means (from the subprocess package) the command could be executed without errors.
        if not result == 0:
            pprint.pprint( result )
            print( "Software 'makeblastdb' wasn't found. Probably you don't have BLAST installed." )
            print( "Install BLAST before running this clustering process.")
            sys.exit()

        command = 'blastp -help'

        result = subprocess.call( command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Zero means (from the subprocess package) the command could be executed without errors.
        if not result == 0:
            print( "Software 'blastp' wasn't found. Probably you don't have BLAST installed." )
            print( "Install BLAST before running this clustering process.")
            sys.exit()


        if errors > 0:
            self.printConfigurationFileExample()
            sys.exit()
        else:
            return True



    # TODO: test, comment
    def getDestination( self ):

        return self.getConfiguration( 'clustering', 'ec_files' )


    # TODO: test, comment
    def getTrackingFilePath( self ):

        destination = self.getDestination()

        return destination + '/' + 'clusteringTracking'



    # TODO: test, comment
    def getFileNameFromList( self, name=None ):

        return self.afs.getFileName( name )


    # TODO: test, comment
    def getEcNumberFromFileName( self, file_name=None ):

        ecNumber = re.sub( '^EC_', '', file_name )
        ecNumber = re.sub( '\.fasta', '', ecNumber )
        ecNumber = str(ecNumber) 

        return ecNumber


    # TODO: test, comment
    def createTrackingFile( self ):
        """
        Create, if it doesn't exist, the tracking file for the clustering.

        """

        trackingFile = self.getTrackingFilePath() 

        # Create the file only if it doesn't exist.
        if not self.afs.fileExists( trackingFile ):

            f = open( trackingFile, 'a' ) 

            files = glob.glob( destination + '/' + '*.fasta' )

            for ecFile in files:
                fileName = self.getFileNameFromList( ecFile )

                ecNumber = self.getEcNumberFromFileName( fileName )

                f.write( ecNumber + "\n" )

            f.close()


    # TODO: test, comment
    def getTrackingFile( self ):

        destination 

    def setConfigurationFile( self, conf_file=None ):
        """
        Set a new configuration file to be used by Anendb.

        Args:
            conf_file(str): Full path for the configuration file.

        """

        self.config.setConfigurationFile( conf_file )
        self.config.loadConfiguration()
        self.conf = self.config.getConfigurations()


    def getConfiguration( self, section=None, value=None ):
        """
        Return configuration from anendb.conf file.

        Args:
            section(str): anendb.conf section (Ex. [kegg2017]).
            value(str): anendb.conf value (Ex. username).

        Returns:
            (str): Value from the section in the anendb.conf file.

        """

        return self.conf.get( section, value )



    def setTrackingFile(self, trackingFile):
        pass

    def setCutoff( self, cutoff ):
        pass


    def getEcBlastResultFileName( self, ecNumber ):
        """
        Concatenates an string using the EC Number and
        the extension '.blastall'
        """
        return str(ecNumber) + ".blastall"


    def getBlastallResultFile( self, ecNumber):
        """
        Returns a file handle representing
        the file where the results of blast
        command are stored.    
        """

        fileName = self.getEcBlastResultFileName( ecNumber )
        blastAllFile = open( self.processingFilesDirectory + '/' + fileName, 'r')

        return blastAllFile


    def similaritiesPgsqlInserts(self, proteinDataList, ecNumber):
        pass


    # TODO: test, comment
    def blastProteins( self, ec_file=None ):

        ecFile = ec_file

        destination = self.getDestination()

        self.purgeTemporaryFiles( ec_file )

        ecBlastResultFileName = ec_file + '.blastall'

        self.log.info( "Creating BLAST DB for: " + ecFile )

        os.popen("makeblastdb -in " + ecFile + " -out " + ecFile + " -dbtype prot")

        self.log.info( "DONE creating BLAST DB for: " + ecFile )

        self.log.info( "BLAST PROTEINS! : " + ecFile )

        os.popen("blastp -query " + ecFile + " -outfmt '6 qseqid sseqid pident ppos score bitscore qstart qend sstart send qlen slen evalue' -evalue 0.1 -num_alignments 10000000 -db " + ecFile + " -out " + ecBlastResultFileName)

        self.log.info( "DONE BLAST PROTEINS! : " + ecFile )

        blastCsv = open( ecBlastResultFileName )
        blastallResultFileCSV = csv.reader( blastCsv , delimiter='\t' )


        self.log.info( "Reading BLAST results: " + ecFile )

        for line in blastallResultFileCSV:

            resultQuery                              = line[0]
            resultSubject                            = line[1]
            resultPercentageIdenticalMatches         = line[2]
            resultPercentageOfPositiveScoringMatches = line[3]
            resultScore                              = line[4]
            resultBitScore                           = line[5]
            resultStartOfAlignmentInQuery            = line[6]
            resultEndOfAlignmentInQuery              = line[7]
            resultStartOfAlignmentInSubject          = line[8]
            resultEndOfAlignmentInSubject            = line[9]
            resultQueryLength                        = line[10]
            resultSubjectLength                      = line[11]
            resultEValue                             = line[12]

            query   = resultQuery
            subject = resultSubject

            score = float( resultScore )
            score = int( score )

            # Most critical element of the whole method:
            # point where the score value is actually determined/tested
            if score >= 120:
                result[query][subject] = 1

                # line added bellow: that's the same effect produced by the balanceHash method, but without the long and overheading 'ifs'.
                result[subject][query] = 1

        self.log.info( "DONE reading BLAST results: " + ecFile )

        return result


    # That's a key function for the clusterization.
    # It takes the values of an specific key of the graph
    # and append to an stack.
    # That stack is where the clusterization will iterate
    # to specify what cluster a protein belongs to.
    def getlist( self, proteinStack=None, result=None ): 
        """
        Append to an stack a list of values from a key of a graph.
        **Returns** a list.
        """

        for key in result.keys():
            if result[key]:
                proteinStack.append(key)

        return proteinStack


    # TODO: test, comment
    def clusterProteins( self, blastResult=None ):

        clusterNumber = 0
        proteinStack = list()
        clusters = {}
        finalClusters = {}

        if clusterNumber == None:
            clusterNumber = 0

        for key in blastResult.keys():

            # if there's a value equals to the key name being iterated
            if key in blastResult[key].keys():
                # starts/increment the cluster id
                clusterNumber = clusterNumber + 1

                # generate/increments the stack of proteins belonging the key
                proteinStack = self.getlist(proteinStack, blastResult[key])

                # remove the whole content of the key (the stack of values was generated already, what's the most important)
                # but anyway is crucial to remove the key to avoid running the 'for' forever (yes, it can happen because the whole logic below) 
                del( blastResult[key] )

                # 'element' is the value we use to put in the cluster
                element = proteinStack.pop()

                # when there's no more elements in the stack we finally reach the end of a cluster
                while element:
                    clusters[element] = clusterNumber
                   
                    # if true there's a relation (cluster) between element and another key in the graph 'blastResult'
                    if element in blastResult[element].keys():

                        # increments stack with the values of the related protein
                        proteinStack = self.getlist( proteinStack, blastResult[element] )

                        # keep removing main related protein (already put in the stack above) from the graph to avoid infinite running
                        del( blastResult[element] )

                    # If there's any value in the stack.
                    # If yes, keep droping the last and testing like above.
                    # If not, go back to the main 'for' loop, and finally we reach a different cluster.
                    if proteinStack:
                        element = proteinStack.pop()
                    else:
                        break

        clusterIndexes = list()
        for clusterValue in clusters:
            clusterIndexes.append( str(clusters[clusterValue]) )

        # remove duplications
        clusterIndexes = set(clusterIndexes)

        # fill the finalClusters with lists indexed by the 'clusters' values. 
        for index in clusterIndexes:
            finalClusters[int(index)] = list()

        # fill the finalClusters with the protein identifications and its cluster index
        for clusterValue in clusters:
            finalClusters[ clusters[clusterValue] ].append( clusterValue )


        self.log.info( "Generating cluster files: " + self.currentEcFile )

        # Write the cluster files
        for cluster in finalClusters:

            clusterDestination = self.getConfiguration( 'clustering', 'cluster_files' )

            if not self.afs.isDirectory( clusterDestination ):
                os.mkdir( clusterDestination )

            clusterFile = clusterDestination + '/' + self.currentEcFile + '_' + str(cluster)

            self.log.info( "Cluster file: " + clusterFile )

            if self.afs.fileExists( clusterFile ):
                self.afs.removeFile( clusterFile )

            f = open( clusterFile, 'a' )

            for protein in finalClusters[ cluster ]:
                f.write( protein + "\n" )

            f.close()

        del finalClusters
        del clusters

        self.log.info( "DONE Generating cluster files: " + self.currentEcFile )

        # Mark the clustering as done.

        self.markClusteringDone()

        self.log.info( "Clustering file marked as DONE: " + self.currentEcFile )



    # TODO: test, comment
    def galperinAnalysis(self, ec_file=None ):

        self.log.info( "Start blasting process for: " + ec_file )
        
        blastResult = self.blastProteins( ec_file )

        self.log.info( "End of blasting process for: " + ec_file )


        self.log.info( "Start clustering for: " + ec_file )
        
        self.clusterProteins( blastResult )

        self.log.info( "End of clustering for: " + ec_file )


    def generateClusters( self ):
        """
        Executes the analysis
        """

        trackingFile = self.getTrackingFilePath()

        self.log.info( "Tracking file is: " + trackingFile )

        doneFiles = self.getDoneFilesList()

        filesDirectory = self.getDestination()

        files = glob.glob( filesDirectory + '/' + '*.fasta' )

        for ecFile in files:

            # Only execute clustering if it wasn't run (and done) before.
            if not ecFile in doneFiles:
                self.currentEcFile = self.afs.getFileName( ecFile )

                self.log.info( "Going to cluster: " + ecFile )

                self.galperinAnalysis( ecFile )
            else:
                self.log.info( "File was already clustered: " + ecFile )


    def executeAnalysis(self, fromScratch=None, ec=None):
        """
        Executes the analysis
        """


        if self.isConfigurationsCorrect():

            self.createLogSystem()

            self.log.info( "START ANALYSIS." )

            clusterDestination = self.getDestination()

            self.log.info( "EC files are in: " + clusterDestination )

            if not self.afs.isDirectory( clusterDestination ):
                os.mkdir( clusterDestination )

            self.generateClusters()

        self.log.info( "DONE ANALYSIS." )


    def getDoneFilesList( self ):

        finished = []

        clusterDestination = self.getConfiguration( 'clustering', 'cluster_files' )

        doneFiles = clusterDestination + '/' + 'done_files'

        if self.afs.fileExists( doneFiles ):
            f = open( doneFiles )

            for line in f:
                finished.append( line )

            f.close()

        return finished


    def markClusteringDone( self ):

        doneFile = self.getDestination() + '/' + 'done_files'
        f = open( doneFile, 'a' )
        f.write( self.currentEcFile + "\n" )
        f.close()


    def purgeTemporaryFiles( self, ec_file=None ):

        destinationDirectory = self.getConfiguration( 'clustering', 'cluster_files' )

        fileMask = ec_file

        files = []

        files.append( fileMask + '.phr'      )
        files.append( fileMask + '.pin'      )
        files.append( fileMask + '.psq'      )
        files.append( fileMask + '.blastall' )

        for fileToRemove in files:
            if self.afs.fileExists( fileToRemove ):
                self.log.info( "Purging old file: " + fileToRemove )
                self.afs.removeFile( fileToRemove )


    # TODO: test, comment
    def purgeClusterFile( self, file_to_purge=None ):

        self.afs.removeFile( file_to_purge )




