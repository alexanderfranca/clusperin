import re
import os
import sys
from collections import defaultdict
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


class Clusperin():
    """
    Check configurations, required softwares and run the clustering of EC numbers Fasta files.
    """


    def __init__( self ):
        # This class is about generating this dictionary and files from it.
        self.result = defaultdict(dict)
        #clusters = defaultdict(dict)
        #finalClusters = defaultdict(dict)
        #proteinStack = list()

        # -------------------------------------------------------------------------------- #
        # This class is completely useless if there's no correct configurations set.       #
        # And it's correct and best not to run without it.                                 #
        # -------------------------------------------------------------------------------- #
        # Configuration options system.
        self.config = Config()
        self.config.loadConfiguration()
        self.conf = self.config.getConfigurations()

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



    def getDestination( self ):
        """
        Returns the configuration that tells where to store result files and where is the EC number Fasta files to be processed.

        Returns:
            (str): Path where this class does its job.

        """

        return self.getConfiguration( 'clustering', 'ec_files' )


    def getTrackingFilePath( self ):
        """
        Returns the path where is the clusteringTracking file.

        In other words, the file that stores all the EC numbers files that this class should cluster.

        Returns:
            (str): Full path for the tracking file.
        """

        destination = self.getDestination()

        return destination + '/' + 'clusteringTracking'



    def getFileNameFromList( self, name=None ):
        """
        Remove the path from a full path file name.

        It have the same effect as using 'basename' shell command.

        Args:
            name(str): A full file name path.

        Returns:
            (str): The file name from the path.

        """

        return self.afs.getFileName( name )


    def getEcNumberFromFileName( self, file_name=None ):
        """
        Returns the flat EC number from a Fasta file name.

        Args:
            file_name(str): A Fasta file name. Ex: EC_1.2.34.3.fasta.

        Returns:
            (str): The EC number without any prefix or sufix. For the example above, the result would be: '1.2.34.3'.

        """

        ecNumber = re.sub( '^EC_', '', file_name )
        ecNumber = re.sub( '\.fasta', '', ecNumber )
        ecNumber = str(ecNumber) 

        return ecNumber


    def createTrackingFile( self ):
        """
        Create, if it doesn't exist, the tracking file for the clustering.

        This file actual isn't useful unless you need to know what is the list of 

        EC numbers that should be clustered.

        There's no process in this class that really depends on this tracking file.

        But... that's important to keep the list of what should be done by Clusperin.

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


    def blastProteins( self, ec_file=None ):
        """
        Execute BLAST software to generate the similarity results.

        Also, and most important, create the dictionary with the hits that fits in the cutoff/score parameter.

        Args:
            ec_file(str): Path for the EC file to be processed.

        Results:
            (dict): 'result' is the proteins that fit in the specified score/cutoff.

        """

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
            systemScore = self.getConfiguration( 'clustering', 'cutoff' )
            if score >= systemScore:
                self.result[query][subject] = 1

                # line added bellow: that's the same effect produced by the balanceHash method, but without the long and overheading 'ifs'.
                self.result[subject][query] = 1

        self.log.info( "DONE reading BLAST results: " + ecFile )

        return self.result


    # That's a key function for the clustering.
    # It takes the values of an specific key of the graph
    # and append to an stack.
    # That stack is where the clustering will iterate
    # to specify what cluster a protein belongs to.
    def getlist( self, proteinStack=None, result=None ): 
        """
        Append to an stack a list of values from a key of a graph.

        Args:
            proteinStack(list): List of protein identifications.
            result(dict): Result created until this moment.

        Returns:
            (list): A new protein stack that now acomplish the new data from 'result'.

        """

        for key in result.keys():
            if result[key]:
                proteinStack.append(key)

        return proteinStack


    def clusterProteins( self, blastResult=None ):
        """
        Actual group/cluster the proteins, write the data into the cluster file and mark the EC file as done.

        Args:
            blastResult(dict): The proteins result to be grouped.

        """

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



    def galperinAnalysis(self, ec_file=None ):
        """
        Call the BLAST software and call the actual clustering method.

        """

        self.log.info( "Start blasting process for: " + ec_file )
        
        blastResult = self.blastProteins( ec_file )

        self.log.info( "End of blasting process for: " + ec_file )


        self.log.info( "Start clustering for: " + ec_file )
        
        self.clusterProteins( blastResult )

        self.log.info( "End of clustering for: " + ec_file )


    def generateClusters( self ):
        """
        Walk through the directory files and call the 'galperingAnalysis' method.

        It also check if the file to be processed was already processed.

        """

        trackingFile = self.getTrackingFilePath()

        self.log.info( "Tracking file is: " + trackingFile )

        doneFiles = self.getDoneFilesList()

        filesDirectory = self.getDestination()

        files = glob.glob( filesDirectory + '/' + '*.fasta' )

        self.totalOfFiles = len(files)

        self.log.info( "Total of " + str( self.totalOfFiles ) + " will be processed." )

        counter = 1

        for ecFile in files:

            self.log.info( "Processing " + str(counter) + " of " + str(self.totalOfFiles) + " total files." )
            counter += 1

            ecFileName = self.afs.getFileName( ecFile )

            # Only execute clustering if it wasn't run (and done) before.
            if not ecFileName in doneFiles:
                self.currentEcFile = self.afs.getFileName( ecFile )

                self.log.info( "Going to cluster: " + ecFile )

                self.galperinAnalysis( ecFile )
            else:
                self.log.info( "File SKIPED. It was already clustered: " + ecFile )


    def executeAnalysis( self ):
        """
        Executes the analysis.

        This is the main method that call all other auxiliary clustering methods.
        """

        self.createLogSystem()

        self.writeMetada()

        self.log.info( "START ANALYSIS." )

        clusterDestination = self.getDestination()

        self.log.info( "EC files are in: " + clusterDestination )

        if not self.afs.isDirectory( clusterDestination ):
            os.mkdir( clusterDestination )

        self.generateClusters()

        self.log.info( "DONE ANALYSIS." )


    def getDoneFilesList( self ):
        """
        Read the 'done_files' and return it files list.

        Returns:
            (list): List of file names.
        """

        finished = []

        clusterDestination = self.getConfiguration( 'clustering', 'ec_files' )

        doneFiles = clusterDestination + '/' + 'done_files'

        if self.afs.fileExists( doneFiles ):
            f = open( doneFiles )

            for line in f:
                # Remove newline
                line = line.rstrip('\r\n')
                finished.append( line )

            f.close()

        return finished


    def markClusteringDone( self ):
        """
        Write the processed EC number into the 'done_files'.

        This file keep tracking of what was already done.

        """

        doneFile = self.getDestination() + '/' + 'done_files'
        f = open( doneFile, 'a' )
        f.write( self.currentEcFile + "\n" )
        f.close()


    def purgeTemporaryFiles( self, ec_file=None ):
        """
        Remove temporary BLAST files. Files that ends with .phr, .pin, .psq and .blastall.

        This method is used when a clustering is interrupted for some reason and old processed files 

        have to be removed: We have to make sure a whole clustering process was executed neat and clean. 
        """

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


    def purgeClusterFile( self, file_to_purge=None ):
        """
        Remove file.
        """

        self.afs.removeFile( file_to_purge )


    def writeMetada( self ):
        """
        This is one of the most important methods.

        Metadata file have to be written since that data will be used as a mark for future relational database insertions.

        Without metadata the relational database wouldn't be able to know what clustering method was used.

        """

        metadataFile = self.getDestination() + '/' + 'clustering_metadata'

        f = open( metadataFile, 'w' )

        label = self.getConfiguration( 'clustering', 'label' )

        f.write( 'label = ' + str(label) + "\n" )
        f.write( 'date = ' + str( datetime.datetime.now() ) + "\n" )

        if self.conf.has_option( 'clustering', 'author' ):
            author = self.getConfiguration( 'clustering', 'author' )
        else:
            author = 'anonymous'

        f.write( 'author = ' + str(author) + "\n" )

        f.close()

