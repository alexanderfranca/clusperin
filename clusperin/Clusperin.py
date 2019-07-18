import os
import sys
from collections import defaultdict
import datetime
import csv
import glob
import subprocess
import logging
import logging.handlers

class Clusperin():
    """
    Group the BLAST result similarities into cluster files.
    """

    def __init__( self, source_directory=None, destination_directory=None, cutoff=None, clustering_label=None, author=None, log_file=None ):

        # This class is about generating this dictionary and files from it.
        self.result = defaultdict(dict)

        # Keep tracking of what fasta file is being clustered.
        self.current_fasta_file = None 

        # Log system
        self.log = None

        # Mandatory parameters.
        # There's no meaning for clustering if those parameters are not set.
        self.source_directory = source_directory
        self.destination_directory = destination_directory
        self.cutoff = cutoff
        self.clustering_label = clustering_label
        self.log_file = log_file
        self.author = author

        self.total_of_files = None


    def execute_analysis(self):
        """
        Executes the analysis.

        This is the main method that call all other auxiliary clustering methods.
        """

        self.create_log_system(self.log_file)

        # Metadata comes from a file to store who and how this clustering is being executed.
        # No, that thing cannot be transient because AnEnDB relies on that data.
        # And Clusperin exists only because AnEnDB and for it.
        # AnEnDB have to know who and how clustered its proteins.
        if not self.metadata_exists():
            self.write_metada(self.clustering_label, self.author, self.cutoff)

        self.log.info("START ANALYSIS.")

        self.log.info("EC files are in: " + self.source_directory)

        if not os.path.isdir(self.destination_directory):
            os.mkdir(self.destination_directory)

        self.generate_clusters()

        self.log.info("DONE ANALYSIS.")


    def generate_clusters(self):
        """
        Walk through the directory files and call the 'galpering_analysis' method.

        It also check if the file to be processed was already processed.

        """

        done_files = self.done_files_list()

        files = glob.glob( self.source_directory + '/' + '*.blastall' )

        self.total_of_files = len(files)

        self.log.info("Total of " + str( self.total_of_files ) + " will be processed.")

        counter = 1

        for fasta_file in files:

            self.log.info( "Processing " + str(counter) + " of " + str(self.total_of_files) + " total files." )
            counter += 1

            fasta_file_name = os.path.basename( fasta_file )

            # Only execute clustering if it wasn't run (and done) before.
            if not fasta_file_name in done_files:
                self.current_fasta_file = fasta_file_name

                self.log.info("Going to cluster: " + fasta_file_name)

                self.galperin_analysis(fasta_file)

            else:
                self.log.info("File SKIPED. It was already clustered: " + fasta_file)


    def galperin_analysis(self, fasta_file=None):
        """
        Read the BLAST result files and group its results.

        """

        self.log.info("Start reading blast result: " + fasta_file)
        
        blast_result = self.read_blast_result(fasta_file)

        self.log.info("End of reading blast result: " + fasta_file)


        self.log.info("Start clustering for: " + fasta_file)
        
        self.cluster_proteins(blast_result)

        self.log.info("End of clustering for: " + fasta_file)


    def read_blast_result( self, fasta_file=None ):
        """
        Read BLAST result files and fill the dictionary used by this software to cluster the proteins.

        Args:
            fasta_file(str): Path for the Fasta file to be processed.

        Results:
            (dict): 'result' is the proteins that fit in the specified score/cutoff.

        """

        blast_result_file_name = fasta_file

        blast_csv = open(blast_result_file_name)
        blastall_result_file_csv = csv.reader( blast_csv , delimiter='\t' )

        self.log.info("Reading BLAST results: " + fasta_file)

        for line in blastall_result_file_csv:

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
            systemScore = self.cutoff
            if score >= int(systemScore):
                self.result[query][subject] = 1

                # line added bellow: that's the same effect produced by the balanceHash method, but without the long and overheading 'ifs'.
                self.result[subject][query] = 1

        self.log.info("DONE reading BLAST results: " + fasta_file)

        return self.result


    def cluster_proteins( self, blast_result=None ):
        """
        Actual group/cluster the proteins, write the data into the cluster file and mark the Fasta file as done.

        Args:
            blast_result(dict): The proteins result to be grouped.

        """

        cluster_number = 0
        protein_stack = list()
        clusters = {}
        final_clusters = {}

        if cluster_number == None:
            cluster_number = 0

        for key in blast_result.keys():

            # if there's a value equals to the key name being iterated
            if key in blast_result[key].keys():
                # starts/increment the cluster id
                cluster_number = cluster_number + 1

                # generate/increments the stack of proteins belonging the key
                protein_stack = self.getlist(protein_stack, blast_result[key])

                # remove the whole content of the key (the stack of values was generated already, what's the most important)
                # but anyway is crucial to remove the key to avoid running the 'for' forever (yes, it can happen because the whole logic below) 
                del( blast_result[key] )

                # 'element' is the value we use to put in the cluster
                element = protein_stack.pop()

                # when there's no more elements in the stack we finally reach the end of a cluster
                while element:
                    clusters[element] = cluster_number
                   
                    # if true there's a relation (cluster) between element and another key in the graph 'blastResult'
                    if element in blast_result[element].keys():

                        # increments stack with the values of the related protein
                        protein_stack = self.getlist( protein_stack, blast_result[element] )

                        # keep removing main related protein (already put in the stack above) from the graph to avoid infinite running
                        del( blast_result[element] )

                    # If there's any value in the stack.
                    # If yes, keep droping the last and testing like above.
                    # If not, go back to the main 'for' loop, and finally we reach a different cluster.
                    if protein_stack:
                        element = protein_stack.pop()
                    else:
                        break

        cluster_indexes = list()
        for cluster_value in clusters:
            cluster_indexes.append( str(clusters[cluster_value]) )

        # remove duplications
        cluster_indexes = set(cluster_indexes)

        # fill the finalClusters with lists indexed by the 'clusters' values. 
        for index in cluster_indexes:
            final_clusters[int(index)] = list()

        # fill the finalClusters with the protein identifications and its cluster index
        for cluster_value in clusters:
            final_clusters[ clusters[cluster_value] ].append( cluster_value )


        self.log.info("Generating cluster files: " + self.current_fasta_file)

        # Write the cluster files
        for cluster in final_clusters:

            if not os.path.isdir(self.destination_directory):
                os.mkdir(self.destination_directory)

            cluster_file = self.destination_directory + '/' + self.current_fasta_file + '_' + str(cluster)

            self.log.info("Cluster file: " + cluster_file)

            if os.path.exists(cluster_file):
                os.remove(cluster_file)

            f = open(cluster_file, 'a' )

            for protein in final_clusters[ cluster ]:
                f.write( protein + "\n" )

            f.close()

        del final_clusters
        del clusters

        self.log.info("DONE Generating cluster files: " + self.current_fasta_file)

        # Mark the clustering as done.
        self.mark_clustering_done()

        self.log.info("Clustering file marked as DONE: " + self.current_fasta_file)


    def create_log_system(self, log_file=None):
        """
        Set all logger parameters (file log path, for example), output format and set the class property that stores the logging system.

        Log the process is extremely important because we don't know how long it will take.

        """

        log = logging.getLogger('')
        log.setLevel(logging.DEBUG)
        format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(format)
        log.addHandler(ch)

        fh = logging.handlers.RotatingFileHandler( log_file , maxBytes=0, backupCount=0)
        fh.setFormatter(format)
        log.addHandler(fh)

        self.log = log


    # That's a key function for the clustering.
    # It takes the values of an specific key of the graph
    # and append to an stack.
    # That stack is where the clustering will iterate
    # to specify what cluster a protein belongs to.
    def getlist( self, protein_stack=None, result=None ): 
        """
        Append to an stack a list of values from a key of a graph.

        Args:
            protein_stack(list): List of protein identifications.
            result(dict): Result created until this moment.

        Returns:
            (list): A new protein stack that now acomplish the new data from 'result'.

        """

        for key in result.keys():
            if result[key]:
                protein_stack.append(key)

        return protein_stack


    def done_files_list( self ):
        """
        Read the 'done_files' and return its files list.

        Returns:
            (list): List of file names.
        """

        finished = []

        done_files = self.destination_directory + '/' + 'done_files'

        if os.path.exists(done_files):
            f = open(done_files)

            for line in f:
                # Remove newline
                line = line.rstrip('\r\n')
                finished.append( line )

            f.close()

        return finished


    def mark_clustering_done( self ):
        """
        Write the processed Fasta file into the 'done_files'.

        This file keep tracking of what was already done.

        """

        done_file = self.destination_directory + '/' + 'done_files'
        f = open( done_file, 'a' )
        f.write( self.current_fasta_file + "\n" )
        f.close()


    def write_metada(self, label=None, author=None, cutoff=None):
        """
        This is one of the most important methods.

        Metadata file have to be written since that data will be used as a mark for future relational database insertions.

        Without metadata the relational database wouldn't be able to know what clustering method was used.

        """

        metadata_file = self.destination_directory + '/' + 'clustering_metadata'

        if not os.path.isdir(self.destination_directory):
            os.mkdir(self.destination_directory)

        f = open( metadata_file, 'w' )

        f.write( 'label = ' + str(label) + "\n" )
        f.write( 'author = ' + str(author) + "\n" )
        f.write( 'cutoff = ' + str(self.cutoff) + "\n" )
        f.write( 'clusters = ' + str(self.destination_directory) + "\n" )
        f.write( 'fasta_origin = ' + str(self.source_directory) + "\n" )
        f.write( 'software = blast' + "\n" )
        f.write( 'date = ' + str( datetime.datetime.now() ) + "\n" )

        f.close()


    def metadata_exists(self):
        """
        Check if the clustering metadata file exists.

        """

        if os.path.exists(self.destination_directory + '/' + 'clustering_metadata' ):
            return True

