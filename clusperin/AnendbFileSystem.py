#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import pprint

class AnendbFileSystem:

    def fileExists( self, file_path=None ):
        """
        Test if the file path exists.
    
        Args:
            file_path(str): Path of the file.

        Returns:
            (boolean)
        """

        try:
            if os.path.exists( file_path ):
                return True 
            else:
                return False

        except:
            print( 'Could not test the file' )


    def directoryExists( self, directory_path=None ):
        """
        Test if the directory exists.

        Args:
            directory_path(str): Path of the directory.

        Returns:
            (boolean)
        """

        try:
            if os.path.exists( directory_path ):
                return True
            else:
                return False

        except:
            print( 'Could not test the directory.' )


    def isDirectory( self, directory_path ):
        """
        Test if the path is a valid directory.

        Args:
            directory_path(str): Path of the directory.

        Returns:
            (boolean)
        """

        try:
            if os.path.isdir( directory_path ):
                return True
            else:
                return False

        except:
            print( 'Coult not test if directory is valid.' )


    def isFile( self, file_path=None ):
        """
        Test if the path is actual a file.

        Args:
            file_path(str): File path.

        Returns:
            (boolean)
        """

        try:
            if os.path.isfile( file_path ):
                return True
            else:
                return False

        except:
            print( "Could not test the file." )

    def isValidFile( self, file_path=None ):
        """
        Validates if the file exists and it's a valid file.

        Args:
            file_path(str): File path.

        Returns:
            (boolean)
        """

        if self.isFile( file_path ) and self.fileExists( file_path ):
            return True
        else:
            raise IOError("File " + file_path + " doesn't exist, you don't have permission or not even it's a file.")


    def isValidDirectory( self, directory_path=None ):
        """
        Validates if the directory exists and it's a valid directory.

        Args:
            directory_path(str): The directory path.

        Returns:
            (boolean)
        """


        if self.isDirectory( directory_path ) and self.directoryExists( directory_path ):
            return True
        else:
            return False



    def openFile( self, file_path=None, mode='r'):
        """
        Opens a file an return its handle.

        Args:
            file_path(str): File Path.
            mode(str): File mode ('r' = reading, 'a' = append, 'w' = writing )


        Returns:
            (file): File handle of the opened file.
        """

        try:
            f = open( file_path, mode )
            return f

        except:
            print( 'Could not open the file: ' + file_path )



    def openFileForReading( self, file_path=None ):
        """
        Opens a file for reading.

        Args:
            file_path(str): File path.

        Returns:
            (file): File handle.
        """

        if self.isValidFile( file_path ):
            return self.openFile( file_path, 'r' )
        else:
            raise IOError( "File " + file_path + " couldn't be opened." )


    def openFileForWriting( self, file_path=None ):
        """
        Opens a file for writing.

        Args:
            file_path(str): File path.

        Returns:
            (file): File handle.
        """

        return self.openFile( file_path, 'w' )


    def openFileForAppend( self, file_path=None ):
        """
        Opens a file for append.

        Args:
            file_path(str): File path.

        Returns:
            (file): File handle.
        """

        return self.openFile( file_path, 'a' )


    def closeFile( self, file_handle=None ):
        """
        Simply closes a file handle.

        Args:
            file_handle(file): File handle to be closed.

        Returns:
            (boolean)

        """

        # make sure the file is really closed. Is that really opened? That kind of stuff.
        try:
            file_handle.close()
            return True

        except:
            print( 'Could not close the file.' )


    def removeFile( self, file_path=None ):
        """
        Remove a file.

        Args:
            file_path(str): File path.

        Returns:
            (boolean)
        """

        if self.isValidFile( file_path ):
            try:
                os.remove( file_path )
                return True
            except:
                print( 'Could not remove the file.' )

        else:
            print( 'Invalid file to be removed.' )
            return False


    def purgeAndCreateFile( self, file_path=None ):
        """
        Create a file but remove previous created (if exists).

        Args:
            file_path(str): File path.

        Returns:
            (file)
        """

        if self.isValidFile( file_path ):
            self.removeFile( file_path )
            f = self.openFileForAppend( file_path ) 
            return f 

        else:
            return False


    def removeDirectory( self, directory_path=None ):
        """
        Remove a directory.

        Args:
            directory_path(str). The directory path.

        Returns:
            (boolean)
        """

        if self.isValidDirectory( directory_path ):
            try:
                shutil.rmtree( directory_path )
                return True
            except:
                print( 'Could not remove the directory: ' + directory_path )
        else:
            return False


    def getDirectoryName( self, path=None ):
        """
        Return the directory path of a file path.

        Args:
            path(str): Full path of some file.

        Returns:
            (str): The directory path without the file name.
        """

        if self.isValidFile( path ):
            return os.path.dirname( path )
        else:
            return None


    def getFileName( self, path=None ):
        """
        Return the name of file from a full path.

        Args:
            path(str): Full path of some file.

        Returns:
            (str): The file name from the path.
        """

        if self.isValidFile( path ):
            return os.path.basename( path )
        else:
            return None

