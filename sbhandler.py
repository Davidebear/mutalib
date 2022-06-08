# ----------------------------------------------------------------------------------
# Functions to handle the unique attributes of the consensus1.mat file in the Gong Lab
# @author: davidebear
# 
# Note: Refer to https://github.com/emiliosalazar/matToPython/blob/main/MatFileMethods.py for a more encompassing MatLabLoader
# ----------------------------------------------------------------------------------

# Needed imports
import numpy as np
from math import pow
import h5py

def load_matlab_file(mat_file, variable_name): # Note the double return statement
    """_summary_
    //UNIQUE: Only one variable for this .mat file
    
    Args:
        mat_file (string): Name of .mat file to intake (path must be given if not in current dir)
        variable_name (string): which variable from the mat file. Only one here because of the .mat file structure.

    Returns:
        h5py_object, data: Overall h5py_object and the specific variable as an nparray
    """
    h5py_object = h5py.File(mat_file, 'r')
    data = h5py_object.get(variable_name)
    data = np.array(data)
    return h5py_object, data

class DNA_SeqBlocks(): # Loaded mat file that serves as a generator of seqblocks and ultimately CleanSeqBlock objects
    """_summary_
    Loads the mat file with DNA sequences for easy iteration through all stored DNA sequence lists of varying coverage
    """
    def __init__(self, h5py_object, data):
        """_summary_
        Args:
            h5py_object (h5py_object): first return variable of load_matlab_file
            data (np.array): second return variable of load_matlab_file
        """
        self.h5py_object = h5py_object
        self.data = np.transpose(data)
        self.size = len(data) 
    def get_seqblock(self, number): # //TODO: make sure the default args work
        """_summary_
        Generates a seqblock and transposes the result. mat files store matrices in their transpose form

        Args:
            number (int): Which of the over 300,000 (self.size) seqblocks do you want? Input with MATLAB indexing

        Returns:
            np.array object: An unclean seqblock with dtype='uint8'. Further procesing to be done in seqblock_parser and into CleanSeqBlock object
        """
        x = self.h5py_object[self.data[number-1, 0]][:, :] # //UNIQUE: 0th column because .mat file, align3 variable only has one row
        return np.transpose(x)
    def get_coverage_count(self, number):
        x = self.h5py_object[self.data[0, number-1]][0,:]
        x = int((len(x)-4)/2)
        return x
    
#//TODO: Create a container of all seqblocks for readability (char arrays). Change seqblock_parser to an int4 array.
# Do I need to make each of the four bases a number like [0.00, 0.33, 0.66, 0.99] to normalize everything?
def seqblock_parser(seqblock):
    """_summary_

    Args:
        seqblock (DNA_SeqBlocks.get_seqblock(number=)): A seqblock generated from the DNA_SeqBlocks dataset. Chooses a specified seqblock

    Returns:
        CleanSeqBlock object: Contains char arrays of cleaned aspects of each seqblock (target sequence, reads, quality, interpretations (3 types)).
        Useable data for training, validation, testing.
    """
    seqblock_parsed = CleanSeqBlock()
    total_columns = len(np.transpose(seqblock)) # is there a more efficient way?
    total_rows = len(seqblock)
    number_of_reads = int((total_rows - 4)/2) # last three rows have intreptations, first row is non-mutated target sequence, and N quality scores for N reads
    
    seqblock_parsed.reads_count = number_of_reads;
    seqblock_parsed.len = total_columns;
    
    for i in range(0,total_rows):
        # print(f" The sequenceblock input {seqblock}") #debug
        # print(f" The m x n size of the seqblock {seqblock.shape}") #debug
        # print(f" The ith row of the seqblock {current_seq}") #debug
        # print(f" The m x n size of the sequence {current_seq.shape}") #debug
        new_format = ""
        
        # print(f" The 0th entry of the initialized char array {new_format[0, 0]}") #debug
        # print(f" Its row size {len(new_format[0])}")
        for j in range(total_columns):
            if seqblock[i,j] == 32: # //Confirmed
                new_format += " "
            else: new_format += chr(seqblock[i,j])
            
         # Neatly divides data into attributes of CleanSeqBlock object based on row number.
        if (i == 0):
            seqblock_parsed.target = new_format
        elif (i == total_rows - 1): # last row is mutations 'x', total_rows is one more than total index
            seqblock_parsed.interp_mutations = new_format
        elif (i == total_rows - 2): # 2nd last is subjective consensus
            seqblock_parsed.interp_consensus = new_format
        elif (i == total_rows - 3): 
            seqblock_parsed.interp_changes = new_format
        elif (i > 0 and i < number_of_reads + 1):
            seqblock_parsed.reads.append(new_format)
        else:
            seqblock_parsed.qscores.append(new_format)
    return seqblock_parsed     

class CleanSeqBlock(): 
    """_summary_
    Houses six key sequence types from mat file and useable, clean data
    Object of this type are initialized by seqblock_parser function.
    """
    def __init__(self):
        self.reads_count = 0;
        self.len = 0;
        self.barcode = 0; # Super important for future dictionary
        
        self.target = 0;
        
        self.reads = [];
        self.qscores = [];
        
        self.interp_changes = 0; # what is this again? 
        self.interp_consensus = 0;
        self.interp_mutations = 0;
    
def pscore(Qscore): 
    return pow(10, -(Qscore-33)/10)