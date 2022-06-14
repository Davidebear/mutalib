
# from sbhandler import *
from math import pow    
# from fastai.vision.all import *

# Needed imports
import numpy as np
import h5py

STOP_CODON = np.array([84, 65, 65, 84, 65, 71]) # T A A T A G

# Source code for functions from seqblock_handler.py
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
class DNA_SeqBlocks():
    """Mat dataset 
    Loads the mat file with DNA sequences for easy iteration through all stored DNA sequence lists of varying coverage
    """
    def __init__(self, h5py_object, data):
        """Initialize DNA_SeqBlocks() object
        Args:
            h5py_object (h5py_object): first return variable of load_matlab_file
            data (np.array): second return variable of load_matlab_file
        """
        self.h5py_object = h5py_object
        self.data = np.transpose(data)
        self.size = len(data) 
        
    def get_seqblock(self, number): # checked
        x = self.h5py_object[self.data[0, number-1]][:, :] # //UNIQUE: 0th column because .mat file, align3 variable only has one row
        return np.transpose(x)
    
    def get_target(self, number):
        x = self.h5py_object[self.data[0, number-1]][:, 0]
        return np.transpose(x)
    
    def get_coverage_count(self, number):
        """Coverage count

        Args:
            number (int): which sequence? matlab idx

        Returns:
            int: Number of reads
        """
        x = self.h5py_object[self.data[0, number-1]][0,:]
        coverage_count = int((len(x)-4)/2)
        return coverage_count
    
    def get_barcode(self, number): # //TODO: Ask Yiyang about a mutated stop codon.
        """Barcode finder
        
        Changes that need to be made: use target sequence to find stop codon (reads may have mutated stop codon like 22776). Use read that doesn't terminate in indels to find the barcode.

        Args:
            number (int): which sequence? matlab indexing

        Returns:
            barcode, true_len: returns barcode uint16 ndarray, last matlab idx of stop codon
        """
        x = self.h5py_object[self.data[0, number-1]][:, 0]
        barcode = ""
        true_len = 0
        total = x.size
        
        idx = 0
        search_stop=True
        
        # This block is concise but had 2x the runtime
        # while idx < total and search_stop: # To identify the stop codon
        #     if (np.all(x[idx:idx+6] == STOP_CODON)):  # Translates to 'TAATAG'
        #         true_len = idx+6
        #         search_stop = False
        #     idx += 1
        
        while idx < total and search_stop:
            if (x[idx] == 84):
                try:
                    if (x[idx+1] == 65):
                        if x[idx+2] == 65:
                            if x[idx+3] == 84:
                                if x[idx+4] == 65:
                                    if x[idx+5] == 71:
                                        true_len = idx+6
                                        search_stop = False
                except: 
                    barcode = None
                    true_len = 0
                    return barcode, true_len
            idx += 1
        barcode = x[true_len:total]
            
        return barcode, true_len
    
    def get_nuc_qscores(self, number, nuc_pos, coverage_count=0): # checked
        if coverage_count == 0:
            cov = self.get_coverage_count(number)
            x = self.h5py_object[self.data[0, number-1]][nuc_pos-1, (1+cov):(1+2*cov)]
            return np.transpose(x)
        else:
            x = self.h5py_object[self.data[0, number-1]][nuc_pos-1, (1+coverage_count):(1+2*coverage_count)]
            return np.transpose(x)
    
    def get_nuc_reads(self, number, nuc_pos, coverage_count=0): # checked
        if coverage_count == 0:
            x = (self.h5py_object[self.data[0,number-1]][nuc_pos-1, 1:self.get_coverage_count(number)+1])
            return np.transpose(x)
        else: 
            x = self.h5py_object[self.data[0,number-1]][nuc_pos-1, 1:(coverage_count+1)]
            return np.transpose(x)
    
    def get_reads(self, number, coverage_count=0, true_len=0): # CHECK THE INDEXING
        if coverage_count == 0 and true_len == 0:
            _, true = self.get_barcode(number)
            x = (self.h5py_object[self.data[0,number-1]][0:true, 1:self.get_coverage_count(number)+1])
            return np.transpose(x)
        else: 
            x = self.h5py_object[self.data[0,number-1]][0:true_len, 1:(coverage_count+1)]
            return np.transpose(x)
    
    def get_qscores(self, number, coverage_count=0, true_len=0):
        if coverage_count == 0 and true_len == 0:
            _, true = self.get_barcode(number)
            cov = self.get_coverage_count(number)
            x = (self.h5py_object[self.data[0,number-1]][0:true, (1+cov):(1+2*cov)])
            return np.transpose(x)
        else: 
            x = self.h5py_object[self.data[0,number-1]][0:true_len, 1:(coverage_count+1)]
            return np.transpose(x)        
    
    def get_interp_changes(self, number, coverage_count=0, true_len=0): # checked
        if coverage_count == 0 and true_len == 0:
            _, true = self.get_barcode(number)
            x = self.h5py_object[self.data[0,number-1]][0:true,(1+2*self.get_coverage_count(number))]
            return np.transpose(x)
        else:
            x = self.h5py_object[self.data[0,number-1]][0:true_len,(1+2*coverage_count)]
            return np.transpose(x)
    
    def get_interp_mutations(self, number, coverage_count=0, true_len=0):
        if coverage_count == 0 and true_len == 0:
            _,true=self.get_barcode(number)
            x = self.h5py_object[self.data[0,number-1]][0:true,(3+2*self.get_coverage_count(number))]
            return np.transpose(x)
        else:
            x = self.h5py_object[self.data[0,number-1]][0:true_len,(3+2*coverage_count)]     
            return np.transpose(x)      
    
    def get_interp_consensus(self, number, coverage_count=0, true_len=0):
        if coverage_count == 0 and true_len == 0:
            _, true= self.get_barcode(number)
            x = self.h5py_object[self.data[0,number-1]][0:true,(2+2*self.get_coverage_count(number))]
            return np.transpose(x) 
        else:
            x = self.h5py_object[self.data[0,number-1]][0:true_len,(2+2*coverage_count)]
            return np.transpose(x) 
    
def seqblock_parser(seqblock):
    seqblock_parsed = CleanSeqBlock()
    total_columns = len(np.transpose(seqblock)) # is there a more efficient way?
    total_rows = len(seqblock)
    number_of_reads = int((total_rows - 4)/2) # last three rows have intreptations, first row is non-mutated target sequence, and N quality scores for N reads
    
    seqblock_parsed.cov = number_of_reads;
    seqblock_parsed.len = total_columns;
    
    # print(total_columns)
    # print(total_rows)
    tag = ""
    tag_length = 10_000
    
    for i in range(total_rows):
        # print(f" The sequenceblock input {seqblock}") #debug
        # print(f" The m x n size of the seqblock {seqblock.shape}") #debug
        # print(f" The ith row of the seqblock {current_seq}") #debug
        # print(f" The m x n size of the sequence {current_seq.shape}") #debug
        
        # new_format = np.chararray(total_columns) //CHANGED to string
        stop_found = False
        # new_format[:] = 'q' #debug
        
        # print(f" The 0th entry of the initialized char array {new_format[0, 0]}") #debug
        # print(f" Its row size {len(new_format[0])}")
        # if (i == 0):
        #     original = ""
        #     for j in range(total_columns):
        #         original += chr(seqblock[i,j])
                
                
        
        # if (i != 0):
  
        if (i == 0 or i == 1):
            new_format = ""
            j=0
        
            while j < total_columns:
                # print(j)
                                        
                if stop_found is False:
                    if (seqblock[i,j] == ord('T') and (j+6) < total_columns):
                        if (seqblock[i,j+1] == ord('A') and seqblock[i, j+2] == ord('A') and seqblock[i, j+3] == ord('T') and seqblock[i, j+4] == ord('A') and seqblock[i, j+5] == ord('G')):
                            new_format += 'T'
                            new_format += 'A'
                            new_format += 'A'
                            new_format += 'T'
                            new_format += 'G'
                            j = j+6
                            stop_found = True
                    if stop_found is False:
                        new_format += chr(seqblock[i,j])
                        j+=1
            
                elif stop_found is True and i == 1:
                    tag += chr(seqblock[i,j])
                    j+=1
                else:
                    j = total_columns
            if (i == 1): 
                seqblock_parsed.barcode = tag
                tag_length = len(tag)
        else:
            new_format = ""
            for j in range(total_columns-tag_length):
                new_format += chr(seqblock[i,j])
                
                
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
    def __init__(self):
        self.cov, self.len, self.target, self.interp_changes, self.interp_mutations, self.interp_consensus, self.barcode = 0,0,0,0,0,0,0
        self.reads, self.qscores = [], []

def pscore(Qscore): 
    return pow(10, -(Qscore-33)/10)