'''
Ad hoc Python library made to handle ASSEMBLY MATRICES specific to Yiyang Gong's lab's MATLAB dataset
@ Davidebear
'''
# --------
# Imports
# --------
from math import pow # Convenience   
from Bio.Seq import Seq # Use Bio.Seq to translate() DNA sequences

import numpy as np 
import pandas as pd

import h5py # Only to handle .mat dataset

STOP_CODON = np.array([84, 65, 65, 84, 65, 71]) # T A A T A G
LINKER = Seq('KLNPPDESGEF') # Used within get_br()
LINKER_CODLEN = 11 # Used within get_br()

LINKER_NUCLEN = 33 # 3 * LINKER_CODLEN
RFP_NUCLEN = 828 # CALCULATED from len(NUC_RFP)
BFP_NUCLEN = 708


# Handling the .mat dataset
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

# For quick conversions from ndarray object to_...
def to_string(arr):
    string_seq = ""
    for i in arr:
        string_seq += chr(i)
    return string_seq
    
def to_seq(targ, string=False): 
    if string: # If it's already a string
        return Seq(targ)
    return Seq(to_string(targ))    

# Used in two DNA_SeqBlocks() funcs: get_br(), br_splitter()
def codon_to_nuc_idx(codon_idx, search_size): return (codon_idx+search_size)*3

def parsed_targ(targ):
    """Convert an np array of ord() DNA nucleotides with indels (45) to a string without indels. Also returns where the indels were.
    Args:
        targ (ndarray): get_target() output
    Returns:
        cleaned_targ, indices_of_indels: removed indels, site of removal
    """
    indxs = (np.where(targ==45)[0]) # Dunno why they made it like this
    return Seq(to_string(targ).replace("-", "")), indxs

def rename(df, i, red=False, blue=False, transpose=False):
    """Syntatic sugar for BFP and RFP pd.DataFrame/s

    Args:
        df (pd.DataFrame): RFP or BFP-parsed seqblock in a df
        i (int): seqblock number
        red (bool, optional): Name as red. Defaults to False.
        blue (bool, optional): Name as blue.. Defaults to False.
        transpose (bool, optional): If transpose wasn't done previously. Defaults to False.

    Returns:
        df: Renamed df. Improved redability
    """
    if transpose == True: df = df.T
    cov = int(0.5*(df.columns.size-4))
    based_on_cov = ['target']
    based_on_cov.extend([f'read{i}' for i in range(1, cov+1)])
    based_on_cov.extend(x for x in [f'q{i}' for i in range(1, cov+1)])
    based_on_cov.extend(['changes', 'contig', 'mutations'])
    df.columns = based_on_cov
    
    to_name = ""
    if cov < 10: to_name += f'0{cov}'
    else: to_name += str(cov)
    
    if blue: to_name += 'B'
    if red: to_name += 'R'
    
    if i < 10: to_name += f'-----{i}'
    elif i < 100: to_name += f'----{i}'
    elif i < 1000: to_name += f'---{i}'
    elif i < 10_000: to_name += f'--{i}'
    elif i < 100_000: to_name += f'-{i}'
    else: to_name += f'{i}'
    df.index.name = to_name
    
    return df

# To house main dataset
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
        
    def get_seqblock(self, number, df=False): # Not really any reason to use this any longer 
        """Grabs entire seqblock
        """
        x = self.h5py_object[self.data[0, number-1]][:, :] # //UNIQUE: 0th column because .mat file, align3 variable only has one row
        if df:
            return x
        return np.transpose(x)
    
    def get_redsb(self, number, red_start, red_end=-1, df=False): # Used within get_br()
        """RFP start known. Get only the RFP

        Args:
            number (int): seqblock #
            red_start (int): index to begin slice
            red_end (int, optional): index to end slice. Defaults to -1.
            df (bool, optional): Transpose not needed if going into a df. Defaults to False.

        Returns:
            h5py_object
        """
        if red_end == -1:
            x = self.h5py_object[self.data[0, number-1]][red_start:, :]
        else:
            x = self.h5py_object[self.data[0, number-1]][red_start:red_end, :]
        if df:
            return x
        return np.transpose(x)
            
    def get_bluesb(self, number, blue_end, df=False): # Used within get_br()
        x = self.h5py_object[self.data[0, number-1]][:blue_end+1, :]
        if df:
            return x
        return np.tranpose(x)
    
    def get_target(self, number, end=-1): # Used within get_br()
        if end != -1:
            x = self.h5py_object[self.data[0, number-1]][:end, 0]
            return np.transpose(x)
        x = self.h5py_object[self.data[0, number-1]][:, 0]
        return np.transpose(x)
    
    def get_len(self, number):
        x = self.h5py_object[self.data[0,number-1]][:,0]
        return len(np.transpose(x))
    
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
    
    def get_barcode_old(self, number): # Crucial. Garbage code but effective.
        """Barcode finder
        
        Changes that need to be made: use target sequence to find stop codon (reads may have mutated stop codon like 22776). Use read that doesn't terminate in indels to find the barcode.

        Args:
            number (int): which sequence? matlab indexing

        Returns:
            barcode, true_len: returns barcode uint16 ndarray, last matlab idx of stop codon
        """
        x = self.h5py_object[self.data[0, number-1]][:, 0] #assess with the targ
        r = self.h5py_object[self.data[0, number-1]][:, 1]
        barcode = 0 #tbr
        true_len = 0 #tbr
        total = x.size
        
        idx = total-1 #set to the end
        search_stop=True
    
        while idx >= 0 and search_stop: # CHANGED: This finds the stop codon if it has no indels
            
            if (x[idx] == 71): #Has to start with a G
                # Best Case Scenario
                if (x[idx-1] == 65):
                    if (x[idx-2] == 84):
                        if (x[idx-3] == 65):
                            if (x[idx-4] == 65):
                                if (x[idx-5] == 84):
                                    true_len = idx+1
                                    search_stop = False
                                    
                #                 elif x[idx-5] == 45: # What if T-AATAG?
                #                     if x[idx-6] == 84:
                #                         true_len = idx+1
                #                         search_stop = False
                                
                #             elif x[idx-4] == 45: # What if TA-ATAG?
                #                 if x[idx-5] == 65:
                #                     if x[idx-6] == 84:
                #                         true_len = idx+1
                #                         search_stop=False
                        
                #         elif x[idx-3] == 45: # What if TAA-TAG?
                #             if x[idx-4] == 65:
                #                 if x[idx-5] == 65:
                #                     if x[idx-6] == 84:
                #                         true_len = idx+1
                #                         search_stop = False
                    
                #     elif x[idx-2] == 45: # What if TAAT-AG?
                #         if x[idx-3] == 84:
                #             if x[idx-4] == 65:
                #                 if x[idx-5] == 65:
                #                     if x[idx-6] == 84:
                #                         true_len = idx+1
                #                         search_stop = False
                
                # elif (x[idx-1] == 45): # What if TAATA-G?
                #     if (x[idx-2] == 65):
                #         if x[idx-3] == 84:
                #             if x[idx-4] == 65:
                #                 if x[idx-5] == 65:
                #                     if x[idx-6] == 84:
                #                         true_len = idx+1
                #                         search_stop = False
    
            idx -= 1
        if true_len != 0:
            
            barcode = r[true_len:true_len+15] # Barcode is only the first 15 nucleotides anyway apparently
            
    
            if (barcode[-1] == 45 and barcode[-2] == 45): # only enters if the current barcode has two indels at the end
                list = [barcode]
                cov = self.get_coverage_count(number)
                if cov == 1: # no other option
                    if np.sum(barcode == 45) == 15:
                        return r[true_len+15:true_len+30], true_len # Only one case where this occurred (SB 236_645)
                    return barcode, true_len
            
                read_num = 2
                while (read_num <= cov):
                    barcode = self.h5py_object[self.data[0, number-1]][true_len:true_len+15, read_num]
                    if barcode[-1] == 45 and barcode[-2] == 45:
                        list.append(barcode)
                        read_num +=1
                    else: 
                        return barcode, true_len # if a better one is found
                
                # if they all have indels at the end
                
                bestindex = 0
                min_indel = 15
                for i in range(len(list)):
                    indel_freq = np.sum(list[i] == 45)
                    if (indel_freq) < min_indel:
                        min_indel = indel_freq
                        bestindex = i
                
                if min_indel == 15: # If the first 15 nucleotides are all indels. Special case scenario found at SB235_194
                    return r[true_len+15:true_len+30], true_len # Not really any reason to find the best. It's good enough to do this.
                return list[bestindex], true_len
                    
                
            return barcode, true_len
        return [],0
    
    def classify_stop(self, number):
        """Classify stop codon type or presence

        Args:
            number (_type_): _description_
        Returns:
            int: [0 == missing, 1 == clean, 2 == mutated]
        """
        targ = self.get_target(number)
        clean_targ, _ = parsed_targ(targ)
        
        if targ.find("TAATAG") != -1:
            return 1
        # This means it wasn't found
        if clean_targ.find("TAATAG") != -1:
            return 0 
        # This means there indeed is a indel_filled one
        return 2 # Further processing TBD to determine if its a true or fake mutation. //TODO
    
    def missing_stop(self, number): # Returns True if at least one read does not have a stop codon
        """Check for missing stop codon

        Args:
            number (no.): which sb?

        Returns:
            bool: True if stop codon is absent. Ignores indels.
        """
        clean_targ, _ = parsed_targ(self.get_target(number))
        
        if clean_targ.find("TAATAG") == -1: 
            return True
        return False # Can be fed into indel_stop to check for indel present in target

    def classify_indel_stop(self, number, run_missing=False): #// TODO: discriminate between a true mutation and a fake mutated stop
        """Check if stop codon has a potential mutation according to the reads. Not confirmed if a true or fake mutation

        Args:
            number (int): which sb?
            
        Returns:
            bool: True if mutation present
        """
        targ = self.get_target(number)
        indel_stop = True # only changes if the stop hasn't been confirmed to be clean or present yet
        
        if run_missing:
            if ~self.missing_stop(number): # if this returns False... meaning there is indeed a stop
                if targ.find("TAATAG") != -1: # Means a clean stop codon is not present...
                    # Confirms that there is an indel....
                    indel_stop = False
        if indel_stop:
            # Discriminate.  Now how do we decide if its a true indel stop or a fake one?
            pass
    
    def get_barcode(self, number): # Will find the stop codon mapping no matter what.
        targ, indels = parsed_targ(self.get_target(number)) # Removes all indels 
        idx = targ.find("TAATAG")
        
        if idx == -1:
            return [], 0
    
        idx += 6 # Gives true_len
        
        if len(indels) != 0:
            for n in indels: # To find the true index... You have to add +1 to the found nuc_idx (final_idx) if an indel was previously in the range...
                if n < idx:
                    idx +=1
        
        
        
        r = self.h5py_object[self.data[0, number-1]][:, 1]
        
        true_len = idx
        
        barcode = r[true_len:true_len+15] # Barcode is only the first 15 nucleotides anyway apparently
        
        if (barcode[-1] != 45 and barcode[-2] != 45):
            return barcode, true_len
        
        # only enters if the current barcode has two indels at the end
        list = [barcode]
        cov = self.get_coverage_count(number)
        if cov == 1: # no other option
            if np.sum(barcode == 45) == 15:
                return r[true_len+15:true_len+30], true_len # Only one case where this occurred (SB 236_645)
            return barcode, true_len
    
        read_num = 2
        while (read_num <= cov):
            barcode = self.h5py_object[self.data[0, number-1]][true_len:true_len+15, read_num]
            if barcode[-1] == 45 and barcode[-2] == 45:
                list.append(barcode)
                read_num +=1
            else: 
                return barcode, true_len # if a better one is found
        
        # if they all have indels at the end
        
        bestindex = 0
        min_indel = 15
        for i in range(len(list)):
            indel_freq = np.sum(list[i] == 45)
            if (indel_freq) < min_indel:
                min_indel = indel_freq
                bestindex = i
        
        if min_indel == 15: # If the first 15 nucleotides are all indels. Special case scenario found at SB235_194
            return r[true_len+15:true_len+30], true_len # Not really any reason to find the best. It's good enough to do this.
        return list[bestindex], true_len
            
    def get_nuc_qscores(self, number, nuc_pos, coverage_count=0): # Before pd.DF use
        if coverage_count == 0:
            cov = self.get_coverage_count(number)
            x = self.h5py_object[self.data[0, number-1]][nuc_pos-1, (1+cov):(1+2*cov)]
            return np.transpose(x)
        else:
            x = self.h5py_object[self.data[0, number-1]][nuc_pos-1, (1+coverage_count):(1+2*coverage_count)]
            return np.transpose(x)
    
    def get_nuc_reads(self, number, nuc_pos, coverage_count=0): # Before pd.DF use
        if coverage_count == 0:
            x = (self.h5py_object[self.data[0,number-1]][nuc_pos-1, 1:self.get_coverage_count(number)+1])
            return np.transpose(x)
        else: 
            x = self.h5py_object[self.data[0,number-1]][nuc_pos-1, 1:(coverage_count+1)]
            return np.transpose(x)
    
    def get_reads(self, number, coverage_count=0, true_len=0): # Before pd.DF use
        if coverage_count == 0 and true_len == 0:
            _, true = self.get_barcode(number)
            x = (self.h5py_object[self.data[0,number-1]][0:true, 1:self.get_coverage_count(number)+1])
            return np.transpose(x)
        else: 
            x = self.h5py_object[self.data[0,number-1]][0:true_len, 1:(coverage_count+1)]
            return np.transpose(x)
    
    def get_qscores(self, number, coverage_count=0, true_len=0): # Before pd.DF use
        if coverage_count == 0 and true_len == 0:
            _, true = self.get_barcode(number)
            cov = self.get_coverage_count(number)
            x = (self.h5py_object[self.data[0,number-1]][0:true, (1+cov):(1+2*cov)])
            return np.transpose(x)
        else: 
            x = self.h5py_object[self.data[0,number-1]][0:true_len, 1:(coverage_count+1)]
            return np.transpose(x)        
    
    def get_interp_changes(self, number, coverage_count=0, true_len=0):# Before pd.DF use
        if coverage_count == 0 and true_len == 0:
            _, true = self.get_barcode(number)
            x = self.h5py_object[self.data[0,number-1]][0:true,(1+2*self.get_coverage_count(number))]
            return np.transpose(x)
        else:
            x = self.h5py_object[self.data[0,number-1]][0:true_len,(1+2*coverage_count)]
            return np.transpose(x)
    
    def get_interp_mutations(self, number, coverage_count=0, true_len=0): # Before pd.DF use
        if coverage_count == 0 and true_len == 0:
            _,true=self.get_barcode(number)
            x = self.h5py_object[self.data[0,number-1]][0:true,(3+2*self.get_coverage_count(number))]
            return np.transpose(x)
        else:
            x = self.h5py_object[self.data[0,number-1]][0:true_len,(3+2*coverage_count)]     
            return np.transpose(x)      
    
    def get_interp_consensus(self, number, coverage_count=0, true_len=0): # Before pd.DF use
        if coverage_count == 0 and true_len == 0:
            _, true= self.get_barcode(number)
            x = self.h5py_object[self.data[0,number-1]][0:true,(2+2*self.get_coverage_count(number))]
            return np.transpose(x) 
        else:
            x = self.h5py_object[self.data[0,number-1]][0:true_len,(2+2*coverage_count)]
            return np.transpose(x) 

    def br_splitter(self, i, end=0): # Expects TRUELENS to be loaded
        """Get nucIDX where BFP terminates and nucIDX where RFP begins. Uses a match with protein sequence of linker 

        Args:
            i (int): sb no. (MATLAB) 
            Internal: Expects TRUELENS to be loaded

        Returns:
            2 indices: nuc_idx where bfp ends, nuc_idx where rfp begins
        """
        if end == 0: # Stop is missing 
            raw_targ = (self.get_target(i)) # give it less to search by specifying an end
        else: 
            raw_targ = self.get_target(i, end=end)
        # raw_targ = self.get_target(i)
        
        targ, indel_coords = parsed_targ(raw_targ) # Removes indels, converts target into a Bio.Seq object
        
            
        p_targ = targ.translate()
        idx = p_targ.find(LINKER) # Looks for the linker... but what if it can't find one?
        
        # No FULL linker found
        if idx == -1: 
            size = len(targ)
            # Let's see if its obviously red...
            if end != 0 and size <= (RFP_NUCLEN + LINKER_NUCLEN): # If it has a stop codon... and its within the nuc_length of the RFP + linker it's RE
                # Red confirmed but let's remove any underlying linkage
                # Doesn't check for indels though ...
                ridx = targ.find("ATGGTG") # The first 6 aas of RFP
                if ridx == -1: # No beginning of red/linker but has STOP
                    return None, 0
                if len(indel_coords) != 0:
                    for n in indel_coords: # To find the true index... You have to add +1 to the found nuc_idx (final_idx) if an indel was previously in the range...
                        if n > ridx:
                            break
                        else:
                            ridx +=1
                return None, ridx
            
            # No stop codon
            
            
            # KEY CONCEPT: Since no linker is present. There can only be BFP or RFP not both
            # First see if we can find the end of blue
            nuc_bfp_end = targ.find("AGATCT")

            # Can find blue_li, lue, ue, etc. But what if end is missing?
            if nuc_bfp_end != -1:
                nuc_bfp_end += 5 # Python index returned I believe...
                return nuc_bfp_end, None
            
            check_blue = targ.find("ATCAGAGG") # Arbitrarily grabbed from middle of targ
            
            if check_blue != -1: # BLUE PARTIAL using middle of NUC_BFP
                return 0, None
            
            check_red = targ.find("AAGAAGA")
            
            if check_red != -1: # RED PARTIAL using middle of NUC_RFP
                
                return None, 0
            
            
            print(f"{i} is either a red partial or a blue partial that couldn't be foudn by nucs in the middle" )
            return None, None # For now
                
        # Rare: Full linker at start of sequence. No blue.
        if idx == 0: 
            nuc_bfp_end = None # DNE
            nuc_rfp_beg = codon_to_nuc_idx( 0, search_size=LINKER_CODLEN)
            if len(indel_coords) != 0:
                for n in indel_coords:
                    if n <= nuc_rfp_beg:
                        nuc_rfp_beg += 1
            return None, nuc_rfp_beg
        
        # BEST CASE SCENARIO... Full linker Present
        codon_bfp_end = idx
        codon_rfp_beg = idx
        
        # Adjust based on indels
        nuc_bfp_end = codon_bfp_end*3-1
        nuc_rfp_beg = codon_to_nuc_idx( codon_rfp_beg, search_size=LINKER_CODLEN) # accounts for the + linker_codlen
    
        # for BFP -> shift RFP as well
        if len(indel_coords) != 0:
            for n in indel_coords: # To find the true index... You have to add +1 to the found nuc_idx (final_idx) if an indel was previously in the range...
                if n <= nuc_bfp_end:
                    nuc_bfp_end+=1
                    nuc_rfp_beg+=1
                elif n <= nuc_rfp_beg: # only gets here if its bigger than nuc_bfp_end   
                    nuc_rfp_beg+=1
        
        return nuc_bfp_end, nuc_rfp_beg

    def get_br(self, i, end): # Expects TRUELENS to be loaded
        """Gives the pd.DF objects that house the RFP and BFP for the specified ith seqblock

        Args:
            i (int): seqblock number

        Returns:
            BFP, RFP (pd.DataFrame): the parsed target with indices aligned to the target within the seqblock
        """
        b, r = self.br_splitter(i, end)
        red_end = end
        
        if b is not None:
            blue = pd.DataFrame(self.get_bluesb(i, blue_end=b, df=True))
            blue = rename(blue, i, blue=True)
        else:
            blue = pd.DataFrame(None)
        
        if r is not None:
            if red_end != 0:
                red = pd.DataFrame(self.get_redsb(i, red_start=r, red_end=red_end, df=True), index=np.arange(r, red_end))
            else:
                insert = self.get_redsb(i, red_start=r, df=True)
                red = pd.DataFrame(insert, index=np.arange(r, len(insert)+r))
            red = rename(red, i, red=True)
        else:
            red = pd.DataFrame(None)
        
        return blue, red

# For conversions
def pscore(Qscore): 
    return pow(10, -(Qscore-33)/10)

# Deprecated functions
def seqblock_parser(seqblock): # Deprecated
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
class CleanSeqBlock(): # Deprecated
    def __init__(self):
        self.cov, self.len, self.target, self.interp_changes, self.interp_mutations, self.interp_consensus, self.barcode = 0,0,0,0,0,0,0
        self.reads, self.qscores = [], []
