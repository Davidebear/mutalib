# Binary Classification Problem: True/False Mutation Detection
_Summer 2022: NeuroToolBox project under Yiyang Gong, PhD (Gong Lab)_ <br>
<br>
Current Objective:  Understand the dataset. Figure out how to extract desired pieces 

## Drivetrain Approach
- _Objective_:  Determine whether a SNP (single nucleotide polymorphism) is truly present or a consequence of sequencing error. Create a dicitonary of mutations associated with a given <br>
- _Levers (Inputs)_:  _Phred_ quality score, probability distribution of random sequencing error from BFP (blue fluorescent protein), error patterns <br>
- _Data needed_:  Randomly parse **high coverage** alignments into training, validation ("parameter tuning"), and test ("restricted") sets <br>
- _Model_:  First, binary classification, then maybe predictive sequence => **Read publications for insight** <br>

## Dataset (.mat file)
whos align3 -> 1 x 355104 <br>
Cell array (align3) of size 1 x 355104 with each column housing a 1x1 cell arary which itself houses a m x n **char array** (more details...) <br>
|- 1 x 355104 cell array (align3) <br>
|--- 1 x 1 cell array  <br>
|----- m x n char array  <br>

_Example_ <br>
align3{1009}(:,end-250:end-100) <br>

### Details about each of the 355104 char arrays
**Important**: **Barcode is all nucleotides after stop codon = 'TAATAG'**

- Row 1: target sequence (non-mutated BFP or RFP DNA sequence)
- Row 2 -> N+1: N reads; N coverage
- Row N+2 -> 2N+1: QV, quality scores, for each N read

Last 3 Rows (index backwards):
- 1st: Subjective corrections (by Gong) 'y = indel, z = base-error'
- 2nd: Subject consensus sequence 
- 3rd (-1): Difference between subj. cons. (2nd last) and target seq. (row 1) 'x = mutation'
