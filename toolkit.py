from collections import defaultdict
from datetime import datetime
from biopytk import sequenceBuilder as sb
from biopytk import fasta_tk as ftk

def toFASTA(ancestry_file, outfile_name = 'output_'+ datetime.now().strftime("%Y%m%d_%H%M%S")+'.fasta', al = 0):
    """
    Parse a standard Ancestry DNA .txt file into FASTA format
    \nNotes: Will parse through file and disregard until the "rsid	chromosome	position	allele1	allele2" line of text is found.
    \n\tBy default, al = 0 will return a fasta file with allele 1 appended by allele 2 in singular file, 1 for allele 1, 2 for allele 2
    \n<- (AncestryDNA.txt): default Ancestry DNA .txt file from download , outfile_name: str, al: int
    \n-> outfile_name.fasta
    """    
    assert al == 0 or al == 1 or al == 2 
    lines = sb.readFile(ancestry_file)
    for line in lines[:]:
        if line[0] == '#':
            lines.remove(line) # remove beginning comments
        if line[0] != '#':
            print(line)
            lines.remove(line) # remove the header line and break
            break
    allele1 = defaultdict(list)
    allele2 = defaultdict(list)
    
    count = 0
    for line in lines:
        count += 1
        stripped = line.strip()
        line_arr = stripped.split('\t')
        chromosome = line_arr[1]
        if chromosome not in allele1.keys():
            count = 0 # reset wrapping for new chr/label 
        allele1[chromosome].append(line_arr[3].replace('0', '-')) # Ancestry deletions are represented by 0
        allele2[chromosome].append(line_arr[4].replace('0', '-'))
        if count >= 60:
            allele1[chromosome].append('\n')
            allele2[chromosome].append('\n')
            count = 0

        

    allele1_out = ""
    for chr, seq in allele1.items():
        allele1_out += ">allele_1_chr_" + chr + '\n'
        allele1_out += ''.join(seq) + '\n'

    allele2_out = ""
    for chr, seq in allele2.items():
        allele2_out += ">allele_2_chr_" + chr + '\n'
        allele2_out += ''.join(seq) + '\n'
      
    
    if al == 0:
        sb.writeFile([allele1_out, allele2_out], outfile_name)
    elif al == 1:
        sb.writeFile(allele1_out, outfile_name)
    elif al == 2:    
        sb.writeFile(allele2_out, outfile_name)

