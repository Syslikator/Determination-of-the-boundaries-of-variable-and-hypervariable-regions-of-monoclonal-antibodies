import sys
import re
import itertools

class FastaReader:
    def __init__(self, fastaFile):
     
        self.fastaFile = fastaFile  
    def readFasta(self):
        seq = ''
        for line in self.fastaFile:
            line = line.strip()
            if line.startswith('>'): 
                if seq != '':
                    yield (head, seq)
                    seq = ''   
                head = line[1:]                
            elif line.isalpha():     
                seq += line     
        yield (head, seq)
        

class CDRsCompiler:
    
    def __init__(self):
        
        self.kappaCDRs = ['W[Y|L|F][Q|L]', #!!!END!!!1st CDR begins with {C}  leng{10-17}  position 24-34
                            '[I|V][Y|F]', #!!!START!!!consists of 7 amino without a signal seq. leng{7}  position 50-56
                            'FG[A-Z]G'] #!!!END!!!3rd CDR begins with {C}   leng{7-11}  position 89-97
        self.lambdaCDRs = ['W[Y|L|F][Q|L]', #!!!END!!!1st CDR begins with {C}  leng{10-17}  position 24-34
                            '[I|V][Y|F]', #!!!START!!!consists of 7 amino without a signal seq. leng{7}  position 50-56
                            'FG[A-Z]G'] #!!!END!!!3rd CDR begins with {C}   leng{7-11}  position 89-97
        self.heavyCDRs = ['CKASG',     #!!!starts with {C} works in murine Vabs/ 'C[A-Z][A-Z][A-Z][A-Z]'  leng{10-12} position 31-35
                            'W[V|I|A]', 
                            'L[Q|E]WIG', # leng{16-19} position 50-65
                            '[K|R][L|I|V|F|T|A][T|S|I|A]', # wide range of aminoacids leng{16-19} position 50-65 
                            'C[A-Z]R', # leng{3-25} position 95-102
                            'WG[A-Z]G'] 

    def compileCDRs(self):
        
        self.compiledKappaCDRs = [[re.compile(CDRs)
                                     for CDRs in self.kappaCDRs],
                                    'kappa']
        self.compiledLambdaCDRs = [[re.compile(CDRs)
                                     for CDRs in self.lambdaCDRs],
                                    'lambda']
        self.compiledHeavyCDRs = [[re.compile(CDRs)
                                     for CDRs in self.heavyCDRs],
                                    'heavy']
        self.allCompiledCDRs = [self.compiledKappaCDRs,
                                  self.compiledLambdaCDRs,
                                  self.compiledHeavyCDRs]

        return self.allCompiledCDRs

class FindAntibodies:
    
    rnaCodonTable = {
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '*', 'UGA': '*',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '*', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
    }
    dnaCodonTable = {key.replace('U','T'): value
                      for key, value in rnaCodonTable.items()}
    
    
    def __init__(self, head, seq, CDRs):

        self.sequence = seq.upper()
        self.header = head
        self.allCompiledCDRs = CDRs
        
        self.fwdrevTranslations = self.translateSeq(self.sequence) + self.getRevComp()
       
                
    def translateSeq(self, seq):        

        fwdAAs = [[], [], []]
        
        joinedAAseqs = []
            

        for frame in range(0, 3):
            for i in range(frame, len(seq), 3):

                if 'N' in seq[i: i+3]:
                    thisAA = 'X'
                    fwdAAs[frame].append(thisAA)
                elif seq[i: i+3] in self.dnaCodonTable.keys():
                    thisAA = self.dnaCodonTable.get(seq[i: i+3])
                    fwdAAs[frame].append(thisAA)
            joinedAAseqs.append(''.join(fwdAAs[frame]))
        
        return(joinedAAseqs)
                          

    def getRevComp(self):

        '''make a dictionary of base pairs'''

        bpDict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}

        revComp = ''.join([bpDict[base] for base in reversed(self.sequence)])
       
        revTranslations = self.translateSeq(revComp)
        return(revTranslations)


    def whereAntibodies(self, head, seq):        

        for translation in self.fwdrevTranslations:
            frame = self.fwdrevTranslations.index(translation)
           
            for CDRs, chaintype in self.allCompiledCDRs:
                CDRsMatches = []
                for CDRs in CDRs:
                    CDRsMatches.append(CDRs.findall(translation))
                flattenedCDRsMatches = itertools.chain(*CDRsMatches)
               
                if all(CDRsMatches):                    
                    print("Complete {} chain found in frame "
                              .format(chaintype), end='')
                    #for forward frames
                    if frame in range(0, 3):
                        print("{} of {}! \nCDRs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(frame+1,
                                        self.header,
                                        ', '.join(flattenedCDRsMatches),
                                        self.header,
                                        frame+1))
                    #for reverse frames
                    else:
                        print("{} of {}! \nCDRs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(-frame+2,
                                        self.header,
                                        ', '.join(flattenedCDRsMatches),
                                        self.header,
                                        -frame+2))
                    print("{}\n".format(translation))

                elif any(CDRsMatches):                                               
                     print("Partial {} chain found in frame "
                              .format(chaintype), end='')
                    #for forward frames
                     if frame in range(0, 3):
                        print("{} of {}. \nCDRs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(frame+1,
                                        self.header,
                                        ', '.join(flattenedCDRsMatches),
                                        self.header,
                                        frame+1))
                    #for reverse frames
                     else:
                        print("{} of {}. \nCDRs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(-frame+2,
                                        self.header,
                                        ', '.join(flattenedCDRsMatches),
                                        self.header,
                                        -frame+2))
                     print("{}\n".format(translation))
                    
CDRs = CDRsCompiler().compileCDRs()
myReader = FastaReader(sys.stdin)
for head, seq in myReader.readFasta():
    checkAllSequences = FindAntibodies(head, seq, CDRs)
    checkAllSequences.whereAntibodies(head, seq)

