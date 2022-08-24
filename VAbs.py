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
        

class MotifCompiler:
    
    def __init__(self):
        
        self.kappaMotifs = ['[I|T|V][A-Z][L|M][T|S]Q[S|T|P][P|H|T|S]',
                            'P[S|D|A|V][R|H]F[S|T|R]GS[G|D|N|R]',
                            'TFG[G|Q|A|S|T]?GT']
        self.lambdaMotifs = ['[V|L][T|H]Q[E|S|P][S|P|A][A|S|L][L|A|V][T|S][T|S|F|G]',
                             '[P|D|S][G|D][V|I|L][P|S][A|V|D|P|N]RFSGS[L|K|S]',
                             '[F|W][G|E][G|S|T][G|E|R][T|K][K|R|S][L|V|Y|T][T|L][V|D|W]']
        self.heavyMotifs = ['V[Q|K|N|M]L[V|K|Q|L|H][Q|E][S|P]G',
                            'W[V|I][R|K][Q|K][A-Z][P|N|H]',
                            '[Q|E|A|H]G[T|S][L|T|S|M][V|L][T|A]VS[S|A]']

    def compileMotifs(self):
        
        self.compiledKappaMotifs = [[re.compile(motif)
                                     for motif in self.kappaMotifs],
                                    'kappa']
        self.compiledLambdaMotifs = [[re.compile(motif)
                                     for motif in self.lambdaMotifs],
                                    'lambda']
        self.compiledHeavyMotifs = [[re.compile(motif)
                                     for motif in self.heavyMotifs],
                                    'heavy']
        self.allCompiledMotifs = [self.compiledKappaMotifs,
                                  self.compiledLambdaMotifs,
                                  self.compiledHeavyMotifs]

        return self.allCompiledMotifs

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
    
    
    def __init__(self, head, seq, motifs):


        self.sequence = seq.upper()
        self.header = head
        self.allCompiledMotifs = motifs
                

        self.fwdrevTranslations = self.translateSeq(self.sequence) + self.getRevComp()
       
                
    def translateSeq(self, seq):
        '''Translate DNA sequences assuming a forward 5' to 3' sequence.'''
        

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
           
            for motifs, chaintype in self.allCompiledMotifs:
                motifMatches = []
                for motif in motifs:
                    motifMatches.append(motif.findall(translation))
                flattenedMotifMatches = itertools.chain(*motifMatches)
               
                if all(motifMatches):                    
                    print("Complete {} chain found in frame "
                              .format(chaintype), end='')

                    if frame in range(0, 3):
                        print("{} of {}! \nMotifs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(frame+1,
                                        self.header,
                                        ', '.join(flattenedMotifMatches),
                                        self.header,
                                        frame+1))
                    else:
                        print("{} of {}! \nMotifs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(-frame+2,
                                        self.header,
                                        ', '.join(flattenedMotifMatches),
                                        self.header,
                                        -frame+2))
                    print("{}\n".format(translation))

                elif any(motifMatches):                                               
                     print("Partial {} chain found in frame "
                              .format(chaintype), end='')

                     if frame in range(0, 3):
                        print("{} of {}. \nMotifs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(frame+1,
                                        self.header,
                                        ', '.join(flattenedMotifMatches),
                                        self.header,
                                        frame+1))

                     else:
                        print("{} of {}. \nMotifs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(-frame+2,
                                        self.header,
                                        ', '.join(flattenedMotifMatches),
                                        self.header,
                                        -frame+2))
                     print("{}\n".format(translation))
                    
motifs = MotifCompiler().compileMotifs()
myReader = FastaReader(sys.stdin)
for head, seq in myReader.readFasta():
    checkAllSequences = FindAntibodies(head, seq, motifs)
    checkAllSequences.whereAntibodies(head, seq)
