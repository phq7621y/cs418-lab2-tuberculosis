import numpy as np
import copy

conv = {}
conv["A"] = 0
conv["C"] = 1
conv["G"] = 2
conv["T"] = 3

genetic_code = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

def create_profile_matrix_with_pseudocounts(motifs):
    profile = np.ones((4,len(motifs[0])))
    
    for i in range(len(motifs[0])):
        for j in range(len(motifs)):
            profile[conv[motifs[j][i]]][i] += 1
        
    profile = profile/(len(motifs)+4)
    return profile
    
def randomly_generated_kmers(k, s, profile):
    the_pattern = ""
    
    random_prop = []
    for i in range(len(s)-k+1):
        pattern = s[i:i+k]
        prop = 1
        for j in range(len(pattern)):
            prop *= profile[conv[pattern[j]]][j]
        random_prop.append(prop)
        
    random_prop = np.array(random_prop)/sum(random_prop)
    pos = np.random.choice(len(s)-k+1,p=random_prop)
    
    the_pattern = s[pos:pos+k]
    
    return the_pattern

def score_motifs(motifs):
    score = 0
    t_s = "ACGT"
    length = len(motifs)
    
    for col in range(len(motifs[0])):
        t_col = ""
        for row in range(len(motifs)):
            t_col += motifs[row][col]
                         
        count = []
        for c in t_s:
            count.append(t_col.count(c))
        score += length-np.max(count)
    
    return score

def gibbs_sampler(k,n,lines):
    motifs = []
    t = len(lines)

    for i in range(0,t):
        col = np.random.randint(0,len(lines[i])-k+1)
        motifs.append(lines[i][col:col+k])
     
    best_motifs =  copy.deepcopy(motifs)
    best_score = score_motifs(best_motifs)
    
    for j in range(n):
        i = np.random.randint(0,t)
        motifs.pop(i)
        profile = create_profile_matrix_with_pseudocounts(motifs)
        motifs.insert(i,randomly_generated_kmers(k,lines[i],profile))
            
        score = score_motifs(motifs)
        if score < best_score:
            best_score = score
            best_motifs = copy.deepcopy(motifs)
    
    return best_motifs

def get_general_motif(motifs):
    general_motif = ""
    t_s = "ACGT"
    length = len(motifs)
    
    for col in range(len(motifs[0])):
        t_col = ""
        for row in range(len(motifs)):
            t_col += motifs[row][col]
                         
        count = []
        for c in t_s:
            count.append(t_col.count(c))
        general_motif += t_s[np.argmax(count)]
    
    return general_motif

def amino_acid_translation(s):
    l = len(s)
    i = 0
    amino_acid = ""
    while i < l:
        rna = s[i:i+3]
        if genetic_code[rna] == "STOP":
            break
        else:
            amino_acid += genetic_code[rna]
        i += 3
    
    return amino_acid


with open("upstream250.txt","r") as f:
    lines = f.read().splitlines()
    t_lines = []
    for i in range(len(lines)):
        if i % 2 != 0:
            t_lines.append(lines[i])
    lines = t_lines
    
    #GibbSampling motifs
    k = 20
    n = 200

    motifs = []
    best_motifs = None
    best_score = float("inf")
        
    for i in range(2000):
        motifs = gibbs_sampler(k,n,lines)
        score = score_motifs(motifs)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    
    general_motif = get_general_motif(best_motifs)

    fo = open("result.txt","w")

    for motif in best_motifs:
        fo.write(motif)
        fo.write("\n")
        print(motif)

    fo.write("\n")
    fo.write(general_motif)

    print("")
    print(general_motif)

    fo.close()
    '''
    #convert dna to rna
    rna = list(general_motif)
    for i in range(len(rna)):
            if rna[i] == "T":
                rna[i] = "U"  
    rna = "".join(rna)

    #convert rna to peptide
    #peptide = amino_acid_translation(rna)

    #Output
    fo = open("result.txt","w")
    for motif in best_motifs:
        fo.write(motif)
        fo.write("\n")
        print(motif)

    fo.write("\n")
    fo.write(general_motif)

    fo.write("\n")
    fo.write(rna)

    fo.write("\n")
    fo.write(peptide)

    print("")
    print(general_motif)
    print(rna)
    print(peptide)

    fo.close()
    '''