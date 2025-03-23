from collections import Counter
import numpy as np
DNA_Nucleotide=['A','C','G','T']
DNA_ReverseComplement={'A':'T','T':'A','C':'G','G':'C'}

DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

RNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUA": "I", "AUU": "I", "AUC": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "UAA": "_", "UAG": "_", "UGA": "_"
}

tobacco_biased={
'F':'UUU',
'L':'CUU',
'I':'AUU',
'M':'AUG',
'V':'GUU',
'S':'UCU',
'P':'CCU',
'T':'ACU',
'A':'GCU',
'Y':'UAU',
'H':'CAU',
'Q':'CAA',
'N':'AAU',
'K':'AAG',
'D':'GAU',
'E':'GAA',
'C':'UGU',
'_':'UGA', #stop codon
'W':'UGG',
'R':'AGA',
'G':'GGU'}

def fasta_to_list(fasta_dir,seq_to_codon=False,separate_aa=False,sos_eos=True): 
    """convert FASTA format to a list format of the sequences.

    Parameters
    ----------
        fasta_dir: your fasta file directory.
        seq_to_codon: if your sequence is a DNA or RNA sequence and you wantto convert it
                      to codon format (triple form), you can set it to True.
        codons_separator: the separator to separate codons.
        separate_aa: if your sequence is an amino acid sequence and for your algorithm you
                     need to separate amino acids from each other, you can set it to True.
        sos_eos: if True, start of sequence character (<SOS>) and end ofsequence character
                 (<EOS>) will be added to the start and the end of a sequence.
                     
"""
    codons_separator=" "
    if fasta_dir!=None:
        with open(fasta_dir) as f:
            l=f.readline()
            seq=[]
            sub_seq=''
            while len(l)!=0:
                if '>' in l:
                    l=f.readline()
                    if len(sub_seq)!=0:
                        seq.append(sub_seq)
                        sub_seq=''
                else:
                    sub_seq+=l.strip()
                    l=f.readline()
                    if len(l)==0:
                        seq.append(sub_seq)
        if seq_to_codon:
            seq=seq_to_cds(seq,codons_separator,sos_eos)
        if separate_aa:
            seq=separate_amino_acids(seq)
    return seq

def seq_to_cds(seqs,sep=None,sos_eos=True):
    """seperate nucleotide sequence to codons (triple form).

    Parameters
    ----------
        seqs: the list of sequences for converting to codon format.
        sep:  the separator you choose to separate codons based on it.
              if you don't specify any separator(default) each codons
              will separate by whitespace.
        sos_eos: if True, start of sequence character SOS and end of
                 sequence character EOS will be added to the start 
                 and the end of a sequence.

Example:
import structure as s
cds=s.fasta_to_list("C:\\Users\\farsh\\Downloads\\cds.txt")
codons=s.seq_to_cds(cds)
    """
    seqs_list=[]
    for seq in seqs:
        str_cds=[]
        for pos in range(0,len(seq)-len(seq)%3,3):
            str_cds.append(seq[pos:pos+3])
        if str_cds[-1] in ("UAA","UAG","UGA","TAA","TAG","TGA"): #remove stop codons
            del str_cds[-1]
        if sos_eos:
            str_cds.insert(0,'[SOS]')
            str_cds.append('[EOS]')
        seqs_list.append(sep.join(str_cds))
    return seqs_list

def separate_amino_acids(aa):
    aa_list=[]
    for amino in aa:
        aa_list.append(" ".join(list(amino)))
    return aa_list

def translate_seq_to_aa(seq):
    if 'U' in seq:
        cds_ref=RNA_Codons
    else:
        cds_ref=DNA_Codons
    seq_list=seq_to_cds(seq)
    aa_list=""
    for i in seq_list:
        aa_list+=cds_ref[i]
    return aa_list

def codon_biased(protein_seq=None,fasta_dir=None,codon_sep=True,sep_char=','):
    """With this function, the amino acid sequence (protein) is translated to the mRNA sequence based on codon bias in Nicotiana tabacum.

    Parameters
    ----------
        protein_seq: it can be a string or FASTA format of the proteins.
        fasta_dir: the dir address of your FASTA file.
        codon_sep: if True, your codon seq will be seperated from each other (eq. triplet sequences)
        sep_char: the sep character to join codons together, default: ','. it won't work if the codon_sep=False.
        
    """
    if fasta_dir!=None:
        protein_seq=fasta_to_list(fasta_dir) 
    codons=[]         
    for aa in protein_seq:
        seq=[]
        aa = aa.replace('\n', '')
        for p in aa:
            seq.append(tobacco_biased[p])
        if codon_sep:
            codons.append(seq)
        else:
            codons.append(','.join(seq).replace(',',sep_char))

    return codons

def cai_calc(high_expressions=None,fasta_dir=None):
    if fasta_dir!=None:
        high_expressions=fasta_to_list(fasta_dir)
    seq_list=[]
    for s in high_expressions:
        seq_list.append(seq_to_cds(s,sep=" "))
        
    keys = tobacco_biased.keys()
    aa_dict = {key: {} for key in keys}
    if 'T' not in high_expressions[0]:
        cds_ref=RNA_Codons
    else:
        cds_ref=DNA_Codons
    nucleotide_freq=Counter()
    for cds in seq_list:
        nucleotide_freq.update(cds)
    nucleotide_freq=dict(nucleotide_freq)
    for k,v in nucleotide_freq.items():
        amino_acid=cds_ref[k]
        aa_dict[amino_acid][k]=(v)
    for amino_acid, codons in aa_dict.items():
        max_value = max(codons.values())    
        for codon in codons:
            codons[codon] = (codons[codon],np.round(codons[codon] / max_value,3))
    return aa_dict
    
    
    
    


    
