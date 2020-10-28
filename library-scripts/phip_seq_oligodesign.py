"""General script for designing oligos for PhIP_Seq."""

import Bio.Seq
import Bio.SeqUtils.CodonUsageIndices
import Bio.Alphabet.IUPAC

import argparse
import os.path
import os
import itertools
import pandas as pd

from Bio.Seq import Seq
from Bio import SeqIO
from operator import itemgetter

cai = Bio.SeqUtils.CodonUsageIndices.SharpEcoliIndex

translation_table = {}
codons = []
amino_acids = []
rt_table = {}

for (_nt1, _nt2, _nt3) in itertools.product(Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters, repeat = 3):
    codon = _nt1 + _nt2 + _nt3
    codons.append(codon)
for codon in codons:
    translation_table[codon] = str(Bio.Seq.Seq(codon).translate())
    if translation_table[codon] not in amino_acids:
        amino_acids.append(translation_table[codon])
        rt_table[translation_table[codon]] = [codon]
    else:
        rt_table[translation_table[codon]].append(codon)

# Dictionary of all the synonymous codons encoding an amino acid and each
# codon's associated cai score. Does not include stop codons.
score_dict = {}
for amino_acid in amino_acids:
    if amino_acid != '*':
        scores = []
        for codon in rt_table[amino_acid]:
            scores.append((codon, cai[codon]))
        score_dict[amino_acid] = scores


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Design overlapping peptides'
                                                  'for PhIP-Seq.')
    parser.add_argument('input_dir', type=str, help='Path to directory'
                        'containing sequences from viruses to make oligos'
                        'for. This directory should contain no other fasta files.')
    parser.add_argument('out_dir', type=str, help="Path to directory for"
                        "output oligo files.")
    parser.add_argument('--length', type=int, default=117, help='Length'
                        'of oligos to order.')
    parser.add_argument('--tile', type=int, default=60, help='Number of'
                        'basepairs to tile by.')
    parser.add_argument('--avoid_motifs', type=list,
                        default=['GAATTC', 'AAGCTT'], help='DNA motifs to'
                        'avoid (i.e. restriction sites for cloning)')
    parser.add_argument('--adaptor5', type=str, default='aggaattctacgctgagt',
                        help="5' adaptor for cloning oligos into phage.")
    parser.add_argument('--adaptor3', type=str, default='tgatagcaagcttgcc',
                        help="3' adaptor for cloning oligos into phage.")
    parser.add_argument('input_type', type=str, choices=['protein', 'dna'],
                        help="Input sequence type, either 'protein' or 'dna'")
    return parser.parse_args()


def translate(oligo):
    """Script for translating DNA more quickly than Bio.Seq
    Args:
        `oligo` (str)
            DNA sequence to translate

    Returns:
        `aa_read` (str)
            Translated DNA sequence

    >>> oligo = 'ATGGGC'
    >>> translate(oligo) == 'MG'
    True
    """
    assert len(oligo) % 3 == 0, "Not translatable due to length"
    aa_read = ''.join([translation_table[oligo[i : i + 3]] for i in range(0, len(oligo), 3)])
    return aa_read


def reverse_translate(protein_seq): 
    """ Reverse_translate viral protein sequence into DNA sequence.

    Codon preferences based on cai for E. coli."""
    revtranslated_seq = []

    for aa in protein_seq: #Go through seq aa by aa
        assert aa in rt_table, 'Ambiguous Sequence. Check {0} for accuracy of protein sequence.'.format(protein_sequence)
        max_codon = max(score_dict[aa], key=lambda item:item[1])[0]
        revtranslated_seq.append(max_codon)

    rt_seq = ''.join(revtranslated_seq)

    return rt_seq


def subsequence(genome, oligo_length, tile):

    """ Subsequence viral genome into oligos of length *oligo_length* tiled by *tile*."""
    subsequences = []
    i = 0

    while i < (len(genome) - oligo_length):
        subsequences.append(genome[i : i + oligo_length])
        i = i + tile

    subsequences.append(genome[len(genome)-oligo_length: len(genome)]) #The final oligo goes from the end of the sequence backwards oligo length.

    return subsequences


def remove_rsites(oligos, codon_scores_by_aa, avoidmotifs):

    clean_oligos = []
    rsites_count = 0
    replace = False

    for motif in avoidmotifs:
        for n in range(len(oligos)):
            if motif in oligos[n]:
                aa_withrs = translate(oligos[n])
                start_length = len(oligos[n])
                replace = True
                while replace:
                    for i in range(len(oligos[n]) - len(motif) + 1):
                        seq = oligos[n][i : i + len(motif)] # sequence starting at i
                        if seq == motif:
                            done = False #need to replace a codon in this seq
                            rsites_count += 1
                            if i % 3 == 0: #The restriction site is in frame.
                                for x in codon_scores_by_aa: #For each amino acid...
                                    for y in codon_scores_by_aa[x]: #...And for each (codon, score) tuple that encodes that aa...
                                        if not done: #If the codon at this location has not already been replaced.
                                            if oligos[n][i:i+3] in y: #If the codon we want to replace is in the score tuple y... 
                                                l = codon_scores_by_aa[x] #...Turn all synonymous codons and their scores into a list
                                                l.sort(key=itemgetter(1),reverse=True)
                                                assert oligos[n][i:i+3] == l[0][0], 'The codon we are replacing is not the highest scoring. We are replacing: {0}'.format(oligos[n][i:i+3])
                                                new_codon = l[1][0] # The new codon is the codon that has the second highset score (so is second in the sorted list)
                                                clean_oligo = oligos[n][:i]+new_codon+oligos[n][i+3:]
                                                oligos[n] = clean_oligo
                                                done = True #The codon has been replaced, don't look at this site anymore. 

                            else: #The restriction site is not in frame. Go through the same steps as if it were. 
                                for x in codon_scores_by_aa:
                                    for y in codon_scores_by_aa[x]:
                                        if not done:
                                            if oligos[n][i-(i%3):i-(i%3)+3] in y: #Replace the first codon that contains part of the restriciton site.
                                                l = codon_scores_by_aa[x]
                                                l.sort(key=itemgetter(1), reverse=True)
                                                assert oligos[n][i-(i%3):i-(i%3)+3] == l[0][0], 'The codon we are replacing is not the highest scoring. We are replacing: {0}'.format(oligos[n][i-(i%3):i-(i%3)+3])
                                                new_codon = l[1][0]
                                                clean_oligo = oligos[n][:i-(1%3)]+new_codon+oligos[n][i-(i%3)+3:]
                                                oligos[n] = clean_oligo
                                                done = True #The codon has been replaced, don't look at other possible synonymous codons.

                    if motif not in oligos[n]: #Make sure motif not in oligo after going through replacements
                        replace = False

                assert len(oligos[n]) == start_length, 'Oligo lengths not maintained while removing restriction sites.'         
                assert aa_withrs == translate(oligos[n]), 'The amino acid sequence has been altered by removing restriction site.'

    print(f"{rsites_count} restriction sites were removed with synonymous substitution.")

    return oligos


def main():
    # parse command line arguments
    args = vars(parse_args())

    indir = args['input_dir']
    outdir = args['out_dir']
    input_type = args['input_type']

    oligo_length = args['length']
    tile = args['tile']
    avoidmotifs = args['avoid_motifs']
    adaptor5 = args['adaptor5']
    adaptor3 = args['adaptor3']


    if not os.path.isdir(indir):
        raise ValueError(f"Input directory {indir} does not exist.")

    os.makedirs(outdir, exist_ok=True)

    all_oligos_df = pd.DataFrame()

    in_file_names = []
    for f in os.listdir(indir):
        if '._' == f[0:2]:
            continue
        else:
            in_file_names.append(f)

    for in_file_name in in_file_names:
        print(f"Designing oligos for {in_file_name}.")
        
        for record in SeqIO.parse(f"{indir}/{in_file_name}", 'fasta'):
            in_record = record

        in_seq = in_record.seq
        if input_type == 'protein':
            dna_seq = reverse_translate(in_seq)
            assert in_seq == translate(dna_seq), "reverse translation did not work"
        elif input_type == 'dna':
            dna_seq = in_seq
        else:
            raise ValueError(f"Invalid input type of: {input_type}.")

        subsequences = subsequence(dna_seq, oligo_length, tile)
        cleaned_oligos = remove_rsites(subsequences, score_dict, avoidmotifs)

        oligos_with_translation = []
        oligo_count = 0
        for oligo in cleaned_oligos:
            translated_oligo = translate(oligo)
            final_oligo = adaptor5 + oligo + adaptor3
            if (oligo_count*tile + oligo_length) <= (len(dna_seq)):
                oligo_start = (oligo_count*tile)/3 + 1
            else:
                oligo_start = (len(dna_seq) - oligo_length)/3 + 1
            oligos_with_translation.append((in_file_name[0:-6], final_oligo,
                                            translated_oligo, oligo_count, oligo_start))
            oligo_count += 1 # Keep track of oligo number separately for each file
        
        print(f"{in_file_name} is {len(dna_seq)//3} amino acids long.")
        print(f"{len(oligos_with_translation)} oligos designed for {in_file_name}.\n")
        
        oligos_df = pd.DataFrame(oligos_with_translation,
                                 columns=['Name', 'Oligo', 'Prot', 'Oligo_Num', 'Prot_Start'])

        all_oligos_df = pd.concat([all_oligos_df, oligos_df], ignore_index=True)

    print(f"{len(all_oligos_df)} oligos designed.")
    all_oligos_df.to_csv(f"{outdir}/oligos_all.txt", index_label='Index', 
            columns=['Name', 'Oligo', 'Prot', 'Prot_Start'])
    
    
    all_dup_info = pd.DataFrame(columns = ['Dup_Seq', 'Num', 'Seq_Names'])
    all_oligos_df_dups = all_oligos_df[all_oligos_df.duplicated(['Prot'], keep=False)].reset_index()
    all_oligos_df_dups.to_csv(f"{outdir}/oligos_dups.txt", index_label='Index', 
            columns=['Name', 'Oligo', 'Prot', 'Prot_Start'])
    dup_seqs = all_oligos_df_dups.groupby('Prot')
    for seq, info in dup_seqs:
        dup_info = {'Dup_Seq': seq, 'Dup_Count': len(info),
                    'Seq_Names': [info['Name'].to_list()]}
        all_dup_info = pd.concat([all_dup_info, pd.DataFrame.from_dict(dup_info)], sort=False)
    all_dup_info.to_csv(f"{outdir}/oligos_dup_counts.txt", sep='\t', index_label='Index',
                        columns=['Dup_Seq', 'Dup_Count', 'Seq_Names'])

    all_oligos_df_nodups = all_oligos_df.drop_duplicates(['Prot'], keep='first')
    print(f"{len(all_oligos_df_nodups)} oligos remain after dropping duplicate protein sequences.")
    all_oligos_df_nodups.to_csv(f"{outdir}/oligos_nodups.txt", index_label='Index', 
            columns=['Name', 'Oligo', 'Prot', 'Prot_Start'])



if __name__ == '__main__':
    main()
