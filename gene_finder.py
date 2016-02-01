# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Duncan Hall

"""

import random
# from amino_acids import aa, codons, aa_table    # you may find these useful
# from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


# Some Resources...
nucleotide_complements = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
}

stop_codons = ['TAA', 'TAG', 'TGA']
start_codon = 'ATG'


def get_codons(dna):
    codons = []
    for i in range(0, len(dna), 3):
        codons.append(dna[i:i + 3])
    return codons


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('T')     # these were added for test for all cases
    'A'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')     #
    'C'
    >>> get_complement('*')     # function catches anything which is not ATCG
    '-'

    # - makes the error easily visible in a string
    """
    # DONE: completes doctest

    # if nucleotide == 'A':
    #     return 'T'
    # elif nucleotide == 'T':
    #     return 'A'
    # elif nucleotide == 'G':
    #     return 'C'
    # elif nucleotide == 'C':
    #     return 'G'
    # else:
    #     return '-'

    if nucleotide in nucleotide_complements:
        return nucleotide_complements[nucleotide]
    else:
        return '-'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("CCGCGDTTCA")    # check with error in input dna
    'TGAA-CGCGG'
    """
    # DONE: completes doctest
    complement = ""

    for i in xrange(len(dna)):
        complement_base = get_complement(dna[i])
        complement += complement_base
    return complement[::-1]


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start_codon
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        added tests: for TAG which is out of place, and multiple stop codons

    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGATTAGG")
    'ATGAGATTAGG'
    >>> rest_of_ORF("ATGAGATAGTGA")
    'ATGAGA'
    """
    # DONE: completes doctest

    dna_codons = get_codons(dna)

    for codon in dna_codons:
        if codon in stop_codons:
            return dna[:dna_codons.index(codon) * 3]
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        added test: dna with nested start codon, and with non-default frame

    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATATGGAATGTAGATAGATGTGCCC")
    ['ATGCATATGGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("TATGCATATGGAATGGGGTAGATAGATGTGCCC")
    ['ATGGGG']
    """
    # DONE: completes doctest
    orfs_list = []
    codon_index = 0
    dna_codons = get_codons(dna)

    while codon_index < len(dna_codons):
        if dna_codons[codon_index] == start_codon:
            an_orf = rest_of_ORF(dna[codon_index * 3:])
            orfs_list.append(an_orf)
            codon_index += len(an_orf) / 3
        else:
            codon_index += 1

    return orfs_list

    # while dna:
    #     dna_codons = get_codons(dna)
    #     if dna_codons[index] == start_codon:
    #         an_orf = rest_of_ORF(dna)
    #         orfs_list.append(an_orf)
    #         dna = dna[len(an_orf)::]
    #     else:
    #         dna = dna[3::]

    # return orfs_list


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("-ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("--ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # DONE: completes doctest
    # all_orfs = []

    for offset in xrange(3):
        return sum([find_all_ORFs_oneframe(dna[offset:]) for offset in [0, 1, 2]], [])

    # for offset in xrange(3):
    #     some_orfs = find_all_ORFs_oneframe(offset * "-" + dna)
    #     all_orfs.extend(some_orfs)
    # all_orfs.sort(key=len, reverse=True)  # for test
    # return all_orfs


def find_all_reverse_ORFs(dna):
    reverse_dna = get_reverse_complement(dna)
    return find_all_ORFs(reverse_dna)


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        no more tests added because of extensive testing at lower levels, and
        because this is a basic function
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # DONE: completes doctest
    return sum([find_all_ORFs(strand) for strand in [dna, get_reverse_complement(dna)]], [])
    # much mo' better
    # return find_all_ORFs(dna) + find_all_reverse_ORFs(dna)


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    # doctest.testmod()
    doctest.run_docstring_examples(
        find_all_ORFs_both_strands, globals(), verbose=True)
