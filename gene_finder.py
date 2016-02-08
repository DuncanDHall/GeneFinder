# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Duncan Hall

"""

import random
from amino_acids import aa_table  # , aa, valid_codons
from datetime import datetime
from load import load_seq

dna = load_seq('./data/X73525.fa')
print('started compile  %s' % datetime.now())


def get_shuffled_string(s):
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
    ''' Returns a list of codons corresponding to the dna string passed in
    >>> get_codons('ATGCATGAATGTAG')
    ['ATG', 'CAT', 'GAA', 'TGT']
    '''
    codons = [dna[i:i+3] for i in range(0, len(dna), 3)]
    if len(codons[-1]) != 3:
        return codons[:-1]
    else:
        return codons[:]


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
    >>> get_complement('G')
    'C'
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
    try:
        return nucleotide_complements[nucleotide]
    except:
        raise Exception('DNA input contains invalid character')
        # this isn't working...


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # DONE: completes doctest

    return ''.join([get_complement(nucleotide) for nucleotide in dna[::-1]])

    # complement = ""
    # for i in xrange(len(dna)):
    #     complement_base = get_complement(dna[i])
    #     complement += complement_base
    # return complement[::-1]


def get_rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start_codon
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        added tests: for TAG which is out of place, and multiple stop codons

    >>> get_rest_of_ORF("ATGTGAA")
    'ATG'
    >>> get_rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> get_rest_of_ORF("ATGAGATTAGG")
    'ATGAGATTAGG'
    >>> get_rest_of_ORF("ATGAGATAGTGA")
    'ATGAGA'
    """
    # DONE: completes doctest

    dna_codons = get_codons(dna)

    for codon in dna_codons:
        if codon in stop_codons:
            return dna[:dna_codons.index(codon) * 3]
    return dna


def get_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        added test: dna with nested start codon, and with non-default frame

    >>> get_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> get_all_ORFs_oneframe("ATGCATATGGAATGTAGATAGATGTGCCC")
    ['ATGCATATGGAATGTAGA', 'ATGGAATGTAGA', 'ATGTGCCC']
    >>> get_all_ORFs_oneframe("TATGCATATGGAATGGGGTAGATAGATGTGCCCA")
    ['ATGGGG']
    """
    # DONE: completes doctest
    dna_codons = get_codons(dna)

    return [get_rest_of_ORF(dna[index*3:])
            for index in range(len(dna_codons))
            if dna_codons[index] == start_codon]

    # orfs_list = []
    # for index in range(len(dna_codons)):
    #     if dna_codons[index] == start_codon:
    #         an_orf = get_rest_of_ORF(dna[index*3:])
    #         orfs_list.append(an_orf)
    #         index += len(an_orf) / 3  # also skips the end codon
    # return orfs_list

    # index = 0
    # while index/3 < len(dna_codons):
    #     if dna_codons[index / 3] == start_codon:
    #         an_orf = get_rest_of_ORF(dna[index:])
    #         orfs_list.append(an_orf)
    #         index += len(an_orf)
    #     else:
    #         index += 3

    # return orfs_list

    # while dna:
    #     dna_codons = get_codons(dna)
    #     if dna_codons[index] == start_codon:
    #         an_orf = get_rest_of_ORF(dna)
    #         orfs_list.append(an_orf)
    #         dna = dna[len(an_orf)::]
    #     else:
    #         dna = dna[3::]

    # return orfs_list


def get_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> get_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> get_all_ORFs("-ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> get_all_ORFs("--ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # DONE: completes doctest
    # all_orfs = []

    return sum([get_all_ORFs_oneframe(dna[offset:])
                for offset
                in [0, 1, 2]],
               [])

    # for offset in xrange(3):
    #     some_orfs = get_all_ORFs_oneframe(offset * "-" + dna)
    #     all_orfs.extend(some_orfs)
    # all_orfs.sort(key=len, reverse=True)  # for test
    # return all_orfs


# def get_all_reverse_ORFs(dna):
#     reverse_dna = get_reverse_complement(dna)
#     return get_all_ORFs(reverse_dna)


def get_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        no more tests added because of extensive testing at lower levels, and
        because this is a basic function
    >>> get_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # DONE: completes doctest
    return sum([get_all_ORFs(strand)
               for strand
               in [dna, get_reverse_complement(dna)]], [])
    # much mo' better
    # return get_all_ORFs(dna) + get_all_reverse_ORFs(dna)


def get_longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> get_longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # DONE: completes doctest
    # if shuffle_instance != 0 and shuffle_instance % 10 == 0:
    #     print('finding longest ORF in shuffle %s' % shuffle_instance)
    try:
        return max(get_all_ORFs_both_strands(dna), key=len)
    except:
        pass


def get_longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        >>> get_longest_ORF_noncoding('TGTAGAC', 1000)
        'no correct answer'
        """
    # Done: should be safe as well
    shuffled_dna = [get_shuffled_string(dna) for _ in range(num_trials)]
    try:
        long_ORFS = [get_longest_ORF(dna) for dna in shuffled_dna
                     if get_longest_ORF(dna) is not None]
    except:
        raise RuntimeError('shuffle produced no ORFs. Solution: re-run')
    # only records shuffles which result in ORFs
    return max(long_ORFS, key=len)


def get_AA_from_coding_strand(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents a protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> get_AA_from_coding_strand("ATGCGA")
        'MR'
        >>> get_AA_from_coding_strand("ATGCCCGCTTT")
        'MPA'
    """
    # Done: passes doctests

    return ''.join([aa_table[codon] for codon in get_codons(dna)])
    # I have verified that all possible codons are included in aa_table


def get_probable_AA_sequence(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    >>> get_probable_AA_sequence("ATGCCCGCTTT")
    """
    threshold_length = len(get_longest_ORF_noncoding(dna, 1500))
    longest_ORF = get_longest_ORF(dna)
    if len(longest_ORF) >= threshold_length:
        return get_AA_from_coding_strand(longest_ORF)
    else:
        return 'No one gene found with confidence.'

if __name__ == "__main__":
    start_time = datetime.now()
    print(get_probable_AA_sequence(dna))
    print('finished execution in %s' % (datetime.now() - start_time))
    # import doctest
    # # doctest.testmod()
    # doctest.run_docstring_examples(get_reverse_complement,
    #                                globals(), verbose=True,
    #                                name="Jus' Testin'")
