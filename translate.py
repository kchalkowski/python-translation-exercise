#!/usr/bin/env python3

import sys

#if __name__ == '__main__':
#     rna_seq = ("AUG", "UAC", "UGG", "CAC", "GCU", "ACU", "GCU", "CCA", "UAU", "ACU", "CAC", "CAG", "AAU", "AUC", "AGU", "ACA", "GCG")
#    genetic_code = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
#    "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
#    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
#    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
#    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
#    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
#    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
#    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
#    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
#    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
#    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
#    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
#    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
#    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
#    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
#    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

#rna_seq = ("AUG", "UAC", "UGG", "CAC", "GCU", "ACU", "GCU", "CCA", "UAU", "ACU", "CAC", "CAG", "AAU", "AUC", "AGU", "ACA", "GCG")
#genetic_code = 
#    {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
#    "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
#    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
#    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
#    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
#    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
#    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
#    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
#    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
#    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
#    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
#    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
#    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
#    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
#    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
#    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.
    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.
    """
    sequence = rna_sequence.upper()
    genetic_code = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L", "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP", "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P", "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M", "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T", "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R", "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A", "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
    rna_list = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

    if len(sequence) < 3:
        return ""
    else:
        for i in rna_list:
            if rna_list[0] == "UAA" or rna_list[0] == "UAG" or rna_list[0] == "UGA":
                return ""
            else:
                while rna_list[i] != "UAA" or "UAG" or "UGA":
                    return genetic_code[i]
#for i in rna_list:
#        if len(sequence) < 3 or rna_list[0] == "UAA" or rna_list[0] == "UAG" or rna_list[0] == "UGA":
#            return ""
#        else:
#            return genetic_code[i]

pass


#####################################################################################################################
def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.
    """
    sequence = rna_sequence.upper
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

    rna_1 = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    rna_2 = [sequence[i:i+3] for i in range(1, len(sequence), 3)]
    rna_3 = [sequence[i:i+3] for i in range(2, len(sequence), 3)]

    if "AUG" not in rna_1:
        return ""
    else:
        for i in rna_1[rna_1.index("AUG")::]:
            return genetic_code[i]
    if "AUG" not in rna_2:
        return ""
    else:
        for i in rna_2[rna_2.index("AUG")::]:
            return genetic_code[i]
    if "AUG" not in rna_3:
        return ""
    else:
        for i in rna_3[rna_3.index("AUG")::]:
            return genetic_code[i]

pass

#####################################################################################################################
def get_reverse(rna_sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.
    """
    sequence = rna_sequence.upper()
    rev_seq = sequence[::-1]

    if len(sequence) == 0:
        return ""
    else:
        return rev_seq

pass

#####################################################################################################################
def get_complement(rna_sequence):
    """Get the complement of `sequence`.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.
    """
#    sequence = rna_sequence.upper()
#    rna_list = [sequence[i:i+1] for i in range(0, len(sequence), 3)]
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
#
#    if len(sequence) == 0:
#        return ""
#    else:
#        for i in rna_list:
#            return complement[i]

#jamies code
    comp_seq = ""
    for c in rna_sequence:
        comp_seq += complement[c]
        return comp_seq
pass

#####################################################################################################################
def reverse_and_complement(rna_sequence):
    """Get the reversed and complemented form of `sequence`.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.
    """

    comp = get_complement(rna_sequence)
    return get_reverse(comp)

pass

#####################################################################################################################
def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (three reading frames of the
    current orientation, and the reversed and complemented form) and return (as
    a string) the longest sequence of amino acids that it encodes, according to
    the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty list is returned.
    """

    sequence = rna_sequence.upper()
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

    rna_a = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    rna_b = [sequence[i:i+3] for i in range(1, len(sequence), 3)]
    rna_c = [sequence[i:i+3] for i in range(2, len(sequence), 3)]

    rev_rna_a = get_reverse(rna_1)
    rev_rna_b = get_reverse(rna_2)
    rev_rna_c = get_reverse(rna_3)

    rna_1 = get_all_translations(rna_a, genetic_code)
    rna_2 = get_all_translations(rna_b, genetic_code)
    rna_3 = get_all_translations(rna_c, genetic_code)
    rev_rna_1 = get_all_translations(rev_rna_a, genetic_code)
    rev_rna_2 = get_all_translations(rev_rna_b, genetic_code)
    rev_rna_3 = get_all_translations(rev_rna_c, genetic_code)

    if len(sequence) == 0:
        return ""
    elif len(rna_1) > len(rna_2) and len(rna_3) and len(rev_rna_1) and len(rev_rna_2) and len(rev_rna_3):
        return "".join(rna_1)
    elif len(rna_2) > len(rna_1) and len(rna_3) and len(rev_rna_1) and len(rev_rna_2) and len(rev_rna_3):
        return "".join(rna_2)
    elif len(rna_3) > len(rna_2) and len(rna_1) and len(rev_rna_1) and len(rev_rna_2) and len(rev_rna_3):
        return "".join(rna_3)
    elif len(rev_rna_1) > len(rna_2) and len(rna_3) and len(rna_1) and len(rev_rna_2) and len(rev_rna_3):
        return "".join(rev_rna_1)
    elif len(rev_rna_2) > len(rna_2) and len(rna_3) and len(rev_rna_1) and len(rna_1) and len(rev_rna_3):
        return "".join(rna_1)
    elif len(rev_rna_3) > len(rna_2) and len(rna_3) and len(rev_rna_1) and len(rna_1) and len(rev_rna_2):
        return "".join(rna_1)


#longest_peptide = get_longest_peptide(rna_sequence = rna_seq, genetic_code = genetic_code)

#assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
#message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(rna_seq, longest_peptide)
#sys.stdout.write(message)
#if longest_peptide == "MYWHATAPYTHQNISTA":
#    sys.stdout.write("Indeed.\n")

#    pass


#if __name__ == '__main__':
#    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
#    rna_seq = ("AUG"
#            "UAC"
#            "UGG"
#            "CAC"
#            "GCU"
#            "ACU"
#            "GCU"
#            "CCA"
#            "UAU"
#            "ACU"
#            "CAC"
#            "CAG"
#            "AAU"
#            "AUC"
#            "AGU"
#            "ACA"
#            "GCG")
