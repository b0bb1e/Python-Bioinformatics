from reference import *

def translate_RNA(RNA: str) -> str:
    """Translates an RNA string into a protein

    :param RNA: the RNA string to translate
    :type RNA: str
    :returns: the transcribed protein string
    :rtype: str
    """
    
    protein = ''
    for i in range(0, len(RNA), 3):
        protein += RNA_to_amino[RNA[i:i + 3]]
    return protein

def translate_DNA(DNA: str, DNA_len: int=-1) -> str:
    """Translates an DNA string into a protein

    :param DNA: the DNA string to translate
    :type DNA: str
    :param DNA_len: the length of the DNA string
    :type DNA_len: int
    :returns: the transcribed protein string
    :rtype: str
    """

    if DNA_len < 0:
        DNA_len = len(DNA)
    protein = ''
    for i in range(0, DNA_len, 3):
        protein += DNA_to_amino[DNA[i:i + 3]]
    return protein

def rev_comp(pat: str) -> str:
    """Finds the reverse complement of a DNA string

    :param pat: the DNA string to REVerse COMPlement
    :type pat: str
    :returns: the reverse complement of pat
    :rtype: str
    """

    if not pat:
        raise ValueError('Cannot reverse-complement empty string')
    rev = ''
    for base in pat:
        try:
            rev = COMP.get(base) + rev
        except TypeError:
            raise ValueError('Non-DNA base "' + base + '" in given string')
    return rev

def find_hidden_proteins(DNA: str, protein: str, allow_rev_comp: bool) -> list:
    """Finds substrings of DNA that encode a protien

    :param DNA: the DNA string to look through
    :type DNA: str
    :param protein: the protein string to look for
    :type protein: str
    :param allow_rev_comp: whether reverse complementary strings should
                           be considered
    :type allow_rev_comp: bool
    :returns: all substrings in DNA that encode protein
    :rtype: list
    """

    subs = []
    DNA_len = len(DNA)
    sub_len = len(protein) * 3
    # loop for each codon frame separetly
    for frame in range(0, 3):
        start = frame
        while start < DNA_len - sub_len + 1:
            trans = translate_DNA(DNA[start:start + sub_len], sub_len)
            if trans == protein:
                subs.append(DNA[start:start + sub_len])
            else:
                # check acids past first to see if extra jumps are needed
                for amino in range(1, len(trans)):
                    if trans[amino] == protein[0]:
                        break
                    # if the first acid doesn't match, skip forward extra
                    else:
                        start += 3
            start += 3
            
    if allow_rev_comp:
        rev_matches = find_hidden_proteins(rev_comp(DNA), protein, False)
        for rev in rev_matches:
            subs.append(rev_comp(rev))
            
    return subs

if __name__ == '__main__':
    with open('data.txt') as data:
        DNA = data.readline().rstrip()
        protein = data.readline().rstrip()
    subs = find_hidden_proteins(DNA, protein, True)
    with open('output.txt', mode='w') as output:
        for sub in subs:
            output.write(sub)
            output.write('\n')
