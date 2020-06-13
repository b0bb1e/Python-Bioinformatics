BASES = ('A', 'C', 'G', 'T')

def DNA_to_num(DNA: str) -> int:
    """Converts a DNA string to a number

    For all strings of a certain length, lower numbers will correspond
    to an earlier point in the alphabet. Numbers are re-used in-between
    lengths.

    :param DNA: the DNA string of ACGT to convert
    :type DNA: str
    :returns: a number corresponding to the DNA string
    :rtype: int
    """
    
    num = 0
    chars = len(DNA)
    for i in range(chars):
        try:
            # earlier chars are given more weight/value in the number
            num += BASES.index(DNA[i]) * (4 ** (chars - i - 1))
        except ValueError:
            print('Invalid base "' + DNA[i] + '" is ignored')
    return num

def num_to_DNA(num: int, chars: int) -> str:
    """Converts a number to a DNA string

    :param num: the number to convert
    :type num: int
    :param chars: the length of the eventual string
    :type chars: int
    :returns: the DNA string of length chars equivalent to num
    :rtype: str
    """
    
    DNA = ''
    for i in range(chars):
        # the smaller numbers, pulled off first, go in the back
        DNA = BASES[int(num % 4)] + DNA
        num = num // 4
    if num != 0:
        raise ValueError('Number is too large for length given')
    return DNA
