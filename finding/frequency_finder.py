BASES = ('A', 'C', 'G', 'T')

def DNA_to_num(DNA: str) -> int:
    num = 0
    chars = len(DNA)
    print(chars)
    for i in range(chars):
        print(DNA[i])
        try:
            num += BASES.index(DNA[i]) * (4 ** (chars - i - 1))
        except ValueError:
            print('Invalid base "' + DNA[i] + '" is ignored')
        print(num)
    return num

def num_to_DNA(num: int, chars: int) -> str:
    DNA = ''
    for i in range(chars):
        print(BASES[int(num % 4)])
        DNA = BASES[int(num % 4)] + DNA
        num = num // 4
        print(num)
    return DNA

if __name__ == '__main__':
    pats = ('A', 'C', 'G', 'T', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
            'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT')
    for pat in pats:
        print(DNA_to_num(pat))
    
