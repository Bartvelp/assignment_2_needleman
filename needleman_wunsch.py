#!/usr/bin/env python3

"""
Author: Bart Grosman
Student number: 1036073 

Description: This script implements the Needleman-Wunsch algorithm and 
executes it on a number of example strings. The BLOSUM62 scoring matrix is used
to calculate similarity between amino acids.
Assignment 2 for Algorithms in Bioinformatics, 2021-5-20

Usage: python3 needleman_wunsch.py
"""
# no imports needed
# functions between here and __main__
blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""

def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int,parts[1:])))
    return order, blosum_matrix

BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()

def score(res1, res2):
    """Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]

# My own functions

def init_matrix(m: int, n: int, default_value = 0, row0_default = 0, col0_default = 0):
    """ Gives a matrix of size m (rows) * n (columns), with passed default vals

    Keyword arguments:
    m -- Number of rows in matrix
    n -- Number of columns in matrix
    default_value -- value used to fill entire array
    row0_default -- value used to overwrite the first row
    col0_default -- value used to overwrite first column
    """
    matrix = [[default_value for column in range(n)] for row in range(m)]
    # fill first row
    for col_num in range(n):
        matrix[0][col_num] = row0_default
    # fill first column
    for row_num in range(m):
        matrix[row_num][0] = col0_default
    
    return matrix

def fill_matrix(seq1: str, seq2: str, p_gap: int = -8, end_gap = 0):
    num_row = len(seq1) + 1
    num_col = len(seq2) + 1 # n
    matrix = init_matrix(num_row, num_col, 0)
    arrow_matrix = init_matrix(num_row, num_col, [], [2], [1])

    # fill first row
    for col_num in range(num_col):
        matrix[0][col_num] = end_gap * col_num
    # fill first column
    for row_num in range(num_row):
        matrix[row_num][0] = end_gap * row_num

    # Fill in the rest of the matrix
    for row_n in range(1, num_row):
        for col_n in range(1, num_col):
            # first the diagonal
            seq1_AA = seq1[row_n - 1]
            seq2_AA = seq2[col_n - 1]
            diagonal = matrix[row_n - 1][col_n - 1] + score(seq1_AA, seq2_AA)
            # now the gapped variants
            # Take into account the end gap if moving along the axis
            gap_p_vertical = end_gap if col_n == num_col - 1 else p_gap
            gap_p_horizontal = end_gap if row_n == num_row - 1 else p_gap

            gapped_1 = matrix[row_n - 1][col_n] + gap_p_vertical
            gapped_2 = matrix[row_n][col_n - 1] + gap_p_horizontal
            # Now save the best possible outcome
            arrows, best_option = get_directions([diagonal, gapped_1, gapped_2])
            matrix[row_n][col_n] = best_option
            arrow_matrix[row_n][col_n] = arrows
    
    return matrix, arrow_matrix

def traceback(arrow_matrix, seq1, seq2):
    # We need to start at the bottom right, so get the size
    num_row = len(arrow_matrix)
    num_col = len(arrow_matrix[0])
    x, y = [num_row - 1, num_col - 1]

    alignment = [] # Array of tuples [(A, A), (D, -), (W, W), ...]
    
    while x != 0 or y != 0:
        # Traceback our steps to the origin
        arrow = arrow_matrix[x][y][0] # Always get the first entry
        res1 = seq1[x - 1]
        res2 = seq2[y - 1]
        if arrow == 0: # Diagonal (i == 0)
            alignment.insert(0, (res1, res2))
            x += -1
            y += -1
        elif arrow == 1: # Vertical
            alignment.insert(0, (res1, '-'))
            x += -1
        elif arrow == 2: # horizontal   
            alignment.insert(0, ('-', res2))
            y += -1
    return alignment

# Helper functions
def get_directions(direction_points):
    max_points = max(direction_points)
    arrows = [i for i, direction_point in enumerate(direction_points) 
        if direction_point == max_points]
    # Arrows implemented as indices, 0
    return arrows, max_points

def calc_perc_identity(alignment):
    num_identical_res = 0
    for aligned_pair in aligment:
        if (aligned_pair[0] == aligned_pair[1]):
            num_identical_res += 1
    perc_identity = num_identical_res / len(alignment) * 100
    return perc_identity

def print_matrix(matrix: list):
    for row in matrix:
        print(' '.join(map(str, row)))
    
def print_aligment(alignment):
    res1 = [res[0] for res in alignment]
    res2 = [res[1] for res in alignment]
    print(''.join(res1))
    print(''.join(res2))

if __name__ == "__main__":
    seq1 = "THISLINE"
    seq2 = "ISALIGNEDYYY"
    matrix, arrow_matrix = fill_matrix(seq1, seq2, -5, -10)
    print_matrix(matrix)
    print_matrix(arrow_matrix)
    aligment = traceback(arrow_matrix, seq1, seq2)
    print_aligment(aligment)
    print("Score", matrix[len(seq1)][len(seq2)])
    print("PI", calc_perc_identity(aligment))
    # seq3: GPA1_ARATH
    seq3 = ("MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQ"
    "TGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEKLSEIGGRLDYPRLTKD"
    "IAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYARVRTTGVVEIQFSPVG"
    "ENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPC"
    "FEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRV"
    "FKIYRTTALDQKLVKKTFKLVDETLRRRNLLEA")
    # seq4: GPA1 BRANA
    seq4=(\
    "MGLLCSRSRHHTEDTDENAQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQASS"
    "DKRKIIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDPAKYTLSSEN"
    "MAIGEKLSEIGARLDYPRLTKDLAEGIETLWNDPAIQETCSRGNELQVPDCTKYLMENLK"
    "RLSDVNYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHL"
    "FEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSIMLFLNKFDIF"
    "EKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYR"
    "TTALDQKLVKKTFKLVDETLRRRNLLEAGLL")

    matrix, arrow_matrix = fill_matrix(seq3, seq4, -5, -1)
    print(matrix[len(seq3)][len(seq4)])

    aligment = traceback(arrow_matrix, seq3, seq4)
    print_aligment(aligment)
    print(calc_perc_identity(aligment))
    print("Score", matrix[len(seq3)][len(seq4)])
    print("PI", calc_perc_identity(aligment))
