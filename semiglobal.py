import numpy as np

scoring_matrix = {
    'A': {'A': 2, 'T': -1, 'G': -1, 'C': -1, '-': -2},
    'T': {'A': -1, 'T': 2, 'G': -1, 'C': -1, '-': -2},
    'G': {'A': -1, 'T': -1, 'G': 2, 'C': -1, '-': -2},
    'C': {'A': -1, 'T': -1, 'G': -1, 'C': 2, '-': -2},
    '-': {'A': -2, 'T': -2, 'G': -2, 'C': -2, '-': 0}
}

def do_semi_global_alignment(seq1, seq2, scoring_matrix, gap_penalty=-2):
    n, m = len(seq1), len(seq2)
    alignment_matrix = [[0 for _ in range(m+1)] for _ in range(n+1)]
    paths = [[[] for _ in range(m)] for _ in range(n)]
    for i in range(n+1):
        for j in range(m+1):
            if i == 0 or j == 0:
                alignment_matrix[i][j] = 0
            else:
                lambdaa = scoring_matrix[seq1[i-1]][seq2[j-1]]
                options = [
                    alignment_matrix[i-1][j-1] + lambdaa,
                    alignment_matrix[i-1][j] + gap_penalty,
                    alignment_matrix[i][j-1] + gap_penalty
                ]
                alignment_matrix[i][j] = max(options)
                max_index = [k for k, v in enumerate(options) if v == alignment_matrix[i][j]]
                paths[i-1][j-1] = max_index
    return alignment_matrix, paths

def conv_path_to_seq(alignment_matrix, paths, seq1, seq2, n, m):
    last_row = alignment_matrix[n]
    last_col = [row[m] for row in alignment_matrix]
    max_last_row = max(last_row)
    max_last_col = max(last_col)
    if max_last_row >= max_last_col:
        j = last_row.index(max_last_row)
        start_i, start_j = n, j
    else:
        i = last_col.index(max_last_col)
        start_i, start_j = i, m
    score = alignment_matrix[start_i][start_j]
    align1, align2 = '', ''
    i, j = start_i, start_j
    while i > 0 and j > 0:
        if 0 in paths[i-1][j-1]:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif 1 in paths[i-1][j-1]:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
    if i > 0:
        align1 = seq1[:i] + align1
        align2 = '-'*i + align2
    elif j > 0:
        align1 = '-'*j + align1
        align2 = seq2[:j] + align2
    return align1, align2, score

def main(_seq_1, _seq_2, scoring_matrix, gap_penalty=-2):
    seq1, seq2 = _seq_1, _seq_2
    if len(seq1) < len(seq2):
        seq1, seq2 = seq2, seq1
        switch_mark = True
    else:
        switch_mark = False
    n, m = len(seq1), len(seq2)
    alignment_matrix, paths = do_semi_global_alignment(seq1, seq2, scoring_matrix, gap_penalty)
    align1, align2, score = conv_path_to_seq(alignment_matrix, paths, seq1, seq2, n, m)
    if switch_mark:
        align1, align2 = align2, align1
    return align1, align2, score

if __name__ == "__main__":
    seq1 = "AGTACGCA"
    seq2 = "TATGC"
    align1, align2, score = main(seq1, seq2, scoring_matrix)
    print(f"Score: {score}")
    print(align1)
    print(align2)