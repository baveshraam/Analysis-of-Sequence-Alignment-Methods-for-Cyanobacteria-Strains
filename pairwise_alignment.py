import numpy as np

scoring_matrix = {
    'A': {'A': 2, 'T': -1, 'G': -1, 'C': -1, '-': -2},
    'T': {'A': -1, 'T': 2, 'G': -1, 'C': -1, '-': -2},
    'G': {'A': -1, 'T': -1, 'G': 2, 'C': -1, '-': -2},
    'C': {'A': -1, 'T': -1, 'G': -1, 'C': 2, '-': -2},
    '-': {'A': -2, 'T': -2, 'G': -2, 'C': -2, '-': 0}
}

def smith_waterman(seq1, seq2, scoring_matrix, gap_penalty=-2):
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((n+1, m+1))
    traceback = np.zeros((n+1, m+1), dtype=int)
    max_score = 0
    max_pos = (0, 0)
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score_matrix[i-1][j-1] + scoring_matrix[seq1[i-1]][seq2[j-1]]
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score = max(0, match, delete, insert)
            score_matrix[i][j] = score
            if score > max_score:
                max_score = score
                max_pos = (i, j)
            traceback[i][j] = np.argmax([0, match, delete, insert])
    align1, align2 = '', ''
    i, j = max_pos
    while traceback[i][j] != 0 and i > 0 and j > 0:
        if traceback[i][j] == 1:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif traceback[i][j] == 2:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
    return align1, align2, max_score

def needleman_wunsch(seq1, seq2, scoring_matrix, gap_penalty=-2):
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((n+1, m+1))
    traceback = np.zeros((n+1, m+1), dtype=int)
    for i in range(1, n+1):
        score_matrix[i][0] = score_matrix[i-1][0] + gap_penalty
    for j in range(1, m+1):
        score_matrix[0][j] = score_matrix[0][j-1] + gap_penalty
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score_matrix[i-1][j-1] + scoring_matrix[seq1[i-1]][seq2[j-1]]
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score = max(match, delete, insert)
            score_matrix[i][j] = score
            traceback[i][j] = np.argmax([match, delete, insert]) + 1
    align1, align2 = '', ''
    i, j = n, m
    while i > 0 or j > 0:
        if traceback[i][j] == 1:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif traceback[i][j] == 2:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
    return align1, align2, score_matrix[n][m]

if __name__ == "__main__":
    seq1 = "AGTACGCA"
    seq2 = "TATGC"
    print("Needleman-Wunsch Global Alignment:")
    a1, a2, score = needleman_wunsch(seq1, seq2, scoring_matrix, gap_penalty=-2)
    print(f"Score: {score}")
    print(a1)
    print(a2, "\n")
    print("Smith-Waterman Local Alignment:")
    a1, a2, score = smith_waterman(seq1, seq2, scoring_matrix, gap_penalty=-2)
    print(f"Score: {score}")
    print(a1)
    print(a2)