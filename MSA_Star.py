import copy
import math

S_MATCH = 3
S_MISMATCH = -1
S_GAP = -2
S_GAP_GAP = 0

def global_align(x, y, s_match, s_mismatch, s_gap):
    A = [[0] * (len(x) + 1) for _ in range(len(y) + 1)]
    for i in range(len(y) + 1):
        A[i][0] = s_gap * i
    for j in range(len(x) + 1):
        A[0][j] = s_gap * j
    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            if y[i-1] == '-' or x[j-1] == '-':
                match_score = s_gap
            else:
                match_score = s_match if y[i-1] == x[j-1] else s_mismatch
            A[i][j] = max(
                A[i-1][j-1] + match_score,
                A[i-1][j] + s_gap,
                A[i][j-1] + s_gap
            )
    align_X, align_Y = "", ""
    i, j = len(y), len(x)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and A[i][j] == A[i-1][j-1] + (s_match if (y[i-1] == x[j-1] and y[i-1] != '-') else (s_mismatch if (y[i-1] != x[j-1] and y[i-1] != '-' and x[j-1] != '-') else s_gap)):
            align_X = x[j-1] + align_X
            align_Y = y[i-1] + align_Y
            i -= 1
            j -= 1
        elif i > 0 and A[i][j] == A[i-1][j] + s_gap:
            align_X = '-' + align_X
            align_Y = y[i-1] + align_Y
            i -= 1
        else:
            align_X = x[j-1] + align_X
            align_Y = '-' + align_Y
            j -= 1
    return align_X, align_Y, A[len(y)][len(x)]

def calc_MSA_score(_seqs):
    score = 0
    for i in range(len(_seqs[0])):
        for j in range(len(_seqs) - 1):
            for k in range(j+1, len(_seqs)):
                if _seqs[j][i] == _seqs[k][i] and _seqs[j][i] != '-':
                    score += S_MATCH
                elif _seqs[j][i] == _seqs[k][i] and _seqs[j][i] == '-':
                    score += S_GAP_GAP
                elif _seqs[j][i] != _seqs[k][i] and (_seqs[j][i] == '-' or _seqs[k][i] == '-'):
                    score += S_GAP
                else:
                    score += S_MISMATCH
    return score

def calc_MSA_seqs(_seqs):
    seq_num = len(_seqs)
    if seq_num < 2:
        return _seqs, calc_MSA_score(_seqs)
    # Pairwise similarities
    pairwise_similarities_matrix = [[0 for _ in range(seq_num)] for _ in range(seq_num)]
    for i in range(seq_num):
        for j in range(i + 1, seq_num):
            _, _, score = global_align(_seqs[i], _seqs[j], S_MATCH, S_MISMATCH, S_GAP)
            pairwise_similarities_matrix[i][j] = score
            pairwise_similarities_matrix[j][i] = score
    sum_score_seqs = [sum(row) for row in pairwise_similarities_matrix]
    center_index = sum_score_seqs.index(max(sum_score_seqs))
    center_seq = _seqs[center_index]
    # Align all sequences to center
    aligned_pairs = []
    for i in range(seq_num):
        if i != center_index:
            align_center, align_other, _ = global_align(center_seq, _seqs[i], S_MATCH, S_MISMATCH, S_GAP)
            aligned_pairs.append((align_center, align_other))
    # Pad to max length
    if aligned_pairs:
        max_len = max(len(pair[0]) for pair in aligned_pairs)
        padded_pairs = []
        for align_center, align_other in aligned_pairs:
            pad_length = max_len - len(align_center)
            padded_center = align_center + '-' * pad_length
            padded_other = align_other + '-' * pad_length
            padded_pairs.append((padded_center, padded_other))
    else:
        max_len = len(center_seq)
        padded_pairs = []
    # Build MSA_seqs
    MSA_seqs = ['' for _ in range(seq_num)]
    MSA_seqs[center_index] = padded_pairs[0][0] if padded_pairs else center_seq
    pair_counter = 0
    for i in range(seq_num):
        if i != center_index:
            MSA_seqs[i] = padded_pairs[pair_counter][1]
            pair_counter += 1
    return MSA_seqs, calc_MSA_score(MSA_seqs)

if __name__ == "__main__":
    seqs = ["AGTACGCA", "TATGC", "AGTAGCA"]
    MSA_seqs, MSA_score = calc_MSA_seqs(seqs)
    print("MSA Score:", MSA_score)
    for seq in MSA_seqs:
        print(seq)