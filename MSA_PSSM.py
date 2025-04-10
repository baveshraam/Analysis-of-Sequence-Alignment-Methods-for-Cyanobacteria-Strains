from Bio import SeqIO
import copy, math

pseudocount = 2
window_size = 0

# Load sequences from FASTA 
fasta_file = "sequence.fasta"
try:
    sequences = [str(record.seq)[:100] for record in SeqIO.parse(fasta_file, "fasta")]
except FileNotFoundError:
    print(f"Error: {fasta_file} not found in the current directory.")
    exit(1)

def fill_profile(PSSM_profile, seqs):
    # First part remains the same
    for seq in seqs:
        for i in range(len(seq)):
            amino = seq[i]
            if i not in PSSM_profile:
                PSSM_profile[i] = {}
            if amino not in PSSM_profile[i]:
                PSSM_profile[i][amino] = 1  # Count occurrence
            else:
                PSSM_profile[i][amino] += 1

    # Calculate score with safety check for log
    scored_profile = {}
    for pos in range(len(seqs[0])):
        if pos not in PSSM_profile:
            continue
        scored_profile[pos] = {}
        total_count = sum(PSSM_profile[pos].values()) + pseudocount * len(PSSM_profile[pos])
        for amino, count in PSSM_profile[pos].items():
            # Add pseudocount to avoid zeros
            prob = (count + pseudocount) / total_count
            # Avoid log of zero or negative numbers
            if prob > 0:
                background = 0.25  # Assuming equal background probability for nucleotides
                scored_profile[pos][amino] = round(math.log(prob/background, 2), 3)
            else:
                scored_profile[pos][amino] = -999.0  # Very low score for impossible cases
    
    return scored_profile

def calc_powerset(fullset):
  listsub = list(fullset)
  subsets = []
  for i in range(2**len(listsub)):
    subset = []
    for k in range(len(listsub)):            
      if i & 1<<k:
        subset.append(listsub[k])
    subsets.append(subset)        
  return subsets

def calc_score(sub_seq, pssm_profile):
    score = 0
    for i in range(len(sub_seq)):
        score += pssm_profile[sub_seq[i]][i]
    return score

def find_best_subseq(target_seq, pssm_profile): 
    best_score = -1*math.inf
    best_sub_target = None
    power_set = calc_powerset(set([i for i in range(window_size)]))
    for gap_num in range(window_size):
        char_num = window_size - gap_num
        for i in range(len(target_seq) - char_num + 1):
            tar_part = target_seq[i:i+char_num]
            # insert gaps
            for sub_p_s in power_set:
                if len(sub_p_s) == gap_num:
                    tar_part_temp = copy.deepcopy(tar_part)
                    for sp in sub_p_s:
                        tar_part_temp = tar_part_temp[:sp] + '-' + tar_part_temp[sp:]
                        score = calc_score(tar_part_temp, pssm_profile)
                        if score > best_score:
                            best_score = score
                            best_sub_target = tar_part_temp
    return best_sub_target, best_score

if __name__=="__main__":
    seqs = sequences[:-1]  # Use all but the last sequence for profile
    target_seq = sequences[-1]  # Use the last sequence as target
    window_size = len(seqs[0])  # Set window size based on sequence length

    PSSM_profile = {}
    pssm_profile = fill_profile(PSSM_profile, seqs)

    best_sub_target, best_score = find_best_subseq(target_seq, pssm_profile)
    print(f"Best Subsequence: {best_sub_target}")
    print(f"Best Score: {best_score}")