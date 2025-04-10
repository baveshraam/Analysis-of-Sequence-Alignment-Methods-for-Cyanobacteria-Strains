import sys
from Bio import SeqIO
import pairwise_alignment as pair_mod
import MSA_Star as star_mod
import MSA_PSSM as pssm_mod
import semiglobal as semi_mod

# Scoring parameters for nucleotide sequences
match_score = 2
mismatch_score = -1
gap_penalty = -2
MAX_LINES = 20  # Limit to first 300 lines of sequence data per file

def load_sequences(filenames):
    """Load sequences from FASTA files, taking the first 300 lines per record."""
    sequences = []
    for filename in filenames:
        try:
            seq_data = []
            with open(filename, 'r') as file:
                lines = file.readlines()
                current_seq = ""
                line_count = 0
                for line in lines:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_seq:  # Save previous sequence
                            seq_data.append(current_seq)
                        current_seq = ""
                        line_count = 0
                    elif line and line_count < MAX_LINES:  # Append sequence data
                        current_seq += line
                        line_count += 1
                        if line_count >= MAX_LINES:
                            break
                if current_seq:  # Save last sequence
                    seq_data.append(current_seq)
            if not seq_data:
                print(f"Warning: No sequences found in {filename}.")
            else:
                print(f"Loaded {len(seq_data)} sequence(s) from {filename}.")
            sequences.extend(seq_data)
        except FileNotFoundError:
            print(f"Error: File '{filename}' not found in the current directory.")
        except Exception as e:
            print(f"Error: Failed to read '{filename}': {e}")
    return sequences

def compute_bayesian_posterior(aligned_seqs, prior=0.1):
    """Simplified Bayesian posterior computation for conserved regions."""
    if not aligned_seqs:
        return []
    length = len(aligned_seqs[0])
    posteriors = []
    for pos in range(length):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0}
        for seq in aligned_seqs:
            if pos < len(seq):  # Make sure we don't go out of bounds
                base = seq[pos]
                if base in counts:
                    counts[base] += 1
                else:  # Handle non-standard bases
                    counts['-'] += 1
        total = sum(counts.values())
        if total == 0:
            posteriors.append(0)
            continue
        conservation = max(counts.values()) / total
        posterior = (conservation * (1 - prior) + prior)  # Bayesian update
        posteriors.append(posterior)
    return posteriors

def train_hmm(aligned_seqs):
    """Simplified HMM training simulation."""
    if not aligned_seqs:
        return None
    states = ['M', 'I', 'D']  # Match, Insert, Delete
    hmm_model = {'transitions': {}, 'emissions': {}}
    for state in states:
        hmm_model['transitions'][state] = {s: 0.33 for s in states}
    
    if not aligned_seqs or len(aligned_seqs) == 0:
        return hmm_model
        
    length = len(aligned_seqs[0])
    for pos in range(length):
        emissions = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0}
        for seq in aligned_seqs:
            if pos < len(seq):  # Make sure we don't go out of bounds
                base = seq[pos]
                if base in emissions:
                    emissions[base] += 1
                else:  # Handle non-standard bases
                    emissions['-'] += 1
        total = sum(emissions.values())
        if total > 0:
            for base in emissions:
                emissions[base] = emissions[base] / total
        hmm_model['emissions'][pos] = emissions
    return hmm_model

def main():
    # Get input files from user
    print("Please provide FASTA file(s) containing sequences.")
    fasta_input = input("Enter sequence file names (comma-separated, e.g., sequence.fasta, sequence2.fasta): ")
    if not fasta_input.strip():
        print("Error: No filenames provided. Please enter at least one filename.")
        sys.exit(1)
    
    fasta_files = [f.strip() for f in fasta_input.split(",") if f.strip()]
    if not fasta_files:
        print("Error: No valid filenames provided.")
        sys.exit(1)

    # Load sequences
    sequences = load_sequences(fasta_files)
    num_sequences = len(sequences)
    if num_sequences == 0:
        print("Error: No valid sequences loaded. Please check your FASTA files and their locations.")
        sys.exit(1)
    if num_sequences < 2:
        print("Error: At least two sequences are required for alignment. Check your FASTA files.")
        sys.exit(1)
    print(f"\nTotal sequences loaded: {num_sequences}")

    # Display options and get user choice
    print("\nSelect an alignment method:")
    print("1. Needleman-Wunsch (Global Alignment)")
    print("2. Smith-Waterman (Local Alignment)")
    print("3. Star Alignment (Multiple Sequence Alignment)")
    print("4. PSSM Analysis")
    print("5. Semi-global Alignment")
    print("6. Bayesian Conservation Analysis")  # New option
    print("7. HMM Training")                    # New option
    print("8. All of the above")                # Updated all option
    choice_input = input("Enter your choice(s) (e.g., 1 or 1,3 or 8): ")

    # Parse and validate choices
    try:
        choices = [int(c.strip()) for c in choice_input.split(",")]
        if not all(1 <= c <= 8 for c in choices):
            raise ValueError
    except ValueError:
        print("Error: Invalid choice. Please enter numbers between 1 and 8, separated by commas.")
        sys.exit(1)

    msa_aligned_seqs = None  # Store MSA results for reuse

    # Execute selected analyses
    for choice in sorted(set(choices)):  # Remove duplicates and sort for consistency
        try:
            if choice == 1 or choice == 8:
                print("\nPerforming Needleman-Wunsch Global Alignment...")
                for i in range(num_sequences):
                    for j in range(i + 1, num_sequences):
                        seq1, seq2 = sequences[i], sequences[j]
                        if not seq1 or not seq2:
                            print(f"Warning: Pair {i+1}-{j+1} skipped due to empty sequence.")
                            continue
                        aligned1, aligned2, score = pair_mod.needleman_wunsch(seq1, seq2, pair_mod.scoring_matrix, gap_penalty)
                        print(f"Pair {i+1}-{j+1}: Score = {score}")
                        print(f"\n Seq1: {aligned1}")
                        print(f"\n Seq2: {aligned2}")

            if choice == 2 or choice == 8:
                print("\nPerforming Smith-Waterman Local Alignment...")
                for i in range(num_sequences):
                    for j in range(i + 1, num_sequences):
                        seq1, seq2 = sequences[i], sequences[j]
                        if not seq1 or not seq2:
                            print(f"Warning: Pair {i+1}-{j+1} skipped due to empty sequence.")
                            continue
                        aligned1, aligned2, score = pair_mod.smith_waterman(seq1, seq2, pair_mod.scoring_matrix, gap_penalty)
                        print(f"Pair {i+1}-{j+1}: Score = {score}")
                        print(f"\n Seq1: {aligned1}")
                        print(f"\n Seq2: {aligned2}")

            if choice == 3 or choice == 8:
                print("\nPerforming Star Alignment (MSA)...")
                msa_aligned_seqs, msa_score = star_mod.calc_MSA_seqs(sequences)
                if not msa_aligned_seqs:
                    print("Warning: MSA failed to produce alignments.")
                    continue
                print(f"MSA Score: {msa_score}")
                for idx, seq in enumerate(msa_aligned_seqs, 1):
                    print(f"\n Seq{idx}: {seq}")

            if choice == 4 or choice == 8:
                print("\nPerforming PSSM Analysis...")
                try:
                    if msa_aligned_seqs is None:
                        msa_aligned_seqs, _ = star_mod.calc_MSA_seqs(sequences)
                    if not msa_aligned_seqs:
                        print("Warning: Cannot perform PSSM: MSA produced no alignments.")
                        continue
                    
                    # Create PSSM profile safely
                    PSSM_profile = {}
                    try:
                        # Create position-based PSSM profile
                        for seq in msa_aligned_seqs:
                            for i in range(len(seq)):
                                if i not in PSSM_profile:
                                    PSSM_profile[i] = {}
                                amino = seq[i]
                                if amino not in PSSM_profile[i]:
                                    PSSM_profile[i][amino] = 1
                                else:
                                    PSSM_profile[i][amino] += 1
                        
                        # Calculate scores with safety checks for log
                        scored_profile = {}
                        pseudocount = 0.1  # Small pseudocount
                        for pos in range(len(msa_aligned_seqs[0])):
                            if pos not in PSSM_profile:
                                continue
                            scored_profile[pos] = {}
                            total = sum(PSSM_profile[pos].values())
                            for amino, count in PSSM_profile[pos].items():
                                # Use pseudocount to avoid zeros
                                frequency = (count + pseudocount) / (total + pseudocount * len(PSSM_profile[pos]))
                                # Use a safer calculation - avoid log of zero
                                background = 0.25  # Equal background for nucleotides
                                if frequency > 0:
                                    import math
                                    try:
                                        scored_profile[pos][amino] = round(math.log2(frequency/background), 3)
                                    except (ValueError, OverflowError):
                                        scored_profile[pos][amino] = -5.0  # Default low score
                                else:
                                    scored_profile[pos][amino] = -5.0  # Default low score
                        
                        print("PSSM Profile (first 5 positions):")
                        for pos in range(min(5, len(msa_aligned_seqs[0]))):
                            print(f"Pos {pos}:")
                            if pos in scored_profile:
                                for amino, score in scored_profile[pos].items():
                                    print(f"  {amino}: {score}")
                            else:
                                print("  No data available")
                                
                    except Exception as e:
                        print(f"Error in PSSM calculation: {e}")
                        import traceback
                        traceback.print_exc()
                except Exception as e:
                    print(f"PSSM Analysis failed: {e}")
                    import traceback
                    traceback.print_exc()

            if choice == 5 or choice == 8:
                print("\nPerforming Semi-global Alignment on first two sequences...")
                if num_sequences < 2:
                    print("Warning: Not enough sequences for semi-global alignment.")
                    continue
                seq1, seq2 = sequences[0], sequences[1]
                if not seq1 or not seq2:
                    print("Warning: Skipping semi-global alignment due to empty sequence.")
                    continue
                aligned1, aligned2, score = semi_mod.main(seq1, seq2, semi_mod.scoring_matrix, gap_penalty)
                print(f"Score: {score}")
                print(f"\n Seq1: {aligned1}")
                print(f"\n Seq2: {aligned2}")
                
            if choice == 6 or choice == 8:  # Bayesian analysis
                print("\nPerforming Bayesian Conservation Analysis...")
                if msa_aligned_seqs is None:
                    msa_aligned_seqs, _ = star_mod.calc_MSA_seqs(sequences)
                if not msa_aligned_seqs:
                    print("Warning: Cannot perform Bayesian analysis: MSA produced no alignments.")
                    continue
                
                posteriors = compute_bayesian_posterior(msa_aligned_seqs)
                print("Bayesian Posterior Probabilities (first 10 positions):")
                for i, prob in enumerate(posteriors[:10]):
                    print(f"Position {i}: Conservation probability = {prob:.3f}")
                
                # Identify highly conserved regions
                conserved_threshold = 0.7  # You can adjust this threshold
                conserved_regions = []
                current_region = []
                
                for i, prob in enumerate(posteriors):
                    if prob >= conserved_threshold:
                        current_region.append(i)
                    elif current_region:
                        conserved_regions.append((min(current_region), max(current_region)))
                        current_region = []
                
                if current_region:  # Don't forget the last region
                    conserved_regions.append((min(current_region), max(current_region)))
                
                print(f"\nFound {len(conserved_regions)} highly conserved regions (conservation > {conserved_threshold}):")
                for start, end in conserved_regions:
                    region_length = end - start + 1
                    print(f"Region from position {start} to {end} (length: {region_length})")

            if choice == 7 or choice == 8:  # HMM training
                print("\nPerforming HMM Training...")
                if msa_aligned_seqs is None:
                    msa_aligned_seqs, _ = star_mod.calc_MSA_seqs(sequences)
                if not msa_aligned_seqs:
                    print("Warning: Cannot train HMM: MSA produced no alignments.")
                    continue
                
                hmm_model = train_hmm(msa_aligned_seqs)
                
                # Display HMM model parameters
                print("HMM Model Parameters:")
                print("Transition Probabilities:")
                for state, transitions in hmm_model['transitions'].items():
                    print(f"  From state {state}:")
                    for next_state, prob in transitions.items():
                        print(f"    To state {next_state}: {prob:.2f}")
                
                print("\nEmission Probabilities (first 5 positions):")
                for pos in range(min(5, len(msa_aligned_seqs[0]))):
                    print(f"Position {pos}:")
                    if pos in hmm_model['emissions']:
                        for base, prob in hmm_model['emissions'][pos].items():
                            print(f"  {base}: {prob:.2f}")
                    else:
                        print("  No data available")

        except Exception as e:
            print(f"Error during analysis {choice}: {e}")
            import traceback
            traceback.print_exc()
            continue

    print("\nAnalysis completed.")

if __name__ == "__main__":
    main()