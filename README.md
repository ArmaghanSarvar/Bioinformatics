## Bioinformatics
***This repository includes the implementations of the following Algorithms:***
### Pairwise Sequence Alignment
Global Sequence Alignment of two input sequences using affine gap penalty function.
(Needleman-Wunch dynamic programming algorithm)

### Multiple Sequence Alignment
Sequence Alignment of the input sequences using the **Star** progressive alignment algorithm.

### Profile HMM for Sequence Alignment
Given the number of input sequences, followed by a column gaps threshold θ, followed by a multiple sequence alignment (either Protein or DNA sequences), followed by a test sequence,
* The transition and emission probabilities of the profile HMM(Alignment, θ) are initiated.
* The forward and backward algorithms are implemented.
* The Hidden Markov model is trained using the Baum-Welch learning algorithm and the input MSA.
* The most likely sequence of hidden states (Viterbi path) is found using the Viterbi algorithm and the aligned output sequence is finally returned. 

### Drawing Phylogenetic trees 
Using Biopython Phylo Package to work with Phylogenetic trees 
