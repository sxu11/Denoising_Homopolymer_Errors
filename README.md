# homopolymer_error


DNA sequencing denoising based on Homopolymer error correction and sequence merging

A homopolymer segment is a string of identical nucleotides appearing in a row, e.g.
‘AAAA’ (homopolymer of ‘A’s of length 4). In 454 sequencing, each homopolymer segment
is called in a single flow indicated by a light signal. The brightness of the light indicates the
length of the homopolymer. When the same nucleotide appears several times in a row, it may
be hard to distinguish the exact brightness of the signal, resulting in erroneous measurement
of the homopolymer length. 
![Alt text](figs/homopolymer_light1.png?raw=true "Optional Title")

Such homopolymer errors lead to extra indels (insertion and
deletion) in the sequence reads and affect the alignments of sequences to the genome data.
The true counts of sequences identified for specific clones will be affected. 

Based on the Levenshtein distance, an intuitive way to handle the more frequent homopolymer
errors is to assign a smaller penalty score (e.g. a positive value between 0 and
1) to the homopolymer-induced indel errors than that of the regular indels (with penalty 1).
However, it is unclear what penalty score should be chosen. Some researches chose
to treat it as normal indel error (assign penalty 1) or completely ignore the homopolymerlength
information rather than use it. For example, some algorithms filter CATAAAG as
CATAG, where the homopolymer length information is lost. This can lead to overestimation
of sequence difference when there is a substitution in the homopolymer segment, for
example TTTC (filtered as TC) and TCTC (filtered as TCTC). Other options include coding
CATAAAG as C1A1T1A3G1 (nucleotide + homopolymer lengths) or as C1A1T1A3G1
(A1 and A3 are considered as different bases). This treatment does not necessarily
decrease the computational complexity of the distance algorithm and the algorithm will have
to be modified substantially for the extra set of characters.

Instead of using only sequence frequency information, we use the following information of
a sequence as criteria to evaluate the quality of the read: mapping type, initial counts
(frequency), and read length. There are three mapping types of sequences obtained by
using BLAT and GMAP algorithm: Single, which means that the sequence was
mapped to a single site on the genome; MultiHits, which means that the sequence was
mapped to multiple genome sites; NoHits, which means that the sequence was not mapped
to any genome site. The initial counts were obtained by simply counting the number of
mapped sequences and those similar to the mapped ones by treating homopolymer indels as
regular indels. The read length of a sequence is between 25 bp (below which the sequence
read is considered unreliable) and 500 bp (above which the sequencing efficiency dropped
significantly). Generally, the longer the read is, the more reliable the sequence is
considered to be.

To calculate the distances between any two sequences, the basic idea is to use a modified
Levenshtein algorithm with an adaptive penalty for the homopolymer indels. One concern
about the Levenshtein-based method is the time complexity; based on the current need
of merging 103 − 104 unmapped sequences to 102 − 103 mapped sequences (representing
several hundreds of distinguishable clones) and some obvious optimization techniques, it is a
realistic approach. We collect indel mismatches among 95% similar sequences from the postalignment
data. This way, we get a rough sense of the occurrence rates for different types
of errors. Although it is likely that some errors are not counted when some genuine copies
are below 95% similarity from the reference segment, relative ratios among different types
of errors should not be affected much. In Figure 2.1(a), the total homopolymer indel rate is
about three times that of the normal indel, we thus set the penalty score for homopolymer
indel as 1/3 . By only focusing on the homopolymer segments of these sequences, we plot
Figure 2.1(b) to show how the error rate increases with the homopolymer length. Intuitively
it is more likely to have homopolymer errors in a segment ‘AAAAA’ than in ‘AA’. When the
homopolymer length is larger than 7, there is more than 50% probability that more than 1
indel errors will occur, so we set the penalty in this region to be 0.

![Alt text](figs/stats.png?raw=true "Optional Title")

We assume two situations that sequences might be similar to each other: First, a sequence is
a ‘genuine’ copy of a genome segment with a small number of mutations, insertions, or deletions.
Second, a sequence is a ‘fake’ copy of a genome segment with a substantial number of
overlapping nucleotides as a result of the functional redundancy of primate genome. Ideally,
‘fake’ sequences from the second group is farther away than the ‘genuine’ ones from the first
group. Our aim is to find a proper threshold that can distinguish the two groups of similar
sequences. For this purpose, we sampled a few genome sequences that were successfully
mapped and plot their distances with all raw sequences in Figure 2.2. It seems that as long
9
as we choose a threshold between 92%−96% (currently we choose 95%), ‘genuine’ and ‘fake’
copies of the genome segments can be well separated.

![Alt text](figs/scores.png?raw=true "Optional Title")

(0) Input filename
(1) Pre-processing e.g. eliminate length < 25 bp sequences, delete reads of ’N’, cut the piece
with linker sequences, etc.
(2) If want to achieve time-efficiency and care only about recovering the counts of mapped
sequences, go to (2.1); otherwise (also recover the counts of unmapped sequences), go to
(2.2).
(2.1) Unsupervised clustering
(2.1.1) Calculate the pair-wise similarities of all sequences
(2.1.2) Group all sequences (e.g. if A=B, B=C, but A6=C, still group {A, B, C}), output
grouping results
(2.1.3) Rank sequences and elect the highest rank one as the representative
(2.1.4) If want to increase sensitivity/decrease specificity, go to (2.1.4.1), otherwise go to

... ... ... ... ...
10 Single - 98 CTAGGAAAACGA−TTATAGCTGCACAAAC−A−CTTT−GT−CTCG...
10 R|NoHits775 149 CTAG−AAA−−GA−TTATAGCTGCACAAACGA−CTTT−GT−CTCG...
10 R|NoHits795 41 CTAGGAAA−−GAATTATAGCTGCACAAAC−AACTTTTGTTCTCG
10 R|NoHits786 154 CTAGGAAAA−GAATTATAGCTGCACAAAC−A−CTTTTGTTCTCG...
... ... ... ... ...
(2.1.4.1.1) Merge all other sequences in a group into the representative one
(2.1.4.1.2) Output merged sequences
(2.1.4.2) exclusive merging Example: If INCLUSIVE is False, we have (NoHits775 got separated
out):
... ... ... ... ...
10 Single - 98 CTAGGAAAACGA−TTATAGCTGCACAAAC−A−CTTT−GT−CTCG...
-10 R|NoHits775 149 CTAG−AAA−−GA−TTATAGCTGCACAAACGA−CTTT−GT−CTCG...
10 R|NoHits795 41 CTAGGAAA−−GAATTATAGCTGCACAAAC−AACTTTTGTTCTCG
10 R|NoHits786 154 CTAGGAAAA−GAATTATAGCTGCACAAAC−A−CTTTTGTTCTCG...
... ... ... ... ...
(2.1.4.2.1) In each group, for those 95% similar to the representative, merge to the representative;
otherwise, separate and mark the unqualified sequences
(2.1.4.2.2) Output sequences
(2.2) Supervised clustering
(2.2.1) Separate master (MultiHits/Single) sequences and slave (NoHits) sequences
(2.2.2) Calculate the pair-wise similarities between masters and slaves
(2.2.3) for those slaves that is similar to a unique master, merge and output
(2.2.4) for those that is similar to multiple masters, if randomly assign to a sequence, go to
(2.2.4.1); if output information, go to (2.2.4.2)
(2.2.4.2) Mark and output
Example: When RANDOM ASSIGN is False, the suspicious NoHits2743 is separated and

marked, merging to neither Single 197 nor Single 238:
... ... ... ... ...
197 Single + 411 GTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTGG...
... ... ... ... ...
238 Single - 112 GTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTGG...
... ... ... ... ...
197 238 R|NoHits2743 39 GTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTCG
... ... ... ... ...
(2.2.4.1) Randomly merge to a similar master, output Example: When RANDOM ASSIGN
is True: NoHits2743 is merged to Single 197.
(2.2.5) for those slaves that is similar to no master, directly output

Homopolymer penalty (use penalty 1/6 as an example). The penalty value for homopolymer
errors
Example: If HOMO PENALTY is 1/6, the distance of the following two Seqs is 1/6+1/6+
1/6 + 1 = 1.5; if HOMO PENALTY is 1 (same as normal indel), the distance is 4.
Seq 1: TATACTTGGGGCTATTT−CACAAT−GGAAATAATTAGCCCGT
Seq 2: TATACTTGGG−CTATTTTCACAATTGGAAATAATT−GCCCGT
End gap (default: False). Whether to consider the end gap when comparing two sequences
Example: If END GAP is False, the following two Seqs are considered the same. Otherwise,
they are not considered the same:
Seq 1: GCCTGCCCCCCGAGCTCTCCCGTGTGGATCCCGCA
Seq 2: GCCTGCCCCCTGAGCTCTCCCGTGT
Similarity Threshold (default: 0.95, range: [0, 1]). The similarity threshold for two
sequences to be considered the same.
Example: (With HOMO PENALTY 1/6, END GAP is False) The distance of the following

two Seqs is 1/6 + 1 + 1 + 1/6 = 2.33, the similarity is 1 − 2.33/40 = 0.942 < 0.95; they are
not considered the same:
Seq 1: GTAG−AGCTTTTAC−ATACTTACAGGCATATGCACAG−CAA−TC
Seq 2: GTAGGAGCTTTTACTATACTTACAGGCATATGCACAGACAAAGT
Score Weights (default: [0.5, 0.3, 0.2], range: 0 ≤ pi ≤ 1,
P
pi = 1). The importance
weights among different rules, including hit type, post-alignment count, and sequence length.
Hit-type Weights (default: [0.70, 0.25, 0.05], range: 0 ≤ pi ≤ 1,
P
pi = 1), The
importance weights among various hit types, including MultiHits, Single, and NoHits.
Clustering strategy (default: exclusive). After the first step of calculating all pairwise
distances between mapped and unmapped sequences, we consider two sequences equivalent
if they are of ≥ 95% similarity with each other. However, this may bring inconsistency when
extending the criterion to group more than two sequences, since essentially this definition
of equivalence is not transitive. Specifically, if we use the equal sign A = B to denote that
two sequences A and B are equivalent, we have a problem when A = B, B = C but A 6= C.
The sources of this inconsistency include the intrinsic redundancy nature of different DNA
segments, the shortness of the reading length of sequences, and the impossibility to perfectly
capture the true variations of different sequences with the same set of parameters (e.g.
similarity threshold) across all sequences. We give the option of considering A=C (inclusive
strategy) and considering A6=C (exclusive strategy) when grouping sequences. These two
options give the upper and lower bounds for the estimate of counts for sequence A.

Performance We checked how many valid VISs were recovered by performing our homopolymer
algorithm. We also checked our results with the eye-balling results obtained by
three lab technicians using three weeks. In comparison, our computation tool used less than
20 minutes on a Thinkpad laptop (with 4GB memory and Intel i5@1.6GHZ CPU).

![Alt text](figs/results.png?raw=true "Optional Title")

