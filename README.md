# Screen ssDNA ChIP-chip library for chromosome size bias effects

Our lab has an extensive library of ssDNA ChIP-chip data, used for mapping DSB hotspots
in meiosis. Not all these data have been collected (or even later examined) to detect
changes in chromosome size bias effects in meiotic recombination. This way, it is worth
screening the entire library for such effects.

#### Strategy:
Use a simple statistical test to compare meiotic recombination levels between small and
large chromosomes and then compare the output between different samples. The approach
is the following:
- Convert Log2Ratio signal back to ratio
- Calculate average ChIP-chip signal for each chromosome
- Consider individual chrs as observations in two groups:
	- small chrs (chr I, III, and VI)
	- large chrs (remaining 13 chrs)
- Use a statistical test to compare the two groups

#### For command line usage run:
```
python chr_size_bias_in_ssDNA_chip-chip.py --help
```
