# Arbor Biotechnologies Coding Challenge

**Candidate**: Kunal Shah (kunalshah.me@gmail.com)
**Completion Date**: 03/08/2020

This repository contains code and answers in response to the following coding challenge posed by Arbor Biotechnologies:

**Coding Exercise**

This toy problem is framed in the spirit of discovering proteins via database construction and computational search.  We want to download some genomic data, extract proteins, calculate features to describe each protein, and then demonstrate how a query can "find" two instances of Cas9.

1. Download the following files from ncbi:
 - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/239/625/GCF_001239625.1_7068_7_24/GCF_001239625.1_7068_7_24_genomic.gbff.gz
 - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/014/815/GCF_002014815.1_ASM201481v1/GCF_002014815.1_ASM201481v1_genomic.gbff.gz

2. Uncompress them, and look at the contents to get some familiarity with how this (plain text) file format annotates sequence data.  In particular, each gbff file is a concatenation of sequence data, where the beginning of each is indicated by the line that starts with "LOCUS".  Also, you will see a list of "FEATURES" -- predominantly a "gene" or "CDS" (coding sequence).  Next to each coding sequence you can find a list of "qualifiers" such as protein_id and translation

3. Use Biopython SeqIO.parse (format = 'genbank') to read the records from each gbff file. You want to pay special attention to these two records, since they contain both CDS and "repeat_region" (CRISPR) features:
 - NZ_CUFU01000032
 - NZ_MUYP01000003

4. Given any protein (CDS), compute a set of features for it.  The set of features should include the size of the protein (length of the translation string), and the distance to the nearest repeat region (or infinite, if the record doesn't contain any).  Distance can be measured in # nucleotides, or as # intervening features.  Feel free to think of and implement additional features (encouraged).

5. Given a set of proteins and their features, can you determine a set of criteria that filters for two instances of Cas9 (i.e. protein_ids WP_053019794.1 and WP_010922251.1)?  Incidentally, these two proteins are SpCas9 and SaCas9 (https://www.ncbi.nlm.nih.gov/pubmed/25830891). Your filter should be designed with the idea that it could be applied to novel proteins to determine if they are Cas9 as well.

6. Create a summary of your findings. Briefly discuss what features you tried, as well as any data exploration that you did. Also discuss how your solution could be scaled up to production on very large data sets.

## Steps 1-5 Solution

#### Running the solution

Note: This has been tested using Python 3.6

Install dependencies: `pip install -r requirements.txt`

Help menu: `python extract_protein_features.py --help`
Output:
```
Usage: extract_protein_features.py [OPTIONS] FILE_PATHS...

Options:
  --gene         TEXT        Gene (e.g. cas9) to display feature metrics for
  --protein_id TEXT         protein_id (e.g. WP_053019794.1) to display feature
                                     metrics for
  --help             Show this message and exit.
```

Example command specifying two local input files, the gene `cas9`, and two `cas9` proteins:
```
python extract_protein_features.py \
	GCF_001239625.1_7068_7_24_genomic.gbff \
	GCF_002014815.1_ASM201481v1_genomic.gbff \
	--gene cas9 \
	--protein_id WP_010922251.1 \
	--protein_id WP_053019794.1
```
Output (converted from Pandas DF to Markdown):

**Overall metrics** 

|       |   min_dist_to_nearest_repeat |   repeat_seqs_w_overlap_ratio |   protein_len |   most_freq_palindrome_freq |   most_freq_palindrome_len |   largest_palindrome_len |   unique_palindromes_count |
|:------|-----------------------------:|------------------------------:|--------------:|----------------------------:|---------------------------:|-------------------------:|---------------------------:|
| count |                        372   |                           372 |      4066     |                         372 |                        372 |               372        |                 372        |
| mean  |                      83774.4 |                             0 |       295.807 |                           1 |                          4 |                 4.08602  |                   1.08602  |
| std   |                      72807.5 |                             0 |       199.508 |                           0 |                          0 |                 0.406309 |                   0.406309 |
| min   |                         30   |                             0 |        20     |                           1 |                          4 |                 4        |                   1        |
| 25%   |                      22700.8 |                             0 |       153     |                           1 |                          4 |                 4        |                   1        |
| 50%   |                      57804.5 |                             0 |       260     |                           1 |                          4 |                 4        |                   1        |
| 75%   |                     142769   |                             0 |       389     |                           1 |                          4 |                 4        |                   1        |
| max   |                     245136   |                             0 |      2045     |                           1 |                          4 |                 6        |                   3        |

**Protein WP_010922251.1 metrics**

|       |   min_dist_to_nearest_repeat |   repeat_seqs_w_overlap_ratio |   protein_len |   most_freq_palindrome_freq |   most_freq_palindrome_len |   largest_palindrome_len |   unique_palindromes_count |
|:------|-----------------------------:|------------------------------:|--------------:|----------------------------:|---------------------------:|-------------------------:|---------------------------:|
| count |                            1 |                             1 |             1 |                           1 |                          1 |                        1 |                          1 |
| mean  |                         1961 |                             0 |          1368 |                           1 |                          4 |                        4 |                          1 |
| std   |                          nan |                           nan |           nan |                         nan |                        nan |                      nan |                        nan |
| min   |                         1961 |                             0 |          1368 |                           1 |                          4 |                        4 |                          1 |
| 25%   |                         1961 |                             0 |          1368 |                           1 |                          4 |                        4 |                          1 |
| 50%   |                         1961 |                             0 |          1368 |                           1 |                          4 |                        4 |                          1 |
| 75%   |                         1961 |                             0 |          1368 |                           1 |                          4 |                        4 |                          1 |
| max   |                         1961 |                             0 |          1368 |                           1 |                          4 |                        4 |                          1 |

**Protein WP_053019794.1 metrics**

|       |   min_dist_to_nearest_repeat |   repeat_seqs_w_overlap_ratio |   protein_len |   most_freq_palindrome_freq |   most_freq_palindrome_len |   largest_palindrome_len |   unique_palindromes_count |
|:------|-----------------------------:|------------------------------:|--------------:|----------------------------:|---------------------------:|-------------------------:|---------------------------:|
| count |                            1 |                             1 |             1 |                           1 |                          1 |                        1 |                          1 |
| mean  |                         2263 |                             0 |          1053 |                           1 |                          4 |                        6 |                          3 |
| std   |                          nan |                           nan |           nan |                         nan |                        nan |                      nan |                        nan |
| min   |                         2263 |                             0 |          1053 |                           1 |                          4 |                        6 |                          3 |
| 25%   |                         2263 |                             0 |          1053 |                           1 |                          4 |                        6 |                          3 |
| 50%   |                         2263 |                             0 |          1053 |                           1 |                          4 |                        6 |                          3 |
| 75%   |                         2263 |                             0 |          1053 |                           1 |                          4 |                        6 |                          3 |
| max   |                         2263 |                             0 |          1053 |                           1 |                          4 |                        6 |                          3 |

**Gene cas9 metrics**

|       |   min_dist_to_nearest_repeat |   repeat_seqs_w_overlap_ratio |   protein_len |   most_freq_palindrome_freq |   most_freq_palindrome_len |   largest_palindrome_len |   unique_palindromes_count |
|:------|-----------------------------:|------------------------------:|--------------:|----------------------------:|---------------------------:|-------------------------:|---------------------------:|
| count |                        2     |                             2 |         2     |                           2 |                          2 |                  2       |                    2       |
| mean  |                     2112     |                             0 |      1210.5   |                           1 |                          4 |                  5       |                    2       |
| std   |                      213.546 |                             0 |       222.739 |                           0 |                          0 |                  1.41421 |                    1.41421 |
| min   |                     1961     |                             0 |      1053     |                           1 |                          4 |                  4       |                    1       |
| 25%   |                     2036.5   |                             0 |      1131.75  |                           1 |                          4 |                  4.5     |                    1.5     |
| 50%   |                     2112     |                             0 |      1210.5   |                           1 |                          4 |                  5       |                    2       |
| 75%   |                     2187.5   |                             0 |      1289.25  |                           1 |                          4 |                  5.5     |                    2.5     |
| max   |                     2263     |                             0 |      1368     |                           1 |                          4 |                  6       |    

#### Solution Explanation

Based on my understanding of the requirements, I performed the following steps:
1. Used BioPython to load the two specified GenBank files
2. Iterated over each record
3. Iterated over each `CDS` and `repeat_region` feature in each record
4. Computed a set of numeric features from each coding sequence
5. Flattened each collection of features into a single iterrable of type `generator[dict]`
6. Displayed metrics relating to the features using `Pandas.DataFrame.describe()`
7. Offerred limited ability to filter and display metrics by: CDS "gene" (e.g. "Cas9") and CDS "protein_id"  (e.g. "WP_010922251.1")

**Feature table**

| Feature                     | Description                                                                           |
|-----------------------------|---------------------------------------------------------------------------------------|
| min_dist_to_nearest_repeat  | Minimum possible distance to nearest repeat region in #-of-bases-away                 |
| repeat_seqs_w_overlap_ratio | Ratio of repeat sequences that overlap with given CDS                                 |
| protein_len                 | Length of protein translation string                                                  |
| most_freq_palindrome_freq   | Frequency of most frequent palindromic sequence in nearest repeat region              |
| most_freq_palindrome_len    | Length of most frequent palindromic sequence in nearest repeat region (in #-of-bases) |
| largest_palindrome_len      | Length of largest palindromic sequence in nearest repeat region                       |
| unique_palindromes_count    | Number of unique palindromic sequences in nearest repeat region                       |

## Step 6 Solution

#### Feature selection 

I added the feature `repeat_seqs_w_overlap_ratio` as well as multiple palindrome-related features. The thought process for the former was to ensure that if there is overlap between an encoding region and repeat region, this is accounted for.

The thought process for the latter was a result of researching the palindromic nature of CRISPR repeat regions. While I am unsure if this approach is at all sensical, I figured it might be interesting to see if repeat regions closest to CDS sequences that encode for Cas9 have longer or more frequent palindromic sequences than others. 

Ultimately I was unable to validate the potential usefulness of these additional features, as there were no repeat regions that overlapped with any CDS sequences, and the sample size of loci with repeat_region features and non-Cas9 CDS sequences was extremely small in the given data files. Given more time, I would like to try using a larger/more diverse dataset.

#### cas9 Feature observations

The average `protein_length` for Cas9 proteins was 1210.50, whereas the overall average was 295.81. Upon first glance, it seems that Cas9 proteins have significantly larger protein lengths. Additionally, the average `min_dist_to_nearest_repeat ` was only 2112.00, as compared to the overall average of 83774.39. 

While this definitely indicates that these two criteria could be used to filter Cas9 proteins (i.e. relatively low distance to nearest repeat regions, and relatively large protein length), there were only 2 Cas9 samples in this dataset, meaning statistical significance cannot be concluded.

#### Developing a scalable solution

The dataset provided was extremely small, making it possible for me to do all feature computations quickly on my local machine. This includes using a brute-force approach to finding palindromic sub-sequences in repeat regions, by looping through a list of all possible sub-sequences. 

My first thoughts, given my lack of experience in the field of the storage and processing of gene sequence data, is to parallelize as many components as possible (e.g. using technologies like Hadoop and Spark). In this scenario, if we are able to collect sufficient labeled data (i.e. labeled Cas9 proteins) develop a mathematical algorithm or Machine Learning model to make automated predictions as to whether or not a new protein is Cas9.

We can then parallel process a large amount of records using the map-reduce paradigm (which has gotten quicker in recent years, with technologies like Spark allowing us to peform in-memory processing on clusters) to extract features from records simultaneously, feed the features through an algorithm/model simultaneously, and reduce the resulting set to entries that have a high probability of being Cas9.
