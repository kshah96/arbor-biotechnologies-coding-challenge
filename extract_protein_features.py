from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import itertools
import click

def compute_sequence_distance_features(seq_a, *seq_list):

    """Computes distances from input seq_a to an arbitrary number of additional sequences passed

    Args:
        seq_a      (Bio.SeqFeature.FeatureLocation)        : Sequence under consideration
        *seq_list  (list of Bio.SeqFeature.FeatureLocation): List of sequences to compare with seq_a

    Returns:
        list of int: List of integers of len(seq_list) denoting distance in min-number-of-bases-away to seq_a
        list of bool: List of bool of len(seq_list) denoting whether or not a sequence overlaps with seq_a
    """

    min_distance_list = []
    is_overlapping_list = []

    # Get integer values representing start and end of seq_a
    seq_a_start = seq_a.start.real
    seq_a_end = seq_a.end.real

    # Do a simple O(n) time comparison against each sequence in seq_list (where n is len(seq_list))
    for seq in seq_list:
       seq_n_start = seq.start.real
       seq_n_end = seq.end.real

       start_to_start_dist = abs(seq_a_start - seq_n_start)
       start_to_end_dist = abs(seq_a_start - seq_n_end)
       end_to_start_dist = abs(seq_a_end - seq_n_start)
       end_to_end_dist = abs(seq_a_end - seq_n_end)

       min_distance = min(start_to_start_dist, start_to_end_dist, end_to_start_dist, end_to_end_dist)

       # NOTE: We will also indicate whether or not sequences overlap as a feature
       is_overlapping = max(seq_a_start, seq_n_start) <= min(seq_a_end, seq_n_end)

       min_distance_list.append(min_distance)
       is_overlapping_list.append(is_overlapping)

    return min_distance_list, is_overlapping_list

def compute_protein_features(record_id, CDS_features, repeat_region_features):
    """Extracts numeric features characteristic of protein from given locus feature information

    Args:
        record_id              (string)                              : Locus GenBank ID
        CDS_features           (iterable[bio.genbank.record.feature]): Locus CDS features
        repeat_region_features (iterable[bio.genbank.record.feature]): Locus repeat_region features

    Returns:
       list(dict): dicts containing numeric features extracted from locus for given CDS region
    """

    protein_features = []

    repeat_region_features = list(repeat_region_features)

    for CDS_feature in CDS_features:
        location = CDS_feature.location
        feature_type = CDS_feature.type

        qualifiers = CDS_feature.qualifiers
        protein_id = qualifiers.get('protein_id')
        gene = qualifiers.get('gene', [None])[0]
        translation = qualifiers.get('translation', [None])[0]

        # TODO: Research the reason as to why some CDS_features do not have translations
        # For now just skip any features missing translation string
        if translation is None:
            continue

        # Translation string is read as an array with 1 index
        protein_length = len(translation)

        min_distance_list, overlapping_list = compute_sequence_distance_features(location, *[rrf.location for rrf in repeat_region_features])

        if not min_distance_list:
            min_dist_to_nearest_repeat = None
            min_dist_index = None
            ratio_overlapping_sequences = None
        else:
            min_dist_to_nearest_repeat = min(min_distance_list)
            min_dist_index = min_distance_list.index(min_dist_to_nearest_repeat)
            ratio_overlapping_sequences = overlapping_list.count(True) / len(overlapping_list)

        if min_dist_index is not None:
            # Calculate statistics relating to the existence of palindromic sub-sequences within the nearest repeat region
            # NOTE: I took a brute force approach to finding ALL possible sub-sequences. This would not be feasible on
            # a large dataset and should be re-thought
            min_dist_repeat_region = repeat_region_features[min_dist_index]
            qualifiers = min_dist_repeat_region.qualifiers
            repeat_str = qualifiers.get('rpt_unit_seq')[0]

            # Compute all contiguos substrings from repeat sequence (the same substring can be repeated)
            # NOTE: Idea obtained from here: https://www.geeksforgeeks.org/python-get-all-substrings-of-given-string/
            repeat_substrings = [repeat_str[i: j] for i in range(len(repeat_str))
                    for j in range(i + 1, len(repeat_str) + 1)]

            # Find all palindromic substrings with base length at least 3
            palindrome_substrings = []
            for substring in repeat_substrings:
                # For now, rule out any base sequences of 1 or 2 as they are probably unimportant
                # TODO: Figure out a more accurate metric based on research done in this field
                if len(substring) < 3:
                    continue

                substring_seq = Seq(substring)
                substring_reverse_str = str(substring_seq.reverse_complement())
                #import pdb;pdb.set_trace()

                if substring_reverse_str == substring:
                    palindrome_substrings.append(substring)

            # TODO: Very inefficient -- compute these values while iterating over substrings initially
            # Extract some potentially meaningful features
            most_freq_palindrome = max(set(palindrome_substrings), key=palindrome_substrings.count)
            len_most_freq_palindrome = len(most_freq_palindrome)
            freq_most_freq_palindrome = len([x for x in palindrome_substrings if x == most_freq_palindrome])
            len_largest_palindrome = len(max(palindrome_substrings, key=len))
            num_unique_palindromes = len(set(palindrome_substrings))
        else:
            len_most_freq_palindrome = None
            freq_most_freq_palindrome = None
            len_largest_palindrome = None
            num_unique_palindromes = None

        # TODO: Better practice to create a nested dict to partition metadata and features,
        # but for simplicity sake creating flat dict
        protein_feature = {
            'metadata.origin_record': record_id,
            'metadata.feature_type': feature_type,
            'metadata.protein_id': protein_id[0],
            'metadata.gene': gene,
            'min_dist_to_nearest_repeat': min_dist_to_nearest_repeat,
            'repeat_seqs_w_overlap_ratio': ratio_overlapping_sequences,
            'protein_len': protein_length,
            'most_freq_palindrome_freq': freq_most_freq_palindrome,
            'most_freq_palindrome_len': len_most_freq_palindrome,
            'largest_palindrome_len': len_largest_palindrome,
            'unique_palindromes_count': num_unique_palindromes
        }

        #protein_features.append(protein_feature)
        #return protein_feature
        yield protein_feature

    #return protein_features

def generate_statistics(protein_features, protein_ids=[], genes=[]):
    """Generates statistics comparing overall average protein feature values to average protein feature values for a particular gene or protein.
    Prints metrics and also returns complete pandas dataframe containing data

    Args:
        protein_features (iterable[dict])  : Iterator over dicts containing generated protein feature information
        protein_ids      (iterable[string]): GenBank CDS protein_ids to generate statistics for (e.g. "WP_002460848.1") (optional)
        genes            (iterable[string]): GenBank CDS genes to generate statistics for (e.g. "cas9") (optional)

    Returns:
        pandas.DataFrame: DataFrame containing result metrics
    """

    # To do aggregations using pandas we will need to load all features into memory...
    # TODO: Explore better ways to do this on large datasets
    df = pd.DataFrame(list(protein_features))

    print('{:*^100}'.format('OVERALL METRICS'))
    print(df.describe().to_markdown())

    for protein_id in protein_ids:
        print('\n{:*^100}'.format('PROTEIN {} METRICS'.format(protein_id)))
        print(df[df['metadata.protein_id'] == protein_id].describe().to_markdown())

    for gene in genes:
        print('\n{:*^100}'.format('GENE {} METRICS'.format(gene)))
        print(df[df['metadata.gene'] == gene].describe().to_markdown())

# NOTE: Currently only supports "genbank" file format
def extract_protein_features_from_file(file_paths, file_format='genbank'):
    """Extracts collections of protein features from each locus in GenBank file(s)

    Args:
        file_paths   (list of str): List of paths corresponding to data files
        file_format  (str)        : Format of files to be read by biopython (should always be "genbank")

    Returns:
        iterable[dict]: Iterator over dicts containing computed protein features and CDS feature metadata
    """

    # Use itertools.chain() method to chain multiple generators
    all_records = itertools.chain()
    for f in file_paths:
        all_records = itertools.chain(all_records, SeqIO.parse(f, file_format))

    for record in all_records:

        features = record.features

        # NOTE: For the purpose of this project, we are concerned with only CDS and repeat_region features
        # Construct two generators for CDS and repeat_region features respectively
        CDS_features = filter(lambda feature: feature.type == 'CDS', features)
        repeat_region_features = filter(lambda feature: feature.type == 'repeat_region', features)

        protein_features = compute_protein_features(record.id, CDS_features, repeat_region_features)

        # Iterate over generator in order to create new generator with ALL entries
        # TODO: Find more conventional way to do this
        for feature in protein_features:
            yield(feature)

@click.command()
@click.argument('file_paths', nargs=-1, required=True)
@click.option('--gene', help='Gene (e.g. cas9) to display feature metrics for', multiple=True)
@click.option('--protein_id', help='protein_id (e.g. WP_053019794.1) to display feature metrics for', multiple=True)
def main(file_paths, gene, protein_id):
    protein_features = extract_protein_features_from_file(file_paths)
    df = generate_statistics(protein_features, genes=gene, protein_ids=protein_id)

if __name__ == '__main__':
    main()
