# Outputting .pkl files with reduced size (omitting df_clonal)
# Command for submitting array job - enter this on the command line
# 5740 genes in the gff file

"""
bsub -J "gene_sums[1-5102]%500" -o output.%J%I -e error.%J%I -M6000 -R "select[mem>6000] rusage[mem=6000] span[hosts=1]" -G team342 '/software/isg/languages/Python-3.9.10/bin/python3 generate_gene_summary.py ${LSB_JOBINDEX}'
"""

# Imports
import sys
import os
import numpy as np
import pandas as pd
import vcf
import collections
from Bio.Seq import Seq
import allel
from pyfasta import Fasta
import pickle
import lzma
import json
import time
from filelock import FileLock

# Array job info
LSB_JOBINDEX = int(sys.argv[1]) - 1

# Get list of all Pf gene ids
gff_fn = '/lustre/scratch124/gsu/legacy/pfalciparum/resources/snpEff/data/Pfalciparum_GeneDB_Feb2020/Pfalciparum_replace_Pf3D7_MIT_v3_with_Pf_M76611.gff'
df_gff = allel.gff3_to_dataframe(gff_fn, attributes=['ID', 'Name'])

# Read IDs in json as into a list
Pf3D7_core_genes_list = pd.read_json('/nfs/users/nfs_e/eu1/mutation-discovery-app/app/files/core_genes.json', orient='index').index.to_list()

# Array job setup
gene_id = Pf3D7_core_genes_list[LSB_JOBINDEX]

# Dictionary to save pipeline logs: sample inclusion/exclusion per gene
gene_logs =  dict()

# In the first gene, log some detals
if LSB_JOBINDEX==1:
    username = os.path.expanduser("~").split(os.path.sep)[-1]
    gene_logs['user'] = username
    current_time = time.strftime("%Y-%m-%d %H:%M:%S")
    gene_logs['run_time']= current_time

# Read in data files
vcf_file_format = "/lustre/scratch124/gsu/legacy/pipelines/builds/pf_70_build/pf_70_internal_release/vcf/%s.pf7.vcf.gz"
samples_fn = '/lustre/scratch126/gsu/team112/pf7/ftp_20221021/Pf7_samples.txt'
ref_genome_fn = '/lustre/scratch124/gsu/legacy/pfalciparum/resources/Pfalciparum.genome.fasta'
gff_fn = '/lustre/scratch124/gsu/legacy/pfalciparum/resources/snpEff/data/Pfalciparum_GeneDB_Feb2020/Pfalciparum_replace_Pf3D7_MIT_v3_with_Pf_M76611.gff'

# Function 1 - determine coding sequence
def determine_cds(gene_id='PF3D7_0709000', gff_fn=gff_fn):
    df_gff = allel.gff3_to_dataframe(gff_fn, attributes=['ID', 'Name'])
    gene_name = df_gff.loc[df_gff['ID'] == gene_id, 'Name'].values[0]
    if (gene_name == '') or (gene_name == '.'):
        gene_name=gene_id
    chrom = df_gff.loc[df_gff['ID'] == gene_id, 'seqid'].values[0]
    strand = df_gff.loc[df_gff['ID'] == gene_id, 'strand'].values[0]
    df_cds = df_gff.loc[
        ( df_gff['ID'].str.startswith(gene_id) )
        & ( df_gff['type'] == 'CDS' )
    ]
    return(
        {
            'gene_name':gene_name,
            'chrom':chrom,
            'strand':strand,
            'starts':df_cds['start'].values,
            'ends':df_cds['end'].values,
        }
    )

# Function 2 - call haplotypes per gene and sample

def call_haplotype(
#     target_name='crt',
    chrom='Pf3D7_07_v3',
    starts=[403222, 403490, 403938, 404283, 404569, 404764, 404936, 405146, 405334, 405539, 405825, 406017, 406241], # Start positions of each exon
    ends=  [403312, 403758, 404110, 404415, 404640, 404839, 405018, 405196, 405390, 405631, 405869, 406071, 406317], # End positions of each exon
    strand='+', # Either '+' or '-'
    vcf_file_format=vcf_file_format, # Location of VCFs, allowing chrom to vary
    samples=None, # The default (None) will call haplotypes for all samples within the VCF
    include_indels=False, # If set to True, this will include indels, and hence potentially output sequences of different lengths
    positions_with_missingness=[403620, 403621], # If there are missing genotypes at these positions, these will not be treated as missing. This is important in CRT as many samples have missing genotypes at Pf3D7_07_v3:403620-403621 due to way haplotypes are treated by HaplotypeCaller
    positions_to_ignore=[], # These positions will be skipped
    positions_to_allow=[403618, 403622], # These allows certain indels to be included even if include_indels is False. Important as common CVIET haplotype in CRT 72-76 is treated by HaplotypeCaller as two 1bp indels, rather than 3 SNPs
    ref_fasta_fn=ref_genome_fn,
    vcf_fn=None, # If supplied, this will overrule what is provided in vcf_file_format
    verbose=False,
    show_genotypes=False, # If true, these will give very verbose output
    use_filter=True, # If true, will only use PASS variants
    first_aa=1, # Important to use if, for example, only creating haplotypes for one exon
    to_stop=True # This argument is used by Bio translate function. If True, sequences will be truncated at premature stop codons, resulting in shorter AA sequences
):
    

    """ This section checks the data types of 'starts' and 'ends' 
   - Is 'starts' neither a list nor a numpy array? If not, check if its an integer. If it is an integer, convert to list. If not, exit
   - Repeat for 'ends' 
   - Are the converted lists of the same length? Exit if not
   - Finally, create a list of exon start positions, relative to the preceeding exon. This is done by calculating diffs between start and ends, 
   getting the cumsum (minus the last entry) and prefixing with 0. """

    if not isinstance(starts, list) and not(isinstance(starts, np.ndarray)):
        if isinstance(starts, int):
            starts = [starts]
        else:
            sys.exit("starts must be a list of integers (or a single integer)")
    if not isinstance(ends, list) and not(isinstance(ends, np.ndarray)):
        if isinstance(ends, int):
            ends = [ends]
        else:
            sys.exit("ends must be a list of integers (or a single integer)")
    
    if len(starts) != len(ends):
        sys.exit("starts and ends must be the same length")
        
    exon_offsets = [0] + list(np.cumsum(np.array(ends)-np.array(starts)+1))[:-1]
    
    """ If vcf_fn is not provided, generate a vcf filename based on the value of the chrom variable using the vcf_file_format template.  """
    if vcf_fn is None:
        vcf_fn = vcf_file_format % chrom
        
    """Generate refs: processing a series of exons specified by starts and ends, 
 combining them into a single sequence, and translating that sequence into an amino acid sequence, with proper handling for the strand direction. """
    ref_sequences = []
    for i in np.arange(len(starts)):
        start = starts[i]
        end = ends[i]
        ref_sequences.append(Seq(Fasta(ref_fasta_fn)[chrom][starts[i]-1:ends[i]]))
    ref_sequence = Seq("")
    for s in ref_sequences:
        ref_sequence += s
    ref_nucleotide_haplotype = str(ref_sequence)
    if strand == '+':
        ref_aa_haplotype = str(ref_sequence.translate(to_stop=to_stop))
    else:
        ref_aa_haplotype = str(ref_sequence.reverse_complement().translate(to_stop=to_stop))

    vcf_reader = vcf.Reader(filename=vcf_fn) # Read the vcf file

#""" Set up default samples = None to read in all samples from vcf. Initialise a set of dictionaries"""
    if samples is None:
        samples = vcf_reader.samples
    sample_sequences = collections.OrderedDict() # will be a set of two sequences to represent the two alt (?) calls for a sample
    sample_offsets = collections.OrderedDict()
    PID = collections.OrderedDict() # phased IDs
    non_phased_het_seen = collections.OrderedDict()
    sample_unphaseable = collections.OrderedDict()
    aa_haplotype = collections.OrderedDict()
    nucleotide_haplotype = collections.OrderedDict()
    ns_changes = collections.OrderedDict()
    num_missing = collections.OrderedDict()
    num_called = collections.OrderedDict()
    swap_phasing = collections.OrderedDict()
    first_PGT = collections.OrderedDict()
    

    for sample in samples:
        #sample_sequences[sample] = [MutableSeq(ref_sequence), MutableSeq(ref_sequence)] # For using Biopython 1.81
        sample_sequences[sample] = [ref_sequence.tomutable(), ref_sequence.tomutable()] # a list with two identical entries.
        sample_offsets[sample] = [0, 0]
        PID[sample] = ''
        non_phased_het_seen[sample] = False
        sample_unphaseable[sample] = ''
        aa_haplotype[sample] = ''
        nucleotide_haplotype[sample] = ''
        num_missing[sample] = 0
        num_called[sample] = 0
        swap_phasing[sample] = False
        first_PGT[sample] = ''

    for i in np.arange(len(starts)):
        start = starts[i]
        end = ends[i]
        exon_offset = exon_offsets[i]
        if verbose:
            print(f"{chrom}:{starts[i]}-{ends[i]}")
        for record in vcf_reader.fetch(chrom, start-1, end): # Find variants in relevant portion of gene
            if ( record.POS >= start) and not record.POS in positions_to_ignore: # Have to restrict to position that are at or after start, to ensure indels in previous positions don't mess things up. Note this is a change from previous behaviour
                if not use_filter or record.FILTER==[]: # Only use PASS variants, unless use_filter is False
                    if verbose:
                        print("\n", record, record.FILTER)
                    for sample in samples:
                        if show_genotypes:
                            print(record.genotype(sample))
                        GT = record.genotype(sample)['GT']
                        AD = record.genotype(sample)['AD']
                        POS = record.POS - start + exon_offset
                        if GT == './.' or GT == '.|.': # Note phased empty genotypes (.|.) have appeared in Pf7. These weren't in Pf6!
                            if not record.POS in positions_with_missingness:
                                num_missing[sample] += 1
                        else:
                            num_called[sample] += 1
                            if GT != '0/0' and GT != '0|0': # Only interested in non-ref genotypes. Hom ref genotypes don't change sequence so can be ignored
                                alleles = [record.REF] + record.ALT
#                                 first_allele_int = int(GT[0])
#                                 second_allele_int = int(GT[2])
                                # Find lengths of alleles (to determine which are indels)
                                REF_len = len(record.REF)
                                first_allele_ALT_len = len(alleles[int(GT[0])]) # not necessarily alts
                                second_allele_ALT_len = len(alleles[int(GT[2])])
                                if (REF_len==first_allele_ALT_len and REF_len==second_allele_ALT_len) or include_indels or record.POS in positions_to_allow:
                                    if (GT[0] != GT[2]) and (alleles[int(GT[0])] != '*') and (alleles[int(GT[2])] != '*'): # heterozygous call at non-spanning deletion
                                        if 'PGT' in record.genotype(sample).data._fields:
                                            is_unphased = (record.genotype(sample)['PGT'] is None)
                                        else:
                                            is_unphased = True
                                        if non_phased_het_seen[sample]:
                                            if verbose:
                                                print("Sample %s unphased het followed by het" % sample)
                                            sample_unphaseable[sample] = '*'
            #                                 genotype[sample] = 'X'
                                        if is_unphased:
                                            non_phased_het_seen[sample] = True
                                        else:
#                                             if 'PID' in record.genotype(sample).data._fields: # Note this won't be the case in het spanning deletions
                                            if PID[sample] == '':
                                                PID[sample] = record.genotype(sample)['PID']
                                                first_PGT[sample] = record.genotype(sample)['PGT']
                                                if AD[int(GT[2])] > AD[int(GT[0])]:
                                                    swap_phasing[sample] = True
                                            else:
                                                if record.genotype(sample)['PID'] != PID[sample]:
                                                    if verbose:
                                                        print("Sample %s two PIDs" % sample)
                                                    sample_unphaseable[sample] = '*'
            #                                         genotype[sample] = '?'
#                                             if record.genotype(sample)['PGT'] == '1|0':
#                                                 GT = GT[2] + '/' + GT[0]
                                        if verbose:
                                            print(f"AD={AD}, GT={GT}")
                                        # Here we use phased genotype (PGT) if in same phasing group as first het
                                        if (
                                            'PID' in record.genotype(sample).data._fields
                                            and ( PID[sample] == record.genotype(sample)['PID'] )
                                        ): # In phasing group with first het call
                                            if verbose:
                                                print(f"PID: {PID[sample]}, GT: {GT}, PGT: {record.genotype(sample)['PGT']}")
                                            if first_PGT[sample] == record.genotype(sample)['PGT']:
                                                if swap_phasing[sample]:
                                                    GT = GT[2] + '/' + GT[0]
                                            else:
                                                if not swap_phasing[sample]:
                                                    GT = GT[2] + '/' + GT[0]
#                                             GT = record.genotype(sample)['PGT']
                                        # Else we ensure we are taking the majority haplotype as the first allele
                                        elif AD[int(GT[2])] > AD[int(GT[0])]:
                                            GT = GT[2] + '/' + GT[0]
                                            if verbose:
                                                print(f"AD={AD}, GT={GT}")

                                for i, sample_offset in enumerate(sample_offsets[sample]):
                                    alleles = [record.REF] + record.ALT
                                    GTint = int(GT[i*2])
                                    REF_len = len(record.REF)
                                    ALT_len = len(alleles[GTint])
                                    if(verbose):
                                        print(f"before: i={i}, offset={sample_offset}, alleles, GT={GTint}, REF={REF_len}, ALT={ALT_len}, {sample_sequences[sample][i]}")
                                    if GTint != 0 and alleles[GTint] != '*' and (REF_len==ALT_len or include_indels or record.POS in positions_to_allow):
                                        sample_sequences[sample][i][POS+sample_offsets[sample][i]:(POS+sample_offsets[sample][i]+REF_len)] = alleles[GTint]
                                        sample_offsets[sample][i] = sample_offset + (len(alleles[GTint]) - len(record.REF))
                                        if(verbose):
                                            print(f"after : i={i}, offset={sample_offset}, alleles, GT={GTint}, REF={REF_len}, ALT={ALT_len}, {sample_sequences[sample][i]}")

    for sample in samples:
        if(show_genotypes):
            print(sample_sequences[sample])
        ###########################
        # populate ns_changes
        ###########################
        if strand == '+':
            aa_0 = (
                str(sample_sequences[sample][0].toseq().translate(to_stop=to_stop))
                if len(sample_sequences[sample][0]) % 3 == 0
                else '!'
            )
            aa_1 = (
                str(sample_sequences[sample][1].toseq().translate(to_stop=to_stop))
                if len(sample_sequences[sample][1]) % 3 == 0
                else '!'
            )
        else:
            aa_0 = (
                str(sample_sequences[sample][0].toseq().reverse_complement().translate(to_stop=to_stop))
                if len(sample_sequences[sample][0]) % 3 == 0
                else '!'
            )
            aa_1 = (
                str(sample_sequences[sample][1].toseq().reverse_complement().translate(to_stop=to_stop))
                if len(sample_sequences[sample][1]) % 3 == 0
                else '!'
            )
#         if aa_0 != '!' and aa_1 != '!' # This was what was used in Pf6
        if aa_0 != '!' and aa_1 != '!' and (len(aa_0) == len(ref_aa_haplotype)) and (len(aa_1) == len(ref_aa_haplotype)): # In Pf7, we have some in-frame indels, hence we also need to check for these, and give "!" if indels exist
            if aa_0 == aa_1:
                ns_changes_list = []
                for i in range(len(ref_aa_haplotype)):
                    if aa_0[i] != ref_aa_haplotype[i]:
                        ns_changes_list.append("%s%d%s" % (ref_aa_haplotype[i], i+first_aa, aa_0[i]))
                ns_changes[sample] = "/".join(ns_changes_list)
            else:
                ns_changes_list_0 = []
                for i in range(len(ref_aa_haplotype)):
                    if aa_0[i] != ref_aa_haplotype[i]:
                        ns_changes_list_0.append("%s%d%s" % (ref_aa_haplotype[i], i+first_aa, aa_0[i]))
                ns_changes_0 = "/".join(ns_changes_list_0)
                ns_changes_list_1 = []
                for i in range(len(ref_aa_haplotype)):
#                     print(len(ref_aa_haplotype), len(aa_1))
                    if aa_1[i] != ref_aa_haplotype[i]:
                        ns_changes_list_1.append("%s%d%s" % (ref_aa_haplotype[i], i+first_aa, aa_1[i]))
                ns_changes_1 = "/".join(ns_changes_list_1)

                if ns_changes_0 == '':
                    ns_changes[sample] = ns_changes_1.lower()
                elif ns_changes_1 == '':
                    ns_changes[sample] = ns_changes_0.lower()
                else:
                    ns_changes[sample] = (",".join([ns_changes_0, ns_changes_1])).lower()
        else:
            if verbose:
                print(aa_0, aa_1, ref_aa_haplotype)
            ns_changes[sample] = '!'

        ################################################
        # populate nucleotide_haplotype and aa_haplotype
        ################################################
        if str(sample_sequences[sample][0]).upper() == str(sample_sequences[sample][1]).upper():
            nucleotide_haplotype[sample] = str(sample_sequences[sample][0].toseq())
            if strand == '+':
                aa_haplotype[sample] = ( str(sample_sequences[sample][0].toseq().translate(to_stop=to_stop))
                                    if len(sample_sequences[sample][0]) % 3 == 0
                                    else '!' )
            else:
                aa_haplotype[sample] = ( str(sample_sequences[sample][0].toseq().reverse_complement().translate(to_stop=to_stop))
                                    if len(sample_sequences[sample][0]) % 3 == 0
                                    else '!' )
        else:
            nucleotide_haplotype[sample] = "%s,%s" % ( str(sample_sequences[sample][0].toseq()), str(sample_sequences[sample][1].toseq()) )
            if strand == '+':
                aa_haplotype[sample] = "%s,%s" % (
                    ( str(sample_sequences[sample][0].toseq().translate(to_stop=to_stop))
                     if len(sample_sequences[sample][0]) % 3 == 0
                     else '!' ),
                    ( str(sample_sequences[sample][1].toseq().translate(to_stop=to_stop))
                     if len(sample_sequences[sample][1]) % 3 == 0
                     else '!' )
                )
            else:
                aa_haplotype[sample] = "%s,%s" % (
                    ( str(sample_sequences[sample][0].toseq().reverse_complement().translate(to_stop=to_stop))
                     if len(sample_sequences[sample][0]) % 3 == 0
                     else '!' ),
                    ( str(sample_sequences[sample][1].toseq().reverse_complement().translate(to_stop=to_stop))
                     if len(sample_sequences[sample][1]) % 3 == 0
                     else '!' )
                )
    
        # Add asterisk to end if sample is unphaseable
        ns_changes[sample] += sample_unphaseable[sample]
        aa_haplotype[sample] += sample_unphaseable[sample]
        nucleotide_haplotype[sample] += sample_unphaseable[sample]

        # Set ns_changes output to missing ('-') if no NS changes and at least one missing genotype call in the region
        if (ns_changes[sample] == '' or ns_changes[sample] == '*') and num_missing[sample] > 0:
            ns_changes[sample] = '-'
                
        # Set sequence output to missing ('-') if any missing genotype call in the gene
        if num_missing[sample] > 0:
            aa_haplotype[sample] = '-'
            nucleotide_haplotype[sample] = '-'
    
    df_all_haplotypes = pd.DataFrame(
        {
            'aa_haplotype': pd.Series(aa_haplotype),
            'nucleotide_haplotype': pd.Series(nucleotide_haplotype),
            'ns_changes': pd.Series(ns_changes),
        }
    )


#     return(df_haplotypes, ref_aa_haplotype, ref_nucleotide_haplotype)
    return(
        {
            'df_all_haplotypes':df_all_haplotypes,
            'ref_aa_haplotype':ref_aa_haplotype,
            'ref_nucleotide_haplotype':ref_nucleotide_haplotype,
        }
    )

# Function 3 - prepare data summaries for plots using both the previous functions and sample metadata
def prepare_plot_data(
    gene_id='PF3D7_0709000',
    haplotype_calls=None,
    samples_fn=samples_fn,
    expected_aa_len=None,
    populations = ['SA', 'AF-W', 'AF-C', 'AF-NE', 'AF-E', 'AS-S-E', 'AS-S-FE', 'AS-SE-W', 'AS-SE-E', 'OC-NG', ],
    # populations = ['AF-W', 'AF-C', 'AF-NE', 'AF-E', ],
    samples_to_show = {
        'PG0563-C': '7G8',
        'PG0564-C': 'GB4',
        'PG0565-C': 'HB3',
        'PG0566-C': 'Dd2',
        'PG0567-C': '3D7',
        'PG0568-C': 'IT',
    },
    sample_to_use_as_background = None, #'PG0566-C', # Dd2
    verbose=True,
    **kwargs
):
    # Haplotpye call only for QC pass samples
    df_samples = pd.read_csv(samples_fn, sep='\t', low_memory=False, index_col=0)
    qc_pass_samples = list(df_samples[df_samples['QC pass'] == True].index)
    # Add lab strain samples in haplotype calling
    ref_samples = list(df_samples[df_samples['Study'] == '1153-PF-Pf3KLAB-KWIATKOWSKI'].index)
#     print(gene_id)
    if gene_id is not None:
        cds_results = determine_cds(gene_id)
        gene_logs[gene_id] = {} 
        gene_start_time = time.time()
#     print(cds_results['gene_name'])
    
        if haplotype_calls is None:
            haplotype_calls = call_haplotype(
                samples=qc_pass_samples + ref_samples,
                chrom=cds_results['chrom'],
                starts=cds_results['starts'],
                ends=cds_results['ends'],
                strand=cds_results['strand'],
                **kwargs
            )
    else:
        haplotype_calls = call_haplotype(
            **kwargs
        )
        
    
    
    if verbose:
        print(f"Read in {df_samples.shape[0]:,} samples, of which {np.count_nonzero(df_samples['QC pass']):,} passed QC")

    df_join = df_samples.join(haplotype_calls['df_all_haplotypes'])
    
    # revert the nan values to '' resulted from the join operation 
    df_join[['aa_haplotype', 'nucleotide_haplotype', 'ns_changes']] = df_join[['aa_haplotype', 'nucleotide_haplotype', 'ns_changes']].fillna('')
    
    if expected_aa_len is None:
        expected_aa_len = len(haplotype_calls['ref_aa_haplotype'])
    
    # check you dont have two haplotypes. Clonal sample aa haplotype will == reference. Mixed will have two+ strings separated by comma.
    df_clonal = df_join.loc[
        (
            ( df_join['QC pass'] == True )
            & ( df_join['aa_haplotype'].astype(str).apply(len) == expected_aa_len )
        ) | ( df_join['Study'] == '1153-PF-Pf3KLAB-KWIATKOWSKI' ),
        [
            'aa_haplotype',
            'nucleotide_haplotype',
            'ns_changes',
            'Population',
        ]
    ]
    # check how many samples excluded/included
    # add number of reference strains which has not been removed even though QC fail
    count_qc_pass = np.count_nonzero(df_samples['QC pass'])
    count_excluded_samples = (count_qc_pass - len(df_clonal) + len(df_join[df_join['Study'] == '1153-PF-Pf3KLAB-KWIATKOWSKI']))
    gene_logs[gene_id]['c_exc_s'] = count_excluded_samples
    gene_logs[gene_id]['c_inc_s'] = df_clonal.shape[0]

    if df_clonal.shape[0] == 0:
        most_common_length = df_join['aa_haplotype'].astype(str).apply(len).value_counts().index[0]
        raise Exception(f"There were no samples with AA haplotype length = {expected_aa_len}. Most common length was {most_common_length}")
    
    df_clonal['ns_changes'].fillna('', inplace=True)
    df_clonal['number_of_mutations'] = df_clonal['ns_changes'].apply(lambda x: len(x.split('/')))
    df_clonal.loc[df_clonal['ns_changes'] == '', 'number_of_mutations'] = 0
    if verbose:
        print(f"Found {df_clonal.shape[0]:,} clonal haplotypes in QC pass samples")

    grouping_columns = [
        'number_of_mutations', 'ns_changes',
    #     'haplotype_name',
    ]
    ascending=[True]*len(grouping_columns)
    ascending.append(False)

    def n_agg(x):
        names = collections.OrderedDict()
        for population in populations:
            names[population] = np.count_nonzero( (x['Population'] == population) )
        return pd.Series(names)

    df_haplotypes = df_clonal.groupby(grouping_columns).apply(n_agg)
    df_haplotypes['Total'] = df_haplotypes.sum(axis=1)
    df_haplotypes.reset_index(inplace=True)
    df_haplotypes.sort_values(grouping_columns + ['Total'], ascending=ascending, inplace=True)
    df_haplotypes['ns_changes_list'] = df_haplotypes['ns_changes'].apply(lambda x: x.split('/'))
    df_haplotypes.sort_values('Total', ascending=False, inplace=True)
    # Determine which samples to show have which haplotypes
    df_haplotypes['sample_names'] = ''
    for sample in samples_to_show.keys():
        ns_changes = df_clonal.loc[sample, 'ns_changes']
        if df_haplotypes.loc[df_haplotypes['ns_changes']==ns_changes, 'sample_names'].values == '':
            df_haplotypes.loc[df_haplotypes['ns_changes']==ns_changes, 'sample_names'] = samples_to_show[sample]
        else:
            df_haplotypes.loc[df_haplotypes['ns_changes']==ns_changes, 'sample_names'] = (
                df_haplotypes.loc[df_haplotypes['ns_changes']==ns_changes, 'sample_names'] + '\n' + samples_to_show[sample]
            )
    if verbose:
        print(f"Found {df_haplotypes[df_haplotypes['Total']>0].shape[0]:,} unique haplotypes in QC pass samples")
        
      
    # Store integers of df_haplotype as np.unit16
    df_haplotypes = df_haplotypes.apply(lambda col: col.astype(np.uint16) if col.dtype == 'int64' else col)
    
    # Breakdown of removed QC pass samples                          
    # Initialize sets to store sample IDs for each rule
    a = set() #removed_missing_genotype
    b = set() #removed_multiple_haplotype_calls
    c = set() #removed_frameshift

    # Iterate over QC pass DataFrame rows
    for index, row in df_join.loc[df_join['QC pass']].iterrows():
        # Check and add sample IDs to the set for each rule
        if '-' in row['nucleotide_haplotype']:
            a.add(row.name)
            df_join.at[index, 'Exclusion reason'] = 'Missing_genotype'
            continue
        elif ',' in row['aa_haplotype']:
            df_join.at[index, 'Exclusion reason'] = 'Het_calls'
            b.add(row.name)
            continue
        elif '!' in row['ns_changes']:
            df_join.at[index, 'Exclusion reason'] = 'Stop_codon'
            c.add(row.name)
            continue
            
    removed_missing = len(a)
    removed_het_calls = len(b)
    removed_stop_codon = len(c)
    
    gene_logs[gene_id]['c_missing'] = removed_missing
    gene_logs[gene_id]['c_het_calls'] = removed_het_calls
    gene_logs[gene_id]['c_stop_codon'] = removed_stop_codon
    
    if verbose:
        print(f"No. of samples removed due to missing genotypes: {removed_missing:,}")
        print(f"No. of samples removed due to having het calls: {removed_het_calls:,}")
        print(f"No. of samples removed due to having a stop codon: {removed_stop_codon:,}")
    
    # Drop columns and index other than 'ns_changes'
    
    df_join = df_join.reset_index().drop(columns=df_join.reset_index().columns.difference(['ns_changes', 'Exclusion reason']))

    # Don't count lab strains as unique haplotypes if there is any, filter by 'Total' >0 since they dont have a population frequency     
    gene_logs[gene_id]['c_unq_h'] = df_haplotypes[df_haplotypes['Total']>0].shape[0]  
    
    # Determine background mutations
    if sample_to_use_as_background in df_clonal.index:
        background_ns_changes = df_clonal.loc[sample_to_use_as_background, 'ns_changes']
    else:
        background_ns_changes = ''
        
    if gene_id is not None:
        gene_name = cds_results['gene_name']
    else:
        gene_name = None
        
    # Save pickle file
    prep_plot_pickle_fn = f'{gene_id}.pkl.xz'
    with lzma.open(prep_plot_pickle_fn, 'wb') as file:
        plot_data = (df_haplotypes, df_join, background_ns_changes, gene_name) # The actual, reduced output to .pkl file
        pickle.dump(plot_data, file)
    
    if verbose:
        print(f"Finished {gene_id}, {gene_name}")

    return(
        {
            'df_haplotypes':df_haplotypes,
            'df_clonal':df_clonal,
            'df_join':df_join,
            'background_ns_changes':background_ns_changes, 
            'haplotype_calls':haplotype_calls, 
            'gene_name':gene_name
        }
    )


# Main code to run for each gene
print(f'running job for gene {gene_id}')
prepare_plot_data(gene_id=gene_id)

# Write gene_log dictionary to a JSON file
# Try to read existing logs from the JSON file

json_file_path = '/nfs/users/nfs_e/eu1/mutation-discovery-app/work/backend/gene_logs.json'
        
lock_path = json_file_path + ".lock"
lock = FileLock(lock_path)

with lock:
    try:
        with open(json_file_path, 'r') as json_file:
            log_file = json.load(json_file)
    except (FileNotFoundError, json.decoder.JSONDecodeError):
        # Create a new JSON file if not found or cannot be decoded
        with open(json_file_path, 'w') as json_file:
            json.dump(gene_logs, json_file, indent=4)
    else:
        # Update existing JSON file with new logs without removing previous data
        log_file.update(gene_logs)
        try:
            with open(json_file_path, 'w') as json_file:
                json.dump(log_file, json_file, indent=4)
        except Exception as e:
            print(f"Error writing to JSON file: {e}")