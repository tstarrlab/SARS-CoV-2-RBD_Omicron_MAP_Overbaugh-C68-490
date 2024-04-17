"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_rulegraph,
            make_summary

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------



# Information on samples and barcode runs for Omicron_BA2 ----------------------
barcode_runs_Omicron_BA2 = pd.read_csv(config['barcode_runs_Omicron_BA2'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs_Omicron_BA2.groupby(['library', 'sample'])) == len(barcode_runs_Omicron_BA2)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *sort_bin*.
sample_vs_expect_Omicron_BA2 = (
    barcode_runs_Omicron_BA2
    .query('experiment_type=="ab_selection"')
    .assign(concentration=lambda x: x['concentration'].astype(int),
            expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'sort_bin']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect_Omicron_BA2['equal'].all(), sample_vs_expect_Omicron_BA2.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_Omicron_BA2_expandR1 = (
    barcode_runs_Omicron_BA2
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )

assert barcode_runs_Omicron_BA2_expandR1['sample_lib'].nunique() == len(barcode_runs_Omicron_BA2_expandR1)
if any(barcode_runs_Omicron_BA2_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_Omicron_BA2_expandR1.query('n_R1 < 1')}")


# Rules -----------------------------------------------------------------------

# this is the target rule (in place of `all`) since it first rule listed
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        rulegraph=os.path.join(config['summary_dir'], 'rulegraph.svg'),
        env='environment_pinned.yml',
        get_mut_bind_expr=config['mut_bind_expr'],
        get_early2020_mut_bind_expr=config['early2020_mut_bind_expr'],
        bind_expr_filters_Omicron_BA2=nb_markdown('bind_expr_filters_Omicron_BA2.ipynb'),
        aggregate_variant_counts_Omicron_BA2=nb_markdown('aggregate_variant_counts_Omicron_BA2.ipynb'),
        variant_counts_Omicron_BA2=config['variant_counts_Omicron_BA2'],
        counts_to_cells_ratio_Omicron_BA2=nb_markdown('counts_to_cells_ratio_Omicron_BA2.ipynb'),
        counts_to_cells_csv_Omicron_BA2=config['counts_to_cells_csv_Omicron_BA2'],
        counts_to_scores_Omicron_BA2=nb_markdown('counts_to_scores_Omicron_BA2.ipynb'),
        escape_scores_Omicron_BA2=config['escape_scores_Omicron_BA2'],
        escape_fracs_Omicron_BA2=config['escape_fracs_Omicron_BA2'],
        call_strong_escape_sites_Omicron_BA2=nb_markdown('call_strong_escape_sites_Omicron_BA2.ipynb'),
        strong_escape_sites_Omicron_BA2=config['strong_escape_sites_Omicron_BA2'],
        escape_profiles_Omicron_BA2=nb_markdown('escape_profiles_Omicron_BA2.ipynb'),
        output_pdbs_Omicron_BA2=nb_markdown('output_pdbs_Omicron_BA2.ipynb'),
        make_supp_data_Omicron_BA2=nb_markdown('make_supp_data_Omicron_BA2.ipynb'),
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the rule graph of the computational workflow:
            ![{path(input.rulegraph)}]({path(input.rulegraph)})

            Here is the Markdown output of each notebook in the workflow:
            
            1. Get prior RBD DMS mutation-level binding and expression measurements and barcode-variant lookup table from the [SARS-CoV-2-RBD_DMS_Omicron repository](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron) and the original DMS library for SARS-CoV-2 (PCR-based mutagenesis) [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS). 

            2. Count variants and then aggregate counts for
               [Omicron_BA2]({path(input.aggregate_variant_counts_Omicron_BA2)})
               to create variant counts files for
               [Omicron_BA2]({path(input.variant_counts_Omicron_BA2)}).

            3. Analyze sequencing counts to cells ratio for 
               [Omicron_BA2]({path(input.counts_to_cells_ratio_Omicron_BA2)})
               this prints a list of any samples where this ratio too low. Also
               creates a CSV for 
               [Omicron_BA2]({path(input.counts_to_cells_csv_Omicron_BA2)}) with the
               sequencing counts, number of sorted cells, and ratios for
               all samples.

            4. Calculate escape scores from variant counts for 
               [Omicron_BA2]({path(input.counts_to_scores_Omicron_BA2)}).

            5. Call sites of strong escape for
               [Omicron_BA2]({path(input.call_strong_escape_sites_Omicron_BA2)}).

            6. Plot escape profiles for 
               [Omicron_BA2]({path(input.escape_profiles_Omicron_BA2)}).

            7. Map escape profiles to ``*.pdb`` files using notebooks here for 
               [Omicron_BA2]({path(input.output_pdbs_Omicron_BA2)}).

            8. Make supplementary data files for 
               [Omicron_BA2]({path(input.make_supp_data_Omicron_BA2)}),
               which are here for
               [Omicron_BA2]({path(config['supp_data_dir_Omicron_BA2'])}). These include
               `dms-view` input files.


            """
            ).strip())


rule make_rulegraph:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'rulegraph.svg')
    shell:
        "snakemake --forceall --rulegraph | dot -Tsvg > {output}"

rule save_pinned_env:
    input:
        dag=os.path.join(config['summary_dir'], 'rulegraph.svg'),
    log:
    	"environment_pinned.yml"
    shell:
        """
        conda env export > {log}
        """


rule make_supp_data_Omicron_BA2:
    input:
        config['escape_profiles_config_Omicron_BA2'],
        config['output_pdbs_config'],
        config['escape_fracs_Omicron_BA2'],
        config['escape_profiles_dms_colors_Omicron_BA2']
    output:
        nb_markdown=nb_markdown('make_supp_data_Omicron_BA2.ipynb'),
        outdir=directory(config['supp_data_dir_Omicron_BA2']),
    params:
        nb='make_supp_data_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule output_pdbs_Omicron_BA2:
    input:
        config['escape_fracs_Omicron_BA2'],
        config['output_pdbs_config'],
    output:
        nb_markdown=nb_markdown('output_pdbs_Omicron_BA2.ipynb'),
        outdir=directory(config['pdb_outputs_dir_Omicron_BA2']),
    params:
        nb='output_pdbs_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule escape_profiles_Omicron_BA2:
    """Make stacked logo plots of antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs_Omicron_BA2'],
        escape_profiles_config=config['escape_profiles_config_Omicron_BA2'],
        site_color_schemes=config['site_color_schemes'],
        wildtype_sequence=config['wildtype_sequence_Omicron_BA2'],
        mut_bind_expr=config['mut_bind_expr'],
        strong_escape_sites=config['strong_escape_sites_Omicron_BA2'],
    output:
        nb_markdown=nb_markdown('escape_profiles_Omicron_BA2.ipynb'),
        escape_profiles_dms_colors=config['escape_profiles_dms_colors_Omicron_BA2'],
    params:
        nb='escape_profiles_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule call_strong_escape_sites_Omicron_BA2:
    """Call sites of strong escape."""
    input:
        escape_fracs=config['escape_fracs_Omicron_BA2'],
    output:
        nb_markdown=nb_markdown('call_strong_escape_sites_Omicron_BA2.ipynb'),
        strong_escape_sites=config['strong_escape_sites_Omicron_BA2'],
    params:
        nb='call_strong_escape_sites_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"


rule counts_to_scores_Omicron_BA2:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts_Omicron_BA2'],
        config['wildtype_sequence_Omicron_BA2'],
        # config['mut_bind_expr'],
        # config['variant_expr'],
        # config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores_Omicron_BA2.ipynb'),
        escape_scores=config['escape_scores_Omicron_BA2'],
        escape_fracs=config['escape_fracs_Omicron_BA2'],
        escape_score_samples=config['escape_score_samples_Omicron_BA2'],
    params:
        nb='counts_to_scores_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

        
rule counts_to_cells_ratio_Omicron_BA2:
    input:
        config['variant_counts_Omicron_BA2'],
        config['barcode_runs_Omicron_BA2'],
        config['wildtype_sequence_Omicron_BA2'],
    output:
        nb_markdown=nb_markdown('counts_to_cells_ratio_Omicron_BA2.ipynb'),
        counts_to_cells_csv_Omicron_BA2=config['counts_to_cells_csv_Omicron_BA2'],
    params:
        nb='counts_to_cells_ratio_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"


rule aggregate_variant_counts_Omicron_BA2:
    input:
        counts=expand(os.path.join(config['counts_dir_Omicron_BA2'],
                                   "{sample_lib}_counts.csv"),
                      sample_lib=barcode_runs_Omicron_BA2_expandR1['sample_lib']),
        fates=expand(os.path.join(config['counts_dir_Omicron_BA2'],
                                  "{sample_lib}_fates.csv"),
                     sample_lib=barcode_runs_Omicron_BA2_expandR1['sample_lib']),
        variant_table=config['bc_variant_lookup_Omicron_BA2'],
        wt_seq=config['wildtype_sequence_Omicron_BA2'],
        barcode_runs=config['barcode_runs_Omicron_BA2'],
    output:
        config['variant_counts_Omicron_BA2'],
        nb_markdown=nb_markdown('aggregate_variant_counts_Omicron_BA2.ipynb')
    params:
        nb='aggregate_variant_counts_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"


        
rule count_variants_Omicron_BA2:
    """Count variants for a specific sample."""
    input:
        variant_table=config['bc_variant_lookup_Omicron_BA2'],
        wt_seq=config['wildtype_sequence_Omicron_BA2'],
        r1s=lambda wildcards: (barcode_runs_Omicron_BA2_expandR1
                               .set_index('sample_lib')
                               .at[wildcards.sample_lib, 'R1']
                               ),
    output:
        counts=os.path.join(config['counts_dir_Omicron_BA2'], "{sample_lib}_counts.csv"),
        fates=os.path.join(config['counts_dir_Omicron_BA2'], "{sample_lib}_fates.csv"),
    params:
        sample_lib="{sample_lib}"
    run:
        # parse sample and library from `sample_lib` wildcard
        lib = params.sample_lib.split('_')[-1]
        sample = params.sample_lib[: -len(lib) - 1]
        assert sample == (barcode_runs_Omicron_BA2_expandR1
                          .set_index('sample_lib')
                          .at[params.sample_lib, 'sample']
                          )
        assert lib == (barcode_runs_Omicron_BA2_expandR1
                       .set_index('sample_lib')
                       .at[params.sample_lib, 'library']
                       )
        # initialize `CodonVariantTable` (used to get valid barcodes)
        wt_seqrecord = Bio.SeqIO.read(input.wt_seq, 'fasta')
        geneseq = str(wt_seqrecord.seq)
        primary_target = wt_seqrecord.name
        variants=dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=input.variant_table,
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target)
        # initialize `IlluminaBarcodeParser`
        parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **config['illumina_barcode_parser_params'])
        # parse barcodes
        counts, fates = parser.parse(input.r1s,
                                     add_cols={'library': lib,
                                               'sample': sample})
        # write files
        counts.to_csv(output.counts, index=False)
        fates.to_csv(output.fates, index=False)
        
rule bind_expr_filters_Omicron_BA2:
    """QC checks on bind & expression filters from DMS data.
    """
    input:
    	file=config['mut_bind_expr']
    output:
        nb_markdown=nb_markdown('bind_expr_filters_Omicron_BA2.ipynb')
    params:
        nb='bind_expr_filters_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"


rule get_mut_bind_expr:
    """Download SARS-CoV-2 mutation ACE2-binding and expression from URL for Wuhan-Hu-1, BA1, and BA2 site-saturation mutagenesis library."""
    output:
        file=config['mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['mut_bind_expr_url'], output.file)
        
rule get_early2020_mut_bind_expr:
    """Download SARS-CoV-2 Wuhan-1 mutation ACE2-binding and expression from URL."""
    output:
        file=config['early2020_mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['early2020_mut_bind_expr_url'], output.file)
