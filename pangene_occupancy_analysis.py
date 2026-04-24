"""
pangene_occupancy_analysis.py

Pan-gene properties by occupancy for Candida auris (6 genomes).
Author: Faiza Khan, MSc Bioinformatics, University of Liverpool
Supervisor: Prof. Andy Jones

Style inspired by Contreras-Moreira et al. (2026) Genome Research 36:226 to 238.

This script computes real pan-gene occupancy (1 to 6) from the GET_PANGENES
-t 0 matrix and visualises gene properties across occupancy classes:
    cloud    (occupancy 1)       species-unique genes
    shell    (occupancy 2 to 4)  accessory genes
    softcore (occupancy 5)
    core     (occupancy 6)       present in all 6 genomes

Inputs:
    1) GET_PANGENES pan-gene matrix (.xlsx or .tsv from -t 0 run)
       Rows = pan-gene clusters, cols = 6 C. auris genomes
    2) Reference gene summary TSV (B8441 annotations from VEuPathDB)

Usage (from VS Code terminal):
    python pangene_occupancy_analysis.py \
        --matrix data/cauris_get_pangenes_t0_run.xlsx \
        --summary data/Cauris_ref_genes_Summary_project.tsv \
        --outdir figures/

Dependencies: pandas, numpy, matplotlib, scipy, openpyxl
"""

import argparse
import os
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats


# Configuration

CLASS_COLORS = {
    'cloud':    '#F8A49A',
    'shell':    '#40C4BE',
    'softcore': '#C4A3D4',
    'core':     '#7AAA3E',
}

GENOME_COLUMNS = {
    'FungiDB-68_CaurisB11220': 'Cauris_B11220',
    'FungiDB-68_CaurisB11245': 'Cauris_B11245',
    'FungiDB-68_CaurisB11243': 'Cauris_B11243',
    'FungiDB-68_CaurisB8441':  'Cauris_B8441',
    'FungiDB-68_CaurisB11221': 'Cauris_B11221',
    'FungiDB-68_Cauris6684':   'Cauris_6684',
}


def classify_occupancy(occ):
    if occ == 1:
        return 'cloud'
    if occ in (2, 3, 4):
        return 'shell'
    if occ == 5:
        return 'softcore'
    if occ == 6:
        return 'core'
    return 'unknown'


# Data loading

def load_pangene_matrix(path):
    """Load GET_PANGENES -t 0 output and compute occupancy."""
    if path.lower().endswith('.xlsx'):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, sep='\t')

    df = df.rename(columns=GENOME_COLUMNS)
    new_genome_cols = list(GENOME_COLUMNS.values())
    df = df.rename(columns={df.columns[0]: 'pangene_id'})

    df = df[df[new_genome_cols].notna().any(axis=1)].reset_index(drop=True)

    def is_present(cell):
        if pd.isna(cell):
            return False
        return str(cell).strip() not in ['-', '', 'nan']

    df['occupancy'] = df[new_genome_cols].apply(
        lambda row: sum(is_present(c) for c in row), axis=1)
    df['class'] = df['occupancy'].apply(classify_occupancy)
    df = df[df['occupancy'] >= 1].reset_index(drop=True)

    print(f"Loaded pan-gene matrix: {len(df):,} clusters")
    print("\nOccupancy distribution:")
    for occ in sorted(df['occupancy'].unique()):
        n = (df['occupancy'] == occ).sum()
        cls = classify_occupancy(occ)
        print(f"  Occupancy {occ} ({cls:>8s}): {n:>5,} clusters")
    return df


def merge_with_reference(pg, summary_path):
    """Merge B8441 reference gene properties onto the pan-gene matrix."""
    ref = pd.read_csv(summary_path, sep='\t')

    for col in ['Protein Length', 'Ortholog count', 'Paralog count',
                'Molecular Weight']:
        ref[col] = pd.to_numeric(ref[col], errors='coerce')

    ref['expr_sensitive'] = pd.to_numeric(ref.iloc[:, 12], errors='coerce')
    ref['expr_resistant'] = pd.to_numeric(ref.iloc[:, 13], errors='coerce')

    ref['has_interpro'] = ref['Interpro ID'].apply(
        lambda x: False if (pd.isna(x) or str(x).strip() in ['N/A', '', 'nan'])
        else True)
    ref['is_annotated'] = ref['Product Description'].apply(
        lambda x: False if (pd.isna(x) or str(x).strip() == 'hypothetical protein')
        else True)

    merged = pg.merge(ref, how='inner',
                      left_on='Cauris_B8441', right_on='Gene ID')
    print(f"\nMerged {len(merged):,} pan-genes with B8441 reference annotations")
    return merged


# Plotting helpers

def setup_style():
    plt.rcParams.update({
        'font.family':       'DejaVu Sans',
        'font.size':         10,
        'axes.spines.top':   False,
        'axes.spines.right': False,
        'axes.linewidth':    1.0,
        'axes.labelsize':    10,
        'axes.titlesize':    11,
        'axes.titleweight':  'bold',
        'xtick.labelsize':   9,
        'ytick.labelsize':   9,
        'figure.facecolor':  'white',
        'axes.facecolor':    'white',
    })


def occupancy_colour_list(df):
    order = sorted(df['occupancy'].unique())
    return [CLASS_COLORS[classify_occupancy(o)] for o in order]


def add_pvalue(ax, groups):
    valid = [g.dropna().values for g in groups if len(g.dropna()) > 1]
    if len(valid) < 2:
        return
    _, p = stats.kruskal(*valid)
    txt = 'p : < 0.0001' if p < 0.0001 else f'p : {p:.3f}'
    y0, y1 = ax.get_ylim()
    yb = y1 * 0.93
    h = (y1 - y0) * 0.015
    ax.plot([0.6, 0.6, len(groups) + 0.4, len(groups) + 0.4],
            [yb, yb + h, yb + h, yb], color='black', lw=0.8)
    ax.text((len(groups) + 1) / 2, yb + h * 1.3, txt,
            ha='center', va='bottom', fontsize=9)
    ax.set_ylim(y0, y1 * 1.03)


def boxplot_by_occupancy(ax, df, column, ylabel, letter):
    order = sorted(df['occupancy'].unique())
    data = [df[df['occupancy'] == o][column].dropna() for o in order]
    colours = occupancy_colour_list(df)

    bp = ax.boxplot(
        data, positions=order, widths=0.65, patch_artist=True,
        medianprops=dict(color='black', linewidth=1.3),
        whiskerprops=dict(linewidth=0.9, color='black'),
        capprops=dict(linewidth=0.9, color='black'),
        flierprops=dict(marker='o', markersize=2.5, alpha=0.45,
                        markerfacecolor='#555555', markeredgecolor='none'),
    )
    for patch, c in zip(bp['boxes'], colours):
        patch.set_facecolor(c)
        patch.set_alpha(0.95)
        patch.set_edgecolor('black')
        patch.set_linewidth(0.9)

    ax.set_xticks(order)
    ax.set_xticklabels([str(o) for o in order])
    ax.set_xlabel('Pangene occupancy')
    ax.set_ylabel(ylabel)
    ax.text(-0.18, 1.04, letter, transform=ax.transAxes,
            fontsize=13, fontweight='bold')
    add_pvalue(ax, data)


def barplot_percentage(ax, df, column, ylabel, letter, ymax=105):
    order = sorted(df['occupancy'].unique())
    vals = [df[df['occupancy'] == o][column].mean() * 100 for o in order]
    colours = occupancy_colour_list(df)

    ax.bar(order, vals, width=0.7, color=colours,
           edgecolor='black', linewidth=0.9)
    ax.set_xticks(order)
    ax.set_xticklabels([str(o) for o in order])
    ax.set_xlabel('Pangene occupancy')
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, ymax)
    ax.text(-0.18, 1.04, letter, transform=ax.transAxes,
            fontsize=13, fontweight='bold')


def barplot_counts(ax, df, letter):
    order = sorted(df['occupancy'].unique())
    vals = [len(df[df['occupancy'] == o]) for o in order]
    colours = occupancy_colour_list(df)
    bars = ax.bar(order, vals, width=0.7, color=colours,
                  edgecolor='black', linewidth=0.9)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v + max(vals) * 0.01,
                f'{v:,}', ha='center', va='bottom',
                fontsize=8.5, fontweight='bold')
    ax.set_xticks(order)
    ax.set_xticklabels([str(o) for o in order])
    ax.set_xlabel('Pangene occupancy')
    ax.set_ylabel('Number of Pan-gene clusters')
    ax.text(-0.18, 1.04, letter, transform=ax.transAxes,
            fontsize=13, fontweight='bold')


# Combined figure

def make_combined_figure(pg_all, merged, outdir):
    merged['log_sensitive']        = np.log2(merged['expr_sensitive'].fillna(0) + 1)
    merged['Molecular Weight_kDa'] = merged['Molecular Weight'] / 1000

    fig = plt.figure(figsize=(11, 11))
    gs = fig.add_gridspec(3, 2, hspace=0.55, wspace=0.32,
                          left=0.09, right=0.97, top=0.93, bottom=0.10)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])
    axE = fig.add_subplot(gs[2, 0])
    axF = fig.add_subplot(gs[2, 1])

    boxplot_by_occupancy(axA, merged, 'Protein Length',
                         'Number of Amino acids', 'A')
    barplot_percentage(axB, merged, 'has_interpro',
                       'Percentage of Pan-genes\nwith an Interpro id', 'B')
    boxplot_by_occupancy(axC, merged, 'Paralog count', 'Paralogs count', 'C')
    boxplot_by_occupancy(axD, merged, 'log_sensitive',
                         'Gene expression\n(in log\u2082(TPM+1))', 'D')
    boxplot_by_occupancy(axE, merged, 'Molecular Weight_kDa',
                         'Molecular Weight (kDa)', 'E')
    barplot_counts(axF, pg_all, 'F')

    fig.suptitle('Exploration of pangenes and their properties by their occupancy\n'
                 'Candida auris: 6 genomes (GET_PANGENES, -t 0 run)',
                 fontsize=13.5, fontweight='bold', y=0.98)

    legend_patches = [
        mpatches.Patch(facecolor=CLASS_COLORS['cloud'],
                       edgecolor='black', label='cloud (1)'),
        mpatches.Patch(facecolor=CLASS_COLORS['shell'],
                       edgecolor='black', label='shell (2 to 4)'),
        mpatches.Patch(facecolor=CLASS_COLORS['softcore'],
                       edgecolor='black', label='softcore (5)'),
        mpatches.Patch(facecolor=CLASS_COLORS['core'],
                       edgecolor='black', label='core (6)'),
    ]
    fig.legend(handles=legend_patches, loc='lower center',
               ncol=4, frameon=False, fontsize=10,
               title='occupancy class', title_fontsize=10,
               bbox_to_anchor=(0.5, 0.02))
    fig.text(0.5, -0.005,
             'Style after Contreras-Moreira et al. Genome Res. 2026; 36:226 to 238',
             ha='center', fontsize=9, style='italic', color='#444')

    out = os.path.join(outdir, 'figure_pangene_occupancy_combined.png')
    plt.savefig(out, dpi=220, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\nCombined figure saved: {out}")


def save_individual_figures(pg_all, merged, outdir):
    merged['log_sensitive']        = np.log2(merged['expr_sensitive'].fillna(0) + 1)
    merged['Molecular Weight_kDa'] = merged['Molecular Weight'] / 1000

    for letter, col, ylabel in [
        ('A', 'Protein Length',       'Number of Amino acids'),
        ('C', 'Paralog count',        'Paralogs count'),
        ('D', 'log_sensitive',        'Gene expression\n(in log\u2082(TPM+1))'),
        ('E', 'Molecular Weight_kDa', 'Molecular Weight (kDa)'),
    ]:
        fig, ax = plt.subplots(figsize=(5.2, 4.2))
        boxplot_by_occupancy(ax, merged, col, ylabel, letter)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir,
                                  f'fig_{letter}_{col.replace(" ", "_")}.png'),
                    dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    barplot_percentage(ax, merged, 'has_interpro',
                       'Percentage of Pan-genes\nwith an Interpro id', 'B')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fig_B_interpro_percentage.png'),
                dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    barplot_counts(ax, pg_all, 'F')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fig_F_cluster_counts.png'),
                dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Pan-gene occupancy analysis for C. auris')
    parser.add_argument('--matrix',  required=True,
                        help='GET_PANGENES matrix file (.xlsx or .tsv)')
    parser.add_argument('--summary', required=True,
                        help='Reference gene summary TSV (B8441)')
    parser.add_argument('--outdir',  default='figures',
                        help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    setup_style()

    pg = load_pangene_matrix(args.matrix)
    merged = merge_with_reference(pg, args.summary)

    make_combined_figure(pg, merged, args.outdir)
    save_individual_figures(pg, merged, args.outdir)

    print(f"\nAll figures saved to: {args.outdir}/")


if __name__ == '__main__':
    main()
