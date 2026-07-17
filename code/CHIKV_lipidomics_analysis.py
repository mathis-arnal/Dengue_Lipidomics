"""
CHIKV Lipidomics Analysis — Lipid Family Level
Generates 7 figures (PDF + JPEG) from protein-normalized lipidomics data.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
import openpyxl
import time

# ── Constants ────────────────────────────────────────────────────────────────
DATA_FILE = 'data/CHIKV-Treatment_normalized_proteins.xlsx'
OUTPUT_DIR = 'analysis/CHIKV'
ND_THRESHOLD = 1e-10

FAMILIES = ['DG', 'CER', 'CHOL', 'HEX2CER', 'HEXCER', 'LPC', 'LPE',
            'PC', 'PE', 'SM', 'TG', 'FA', 'PS', 'PI', 'PG']

COMPARISONS = [
    # 1 vs CTRL1: CHIKV effect in untreated Vero
    ('Vero DMSO CHIKV',        'Vero DMSO',           'CHIKV effect\n(Vero)'),
    # 2 vs 1 and 2 vs 7: HMN214
    ('Vero HMN214 CHIKV',      'Vero DMSO CHIKV',     'HMN214\n(drug effect)'),
    ('Vero HMN214 CHIKV',      'Vero HMN214',         'HMN214\n(CHIKV effect)'),
    # 3 vs 1 and 3 vs 8: Simvastatin
    ('Vero Simvastatin CHIKV', 'Vero DMSO CHIKV',     'Simvastatin\n(drug effect)'),
    ('Vero Simvastatin CHIKV', 'Vero Simvastatin',    'Simvastatin\n(CHIKV effect)'),
    # 4 vs 1 and 4 vs 9: Atorvastatin
    ('Vero Atorvastatin CHIKV','Vero DMSO CHIKV',     'Atorvastatin\n(drug effect)'),
    ('Vero Atorvastatin CHIKV','Vero Atorvastatin',   'Atorvastatin\n(CHIKV effect)'),
    # 5 vs 1 and 5 vs 10: Fenretinide
    ('Vero Fenretinide CHIKV', 'Vero DMSO CHIKV',     'Fenretinide\n(drug effect)'),
    ('Vero Fenretinide CHIKV', 'Vero Fenretinide',    'Fenretinide\n(CHIKV effect)'),
    # 6 vs 1 and 6 vs 11: Proflavine
    ('Vero Proflavine CHIKV',  'Vero DMSO CHIKV',     'Proflavine\n(drug effect)'),
    ('Vero Proflavine CHIKV',  'Vero Proflavine',     'Proflavine\n(CHIKV effect)'),
    # 7–11 vs CTRL1: drug effect in mock (uninfected) Vero
    ('Vero HMN214',            'Vero DMSO',           'HMN214\n(mock)'),
    ('Vero Simvastatin',       'Vero DMSO',           'Simvastatin\n(mock)'),
    ('Vero Atorvastatin',      'Vero DMSO',           'Atorvastatin\n(mock)'),
    ('Vero Fenretinide',       'Vero DMSO',           'Fenretinide\n(mock)'),
    ('Vero Proflavine',        'Vero DMSO',           'Proflavine\n(mock)'),
    # C636 — all vs C636 DMSO mock, plus drug effect (Simvastatin+CHIKV vs CHIKV)
    ('C636 CHIKV',             'C636 DMSO',           'CHIKV effect\n(C636)'),
    ('C636 Simvastatin+CHIKV', 'C636 DMSO',           'Simvastatin\n(C636 inf.)'),
    ('C636 Simvastatin+CHIKV', 'C636 CHIKV',          'Simvastatin\n(C636 drug eff.)'),
    ('C636 Simvastatin',       'C636 DMSO',           'Simvastatin\n(C636 mock)'),
]

DRUGS_VERO = ['HMN214', 'Simvastatin', 'Atorvastatin', 'Fenretinide', 'Proflavine']

FAMILY_COLORS = {
    'DG': '#e6194b', 'CER': '#3cb44b', 'CHOL': '#ffe119', 'HEX2CER': '#4363d8',
    'HEXCER': '#f58231', 'LPC': '#911eb4', 'LPE': '#42d4f4', 'PC': '#f032e6',
    'PE': '#bfef45', 'SM': '#fabed4', 'TG': '#469990', 'FA': '#dcbeff',
    'PS': '#9A6324', 'PI': '#fffac8', 'PG': '#800000'
}

# ── Species-level comparisons ─────────────────────────────────────────────────
# Each entry: ((treatment, control), slug_for_filename, human_readable_label)
SPECIES_COMPARISONS_VERO = [
    (('Vero DMSO CHIKV',        'Vero DMSO'),         'CHIKV+DMSO_vs_DMSO',
     'CHIKV + DMSO vs. DMSO (Vero)'),
    (('Vero Simvastatin',        'Vero DMSO'),         'Simvastatin_vs_DMSO',
     'Simvastatin vs. DMSO (Vero mock)'),
    (('Vero Simvastatin CHIKV',  'Vero DMSO CHIKV'),  'CHIKV+Simvastatin_vs_CHIKV+DMSO',
     'CHIKV + Simvastatin vs. CHIKV + DMSO (Vero)'),
    (('Vero Fenretinide',        'Vero DMSO'),         'Fenretinide_vs_DMSO',
     'Fenretinide vs. DMSO (Vero mock)'),
    (('Vero Fenretinide CHIKV',  'Vero DMSO CHIKV'),  'CHIKV+Fenretinide_vs_CHIKV+DMSO',
     'CHIKV + Fenretinide vs. CHIKV + DMSO (Vero)'),
    (('Vero Proflavine',         'Vero DMSO'),         'Proflavine_vs_DMSO',
     'Proflavine vs. DMSO (Vero mock)'),
    (('Vero Proflavine CHIKV',   'Vero DMSO CHIKV'),  'CHIKV+Proflavine_vs_CHIKV+DMSO',
     'CHIKV + Proflavine vs. CHIKV + DMSO (Vero)'),
]

SPECIES_COMPARISONS_C636 = [
    (('C636 CHIKV',              'C636 DMSO'),         'C636_CHIKV+DMSO_vs_DMSO',
     'CHIKV + DMSO vs. DMSO (C636)'),
    (('C636 Simvastatin',        'C636 DMSO'),         'C636_Simvastatin_vs_DMSO',
     'Simvastatin vs. DMSO (C636 mock)'),
    (('C636 Simvastatin+CHIKV',  'C636 CHIKV'),        'C636_CHIKV+Simvastatin_vs_CHIKV+DMSO',
     'CHIKV + Simvastatin vs. CHIKV + DMSO (C636)'),
]

SPECIES_COMPARISONS = SPECIES_COMPARISONS_VERO + SPECIES_COMPARISONS_C636


def load_data(filepath):
    """
    Load protein-normalized lipidomics xlsx and return:
      df              : DataFrame, rows=lipid families (15), cols=individual sample names
      protein_amounts : Series, index=sample names, values=µg protein
      cond_to_samples : dict, base_condition_name -> [list of 3 sample column names]
      species_dfs     : dict {family_name -> DataFrame(species × samples)}
    """
    wb = openpyxl.load_workbook(filepath, read_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    wb.close()

    cond_row = rows[3]   # condition names
    prot_row = rows[4]   # protein amounts

    def normalize_name(raw):
        """Replace \xa0 and multiple whitespace, strip."""
        s = str(raw).replace('\xa0', ' ')
        return re.sub(r'\s+', ' ', s).strip()

    def base_name(norm):
        """Remove trailing R1/R2/R3 replicate suffix and normalize + signs."""
        b = re.sub(r'\s*R[123]\s*$', '', norm).strip()
        # Handle 'C636 Simvastatin + CHIKVR1' edge case (no space before R1)
        b = re.sub(r'R[123]$', '', b).strip()
        b = re.sub(r'\s*\+\s*CHIKV', '+CHIKV', b)
        return b

    # Build: col_index -> normalized sample name
    col_to_sample = {}
    protein_amounts = {}

    for col_idx in range(2, len(cond_row)):
        raw = cond_row[col_idx]
        if raw is None:
            break
        norm = normalize_name(raw)
        # Rebuild sample name as "BaseName R1/R2/R3" for consistency.
        # Raw names can be malformed (e.g. 'C636 Simvastatin + CHIKVR1' has no
        # space before R1), so we extract the replicate number and reformat.
        m = re.search(r'R([123])\s*$', norm)
        if m:
            clean_name = base_name(norm) + ' R' + m.group(1)
        else:
            clean_name = norm  # Méthanol — no replicate suffix
        col_to_sample[col_idx] = clean_name
        protein_amounts[clean_name] = prot_row[col_idx]

    protein_amounts = pd.Series(protein_amounts)

    # Group: base_condition -> [sample names]
    cond_to_samples = {}
    for col_idx, sample_name in col_to_sample.items():
        base = base_name(sample_name)
        cond_to_samples.setdefault(base, []).append(sample_name)

    # Extract 'Total of [Family]' rows AND individual species rows
    family_data  = {}
    species_raw  = {}   # family -> {species_name -> {sample_name -> value}}
    current_family = None

    for row in rows[7:]:
        if row[0] is not None:
            current_family = row[0]
        if row[1] is None:
            continue
        cell = str(row[1])
        row_data = {}
        for col_idx, sample_name in col_to_sample.items():
            val = row[col_idx]
            if val is None or (isinstance(val, (int, float)) and abs(val) < ND_THRESHOLD):
                row_data[sample_name] = np.nan
            else:
                row_data[sample_name] = float(val)
        if cell.startswith('Total of'):
            family_data[current_family] = row_data
        else:
            species_raw.setdefault(current_family, {})[cell] = row_data

    df = pd.DataFrame(family_data).T   # rows=families, cols=samples
    species_dfs = {
        fam: pd.DataFrame(spp).T        # rows=species, cols=samples
        for fam, spp in species_raw.items()
    }
    return df, protein_amounts, cond_to_samples, species_dfs


def compute_stats(df, cond_to_samples, comparisons):
    """
    Returns:
      means   : dict[group_name] -> pd.Series(family -> mean)
      stds    : dict[group_name] -> pd.Series(family -> std)
      results : dict[(treatment, control)] -> {
                    'label'  : str,
                    'log2fc' : pd.Series(family -> float),
                    'pvalue' : pd.Series(family -> float),
                    'stars'  : pd.Series(family -> str),
                }
    Statistical test: Welch's t-test (two-sided). With n=3 per group the
    minimum achievable p-value is ~0.017; Mann-Whitney cannot reach p<0.05
    with n=3 (min two-sided p = 0.1).
    """
    means = {}
    stds = {}

    for group, samples in cond_to_samples.items():
        valid_cols = [s for s in samples if s in df.columns]
        sub = df[valid_cols].astype(float)
        means[group] = sub.mean(axis=1, skipna=True)
        stds[group]  = sub.std(axis=1, skipna=True, ddof=1)

    results = {}

    for treatment, control, label in comparisons:
        if treatment not in cond_to_samples or control not in cond_to_samples:
            continue

        treat_cols = cond_to_samples[treatment]
        ctrl_cols  = cond_to_samples[control]

        log2fc_dict = {}
        pvalue_dict = {}
        stars_dict  = {}

        for family in df.index:
            treat_vals = df.loc[family, treat_cols].astype(float).dropna().values
            ctrl_vals  = df.loc[family, ctrl_cols].astype(float).dropna().values

            t_mean = treat_vals.mean() if len(treat_vals) > 0 else np.nan
            c_mean = ctrl_vals.mean()  if len(ctrl_vals)  > 0 else np.nan

            if np.isnan(t_mean) or np.isnan(c_mean) or c_mean == 0:
                log2fc_dict[family] = np.nan
                pvalue_dict[family] = np.nan
                stars_dict[family]  = ''
                continue

            log2fc_dict[family] = np.log2(t_mean / c_mean)

            if len(treat_vals) >= 2 and len(ctrl_vals) >= 2:
                _, p = stats.ttest_ind(treat_vals, ctrl_vals, equal_var=False)
                pvalue_dict[family] = p
                stars_dict[family]  = pval_to_stars(p)
            else:
                pvalue_dict[family] = np.nan
                stars_dict[family]  = ''

        results[(treatment, control)] = {
            'label':   label,
            'log2fc':  pd.Series(log2fc_dict),
            'pvalue':  pd.Series(pvalue_dict),
            'stars':   pd.Series(stars_dict),
        }

    return means, stds, results


def compute_species_stats(species_dfs, cond_to_samples, species_comparisons):
    """
    For each (family, comparison) pair, compute per-species log2FC and p-value.

    Parameters
    ----------
    species_dfs : dict {family_name -> DataFrame(species × samples)}
    cond_to_samples : dict {condition_name -> [sample_names]}
    species_comparisons : list of ((treatment, control), slug, label)

    Returns
    -------
    dict keyed by (family, treatment, control) -> {
        'log2fc': pd.Series(species -> float),
        'pvalue': pd.Series(species -> float),
        'stars':  pd.Series(species -> str),
        'label':  str,
        'slug':   str,
    }
    """
    results = {}
    for (treatment, control), slug, label in species_comparisons:
        if treatment not in cond_to_samples or control not in cond_to_samples:
            continue
        treat_cols = cond_to_samples[treatment]
        ctrl_cols  = cond_to_samples[control]

        for family, sp_df in species_dfs.items():
            valid_treat = [s for s in treat_cols if s in sp_df.columns]
            valid_ctrl  = [s for s in ctrl_cols  if s in sp_df.columns]

            log2fc_dict = {}
            pvalue_dict = {}
            stars_dict  = {}

            for species in sp_df.index:
                treat_vals = sp_df.loc[species, valid_treat].astype(float).dropna().values
                ctrl_vals  = sp_df.loc[species, valid_ctrl].astype(float).dropna().values

                t_mean = treat_vals.mean() if len(treat_vals) > 0 else np.nan
                c_mean = ctrl_vals.mean()  if len(ctrl_vals)  > 0 else np.nan

                if np.isnan(t_mean) or np.isnan(c_mean) or c_mean == 0:
                    log2fc_dict[species] = np.nan
                    pvalue_dict[species] = np.nan
                    stars_dict[species]  = ''
                    continue

                log2fc_dict[species] = np.log2(t_mean / c_mean)

                if len(treat_vals) >= 2 and len(ctrl_vals) >= 2:
                    _, p = stats.ttest_ind(treat_vals, ctrl_vals, equal_var=False)
                    pvalue_dict[species] = p
                    stars_dict[species]  = pval_to_stars(p)
                else:
                    pvalue_dict[species] = np.nan
                    stars_dict[species]  = ''

            results[(family, treatment, control)] = {
                'log2fc': pd.Series(log2fc_dict),
                'pvalue': pd.Series(pvalue_dict),
                'stars':  pd.Series(stars_dict),
                'label':  label,
                'slug':   slug,
            }
    return results


def save_figure(fig, output_dir, filename):
    os.makedirs(output_dir, exist_ok=True)
    stem = os.path.join(output_dir, filename)
    fig.savefig(f'{stem}.pdf', bbox_inches='tight')
    fig.savefig(f'{stem}.jpeg', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'  Saved {stem}.pdf / .jpeg')


def pval_to_stars(p):
    if p is None or np.isnan(p):
        return ''
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return 'ns'


def plot_qc(protein_amounts, output_dir):
    """
    Fig 1 — Horizontal bar chart of protein quantity per sample.
    Helps identify outliers in normalization.
    """
    # Drop Méthanol (negative value — instrument artifact)
    pa = protein_amounts.drop('Méthanol', errors='ignore').sort_values()

    # Color by condition group (first two words)
    palette = matplotlib.colormaps['tab20'].colors
    groups = sorted(set(' '.join(s.split()[:2]) for s in pa.index))
    color_map = {g: palette[i % 20] for i, g in enumerate(groups)}
    colors = [color_map[' '.join(s.split()[:2])] for s in pa.index]

    fig, ax = plt.subplots(figsize=(10, max(8, len(pa) * 0.3)))
    ax.barh(range(len(pa)), pa.values, color=colors, edgecolor='white', linewidth=0.5)
    ax.set_yticks(range(len(pa)))
    ax.set_yticklabels(pa.index, fontsize=8)
    ax.set_xlabel('Protein quantity (µg)', fontsize=11)
    ax.set_title('QC — Protein quantity per sample\n(used for normalization)', fontsize=12)
    ax.axvline(pa.median(), color='red', linestyle='--', linewidth=1,
               label=f'Median = {pa.median():.1f} µg')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    median_line = plt.Line2D([0], [0], color='red', linestyle='--', linewidth=1,
                             label=f'Median = {pa.median():.1f} µg')
    patches = [mpatches.Patch(color=color_map[g], label=g) for g in groups]
    ax.legend(handles=[median_line] + patches, bbox_to_anchor=(1.01, 1), loc='upper left',
              fontsize=7, title='Condition group')

    plt.tight_layout()
    save_figure(fig, output_dir, 'Fig1_QC_proteins')


def plot_composition(df, cond_to_samples, output_dir):
    """
    Fig 2 — Stacked bar chart (100%) of lipid family composition per condition group.
    Shows global remodelling of the lipidome across conditions.
    """
    # Build mean abundance per family per group (exclude Méthanol)
    groups = [g for g in cond_to_samples.keys() if g != 'Méthanol']
    comp = {}
    for group in groups:
        cols = [s for s in cond_to_samples[group] if s in df.columns]
        comp[group] = df[cols].mean(axis=1, skipna=True)

    comp_df = pd.DataFrame(comp)   # rows=families, cols=groups
    # Reorder families
    comp_df = comp_df.reindex([f for f in FAMILIES if f in comp_df.index])
    # Fill NaN with 0 for stacking
    comp_df = comp_df.fillna(0)
    # Normalize to 100%
    comp_pct = comp_df.div(comp_df.sum(axis=0), axis=1) * 100

    fig, ax = plt.subplots(figsize=(18, 6))
    bottom = np.zeros(len(groups))
    x = np.arange(len(groups))

    for family in comp_pct.index:
        values = comp_pct.loc[family, groups].values
        ax.bar(x, values, bottom=bottom, color=FAMILY_COLORS.get(family, '#aaaaaa'),
               label=family, width=0.8, edgecolor='white', linewidth=0.3)
        bottom += values

    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Lipid family composition (%)', fontsize=11)
    ax.set_title('Lipid family composition by condition group', fontsize=12)
    ax.set_ylim(0, 100)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8, title='Lipid family')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    save_figure(fig, output_dir, 'Fig2_Lipid_composition')


def plot_heatmap(stats_results, output_dir):
    """
    Fig 3 — Heatmap of log2FC for all 14 comparisons × 15 families.
    Cells annotated with significance stars. Color: blue (down) → red (up).
    """
    # Build matrices
    comparisons_keys = [(t, c) for t, c, _ in COMPARISONS if (t, c) in stats_results]
    col_labels = [stats_results[k]['label'] for k in comparisons_keys]
    families = [f for f in FAMILIES if f in next(iter(stats_results.values()))['log2fc'].index]

    fc_matrix   = pd.DataFrame(index=families, columns=col_labels, dtype=float)
    star_matrix = pd.DataFrame(index=families, columns=col_labels, dtype=str)

    for key, label in zip(comparisons_keys, col_labels):
        r = stats_results[key]
        for fam in families:
            fc_matrix.loc[fam, label]   = r['log2fc'].get(fam, np.nan)
            s = r['stars'].get(fam, '')
            star_matrix.loc[fam, label] = s if s in ('*', '**', '***') else ''

    # Clamp log2FC for color scale
    vmax = np.nanpercentile(np.abs(fc_matrix.values.astype(float)), 95)
    vmax = max(vmax, 1.0)

    fig, ax = plt.subplots(figsize=(max(14, len(col_labels) * 1.1), 7))
    sns.heatmap(
        fc_matrix.astype(float),
        ax=ax,
        cmap='RdBu_r',
        center=0,
        vmin=-vmax, vmax=vmax,
        annot=star_matrix.values,
        fmt='',
        linewidths=0.5,
        linecolor='white',
        cbar_kws={'label': 'log2 Fold Change', 'shrink': 0.6},
        annot_kws={'size': 9},
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha='right', fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)
    ax.set_title('log2 Fold Change — all comparisons\n(* p<0.05  ** p<0.01  *** p<0.001  ns = not significant)',
                 fontsize=11)
    ax.set_xlabel('')
    ax.set_ylabel('Lipid family', fontsize=10)

    plt.tight_layout()
    save_figure(fig, output_dir, 'Fig3_Heatmap_log2FC')


def _grouped_bar_two_groups(ax, df, cond_to_samples, group_a, group_b,
                             means, stds, stars_series,
                             label_a, label_b, color_a, color_b):
    """
    Draw grouped bars (group_a vs group_b) for all FAMILIES on ax.
    Error bars = SD. Stars drawn above the taller bar for significant results.
    """
    n = len(FAMILIES)
    x = np.arange(n)
    width = 0.35

    vals_a = [means[group_a].get(f, np.nan) if group_a in means else np.nan for f in FAMILIES]
    vals_b = [means[group_b].get(f, np.nan) if group_b in means else np.nan for f in FAMILIES]
    err_a  = [stds[group_a].get(f, 0) if group_a in stds else 0 for f in FAMILIES]
    err_b  = [stds[group_b].get(f, 0) if group_b in stds else 0 for f in FAMILIES]

    ax.bar(x - width/2, vals_a, width, yerr=err_a, label=label_a,
           color=color_a, capsize=3, error_kw={'linewidth': 1}, edgecolor='white')
    ax.bar(x + width/2, vals_b, width, yerr=err_b, label=label_b,
           color=color_b, capsize=3, error_kw={'linewidth': 1}, edgecolor='white')

    # Add stars above bars for significant results only
    for i, fam in enumerate(FAMILIES):
        star = stars_series.get(fam, '') if hasattr(stars_series, 'get') else ''
        if not star or star == 'ns':
            continue
        va = (vals_a[i] if not np.isnan(vals_a[i]) else 0) + (err_a[i] or 0)
        vb = (vals_b[i] if not np.isnan(vals_b[i]) else 0) + (err_b[i] or 0)
        y_top = max(va, vb) * 1.08
        ax.text(x[i], y_top, star, ha='center', va='bottom', fontsize=10, color='black')

    ax.set_xticks(x)
    ax.set_xticklabels(FAMILIES, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Normalized abundance\n(µg protein)', fontsize=9)
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def plot_chikv_effect(df, cond_to_samples, stats_results, output_dir):
    """
    Fig 4 — Effect of CHIKV infection on Vero cells.
    Small multiples: one subplot per lipid family, each with its own Y-scale.
    This avoids TG dominating a shared axis and crushing smaller but significant families.
    """
    key = ('Vero DMSO CHIKV', 'Vero DMSO')
    if key not in stats_results:
        print('  WARNING: comparison not found, skipping Fig 4')
        return

    means, stds, _ = compute_stats(df, cond_to_samples, COMPARISONS)
    stars = stats_results[key]['stars']

    n_cols = 5
    n_rows = 3   # 15 families → 3×5 grid
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 10))
    axes = axes.flatten()

    groups  = ['Vero DMSO', 'Vero DMSO CHIKV']
    labels  = ['Mock', 'CHIKV']
    colors  = ['#4878CF', '#D65F5F']
    x       = np.arange(len(groups))
    width   = 0.5

    for i, fam in enumerate(FAMILIES):
        ax = axes[i]
        vals = [means[g].get(fam, np.nan) for g in groups]
        errs = [stds[g].get(fam, 0)       for g in groups]

        # Plot bars replacing NaN with 0 (ND bars remain invisible at height 0)
        plot_vals = [0 if np.isnan(v) else v for v in vals]
        plot_errs = [0 if np.isnan(v) else e for v, e in zip(vals, errs)]
        ax.bar(x, plot_vals, width, yerr=plot_errs, color=colors,
               capsize=4, error_kw={'linewidth': 1}, edgecolor='white')

        ax.set_title(fam, fontsize=9, fontweight='bold',
                     color=FAMILY_COLORS.get(fam, '#333333'))
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel('Abundance\n(norm. protein)', fontsize=6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(bottom=0)

        # ND annotation for bars where mean is not detected
        y_range = ax.get_ylim()[1]
        y_nd = y_range * 0.05 if y_range > 0 else 0.5
        for j, v in enumerate(vals):
            if np.isnan(v):
                ax.text(x[j], y_nd, 'ND', ha='center', va='bottom',
                        fontsize=7, color='gray', style='italic')

        # Significance star above the taller bar
        star = stars.get(fam, '')
        if star and star != 'ns':
            tops = [plot_vals[j] + (plot_errs[j] or 0) for j in range(len(groups))]
            y_top = max(tops)
            ax.text(x.mean(), y_top * 1.12, star, ha='center', va='bottom',
                    fontsize=11, color='black')

    # Hide unused subplots (grid has 15 cells for 15 families — none unused)
    for j in range(len(FAMILIES), len(axes)):
        axes[j].set_visible(False)

    # Shared legend
    patches = [mpatches.Patch(color=c, label=l) for c, l in zip(colors, ['Vero DMSO (mock)', 'Vero DMSO CHIKV (infected)'])]
    fig.legend(handles=patches, loc='lower right', fontsize=9, frameon=True)

    fig.suptitle('Effect of CHIKV infection on Vero cells — per lipid family\n'
                 '(* p<0.05  ** p<0.01  *** p<0.001  — Welch t-test, n=3  |  each panel: independent Y-scale)',
                 fontsize=11, y=1.01)
    plt.tight_layout()
    save_figure(fig, output_dir, 'Fig4_CHIKV_effect_Vero')


def _drug_heatmap(stats_results, all_cols, title, filename, output_dir):
    """
    Heatmap log2FC: rows = 15 lipid families, cols = all_cols.
    all_cols: list of ((treatment, control), label) — ordered column definitions.
    """
    col_labels  = [label for _, label in all_cols]
    fc_matrix   = pd.DataFrame(index=FAMILIES, columns=col_labels, dtype=float)
    star_matrix = pd.DataFrame(index=FAMILIES, columns=col_labels, dtype=str)

    for key, label in all_cols:
        if key not in stats_results:
            continue
        r = stats_results[key]
        for fam in FAMILIES:
            fc_matrix.loc[fam, label]   = r['log2fc'].get(fam, np.nan)
            s = r['stars'].get(fam, '')
            star_matrix.loc[fam, label] = s if s in ('*', '**', '***') else ''

    vmax = np.nanpercentile(np.abs(fc_matrix.values.astype(float)), 95)
    vmax = max(vmax, 1.0)

    fig, ax = plt.subplots(figsize=(max(8, len(col_labels) * 1.2), 7))
    sns.heatmap(
        fc_matrix.astype(float),
        ax=ax,
        cmap='RdBu_r',
        center=0,
        vmin=-vmax, vmax=vmax,
        annot=star_matrix.values,
        fmt='',
        linewidths=0.5,
        linecolor='white',
        cbar_kws={'label': 'log2 Fold Change', 'shrink': 0.7},
        annot_kws={'size': 10},
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha='right', fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)
    ax.set_xlabel('')
    ax.set_ylabel('Lipid family', fontsize=11)
    ax.set_title(title + '\n(* p<0.05  ** p<0.01  *** p<0.001  — Welch t-test, n=3)',
                 fontsize=11)
    plt.tight_layout()
    save_figure(fig, output_dir, filename)


def _drug_supplement_grid(df, cond_to_samples, stats_results, means, stds,
                           drug_conditions, control, color_ctrl, color_drug,
                           title, filename, output_dir):
    """
    Supplement grid: rows = 15 families, cols = 5 drugs.
    Each cell: 2 bars (DMSO vs drug) ± SD, independent Y-scale.
    Width ~12", height ~25" — suitable for PDF supplement.
    """
    n_fams  = len(FAMILIES)
    n_drugs = len(DRUGS_VERO)
    fig, axes = plt.subplots(n_fams, n_drugs, figsize=(13, 28))
    x     = np.array([0, 1])
    width = 0.55

    for col, drug in enumerate(DRUGS_VERO):
        treatment = drug_conditions[drug]
        key       = (treatment, control)

        for row, fam in enumerate(FAMILIES):
            ax = axes[row, col]

            if key not in stats_results or treatment not in cond_to_samples:
                ax.set_visible(False)
                continue

            vals = [means[control].get(fam, np.nan),
                    means[treatment].get(fam, np.nan)]
            errs = [stds[control].get(fam, 0),
                    stds[treatment].get(fam, 0)]

            plot_vals = [0 if np.isnan(v) else v for v in vals]
            plot_errs = [0 if np.isnan(v) else e for v, e in zip(vals, errs)]
            ax.bar(x, plot_vals, width, yerr=plot_errs,
                   color=[color_ctrl, color_drug],
                   capsize=3, error_kw={'linewidth': 0.8}, edgecolor='white')

            ax.set_ylim(bottom=0)
            y_range = ax.get_ylim()[1]
            y_nd = y_range * 0.05 if y_range > 0 else 0.5
            for j, v in enumerate(vals):
                if np.isnan(v):
                    ax.text(x[j], y_nd, 'ND', ha='center', va='bottom',
                            fontsize=6, color='gray', style='italic')

            # PRISM-style significance bracket
            star = stats_results[key]['stars'].get(fam, '')
            if star and star != 'ns':
                tops = [plot_vals[j] + (plot_errs[j] or 0) for j in range(len(vals))]
                y_max = max(tops) if max(tops) > 0 else 1
                y_tick_bot = y_max * 1.08
                y_bracket  = y_max * 1.20
                ax.plot([x[0], x[1]], [y_bracket, y_bracket],
                        color='black', linewidth=0.8)
                ax.plot([x[0], x[0]], [y_tick_bot, y_bracket],
                        color='black', linewidth=0.8)
                ax.plot([x[1], x[1]], [y_tick_bot, y_bracket],
                        color='black', linewidth=0.8)
                ax.text((x[0] + x[1]) / 2, y_bracket, star,
                        ha='center', va='bottom', fontsize=7, color='black')
                ax.set_ylim(top=max(ax.get_ylim()[1], y_bracket * 1.35))

            # Drug name on top row only
            if row == 0:
                ax.set_title(drug, fontsize=8, fontweight='bold')

            # Family name on leftmost column only
            if col == 0:
                ax.set_ylabel(fam, fontsize=7, fontweight='bold',
                              color=FAMILY_COLORS.get(fam, '#333333'),
                              rotation=0, labelpad=40, va='center')

            ax.set_xticks(x)
            ax.set_xticklabels(['ctrl', 'drug'], fontsize=5)
            ax.tick_params(axis='y', labelsize=5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

    patches = [
        mpatches.Patch(color=color_ctrl, label='DMSO (control)'),
        mpatches.Patch(color=color_drug, label='Drug treatment'),
    ]
    fig.legend(handles=patches, loc='lower right', fontsize=9, frameon=True)
    fig.suptitle(title + '\n(supplement — each cell: independent Y-scale  |'
                 '  * p<0.05  ** p<0.01  *** p<0.001)',
                 fontsize=11, y=1.002)
    plt.tight_layout()
    save_figure(fig, output_dir, filename)


def save_fig5_csv(stats_results, all_cols_vero, all_cols_c636, output_dir):
    """
    Export log2FC, p-value and significance stars for Fig 5 comparisons to CSV.
    Format: tidy/long — one row per (lipid_family, comparison).
    Two files: Fig5_log2FC_Vero.csv and Fig5_log2FC_C636.csv
    """
    os.makedirs(output_dir, exist_ok=True)

    for all_cols, fname in [(all_cols_vero, 'Fig5_log2FC_Vero'),
                            (all_cols_c636, 'Fig5_log2FC_C636')]:
        rows = []
        for (treatment, control), label in all_cols:
            key = (treatment, control)
            if key not in stats_results:
                continue
            r = stats_results[key]
            for fam in FAMILIES:
                rows.append({
                    'lipid_family': fam,
                    'comparison':   label.replace('\n', ' '),
                    'treatment':    treatment,
                    'control':      control,
                    'log2FC':       r['log2fc'].get(fam, float('nan')),
                    'pvalue':       r['pvalue'].get(fam, float('nan')),
                    'significance': r['stars'].get(fam, ''),
                })
        out = pd.DataFrame(rows, columns=[
            'lipid_family', 'comparison', 'treatment', 'control',
            'log2FC', 'pvalue', 'significance'
        ])
        csv_path = os.path.join(output_dir, f'{fname}.csv')
        out.to_csv(csv_path, index=False, float_format='%.6g')
        print(f'  Saved {csv_path}')


def plot_drug_effects_infected(df, cond_to_samples, stats_results, output_dir):
    """
    Fig 5_Vero (main) — Heatmap: 15 families × 11 cols (Vero only).
    Fig 5_C636 (main) — Heatmap: 15 families × 4 cols (C636 + Vero DMSO reference).
    Fig 5S (suppl)    — Bar chart grid: 15 rows × 5 cols, drug+CHIKV vs DMSO CHIKV.
    """
    means, stds, _ = compute_stats(df, cond_to_samples, COMPARISONS)
    drug_conditions = {drug: f'Vero {drug} CHIKV' for drug in DRUGS_VERO}

    # ── Vero figure ───────────────────────────────────────────────────────────
    all_cols_vero = [
        (('Vero DMSO CHIKV',         'Vero DMSO'),          'CHIKV + DMSO\nvs. DMSO'),
        (('Vero HMN214 CHIKV',       'Vero DMSO CHIKV'),    'HMN214\nvs. DMSO'),
        (('Vero HMN214 CHIKV',       'Vero HMN214'),        'CHIKV + HMN214\nvs. CHIKV + DMSO'),
        (('Vero Simvastatin CHIKV',  'Vero DMSO CHIKV'),    'Simvastatin\nvs. DMSO'),
        (('Vero Simvastatin CHIKV',  'Vero Simvastatin'),   'CHIKV + Simvastatin\nvs. CHIKV + DMSO'),
        (('Vero Atorvastatin CHIKV', 'Vero DMSO CHIKV'),    'Atorvastatin\nvs. DMSO'),
        (('Vero Atorvastatin CHIKV', 'Vero Atorvastatin'),  'CHIKV + Atorvastatin\nvs. CHIKV + DMSO'),
        (('Vero Fenretinide CHIKV',  'Vero DMSO CHIKV'),    'Fenretinide\nvs. DMSO'),
        (('Vero Fenretinide CHIKV',  'Vero Fenretinide'),   'CHIKV + Fenretinide\nvs. CHIKV + DMSO'),
        (('Vero Proflavine CHIKV',   'Vero DMSO CHIKV'),    'Proflavine\nvs. DMSO'),
        (('Vero Proflavine CHIKV',   'Vero Proflavine'),    'CHIKV + Proflavine\nvs. CHIKV + DMSO'),
    ]

    _drug_heatmap(
        stats_results, all_cols_vero,
        title='Drug effects on CHIKV-infected Vero cells',
        filename='Fig5_Drug_effects_Vero_infected',
        output_dir=output_dir,
    )

    # ── C636 figure (+ Vero DMSO reference column) ───────────────────────────
    all_cols_c636 = [
        (('Vero DMSO CHIKV',         'Vero DMSO'),       'CHIKV + DMSO\nvs. DMSO\n(Vero ref.)'),
        (('C636 CHIKV',              'C636 DMSO'),       'CHIKV + DMSO\nvs. DMSO\n(C636)'),
        (('C636 Simvastatin+CHIKV',  'C636 CHIKV'),      'Simvastatin\nvs. DMSO\n(C636)'),
        (('C636 Simvastatin+CHIKV',  'C636 DMSO'),       'CHIKV + Simvastatin\nvs. CHIKV + DMSO\n(C636)'),
    ]

    _drug_heatmap(
        stats_results, all_cols_c636,
        title='Drug effects on CHIKV-infected C636 cells',
        filename='Fig5_Drug_effects_C636_infected',
        output_dir=output_dir,
    )

    # ── CSV export of log2FC values for Fig 5 ────────────────────────────────
    save_fig5_csv(stats_results, all_cols_vero, all_cols_c636, output_dir)

    # ── Supplement: drug+CHIKV vs DMSO CHIKV in Vero (Fig 5S) ───────────────
    _drug_supplement_grid(
        df, cond_to_samples, stats_results, means, stds,
        drug_conditions=drug_conditions,
        control='Vero DMSO CHIKV',
        color_ctrl='#D65F5F', color_drug='#6ABF69',
        title='Drug effects on CHIKV-infected Vero cells (drug+CHIKV vs CHIKV + DMSO)',
        filename='Fig5S_Drug_effects_Vero_infected_supplement',
        output_dir=output_dir,
    )


def plot_drug_effects_mock(df, cond_to_samples, stats_results, output_dir):
    """
    Fig 6 (main)   — Heatmap log2FC: 15 families × 5 drugs, mock Vero (7–11 vs CTRL1).
    Fig 6S (suppl) — Bar chart grid: 15 rows × 5 cols, independent Y-scale per cell.
    """
    means, stds, _ = compute_stats(df, cond_to_samples, COMPARISONS)
    drug_conditions = {drug: f'Vero {drug}' for drug in DRUGS_VERO}

    all_cols = [
        (('Vero HMN214',       'Vero DMSO'), 'HMN214'),
        (('Vero Simvastatin',  'Vero DMSO'), 'Simvastatin'),
        (('Vero Atorvastatin', 'Vero DMSO'), 'Atorvastatin'),
        (('Vero Fenretinide',  'Vero DMSO'), 'Fenretinide'),
        (('Vero Proflavine',   'Vero DMSO'), 'Proflavine'),
    ]

    _drug_heatmap(
        stats_results, all_cols,
        title='Drug effects on mock (non-infected) Vero cells (vs DMSO mock)',
        filename='Fig6_Drug_effects_Vero_mock',
        output_dir=output_dir,
    )
    _drug_supplement_grid(
        df, cond_to_samples, stats_results, means, stds,
        drug_conditions=drug_conditions,
        control='Vero DMSO',
        color_ctrl='#4878CF', color_drug='#E8A838',
        title='Drug effects on mock (non-infected) Vero cells (vs DMSO mock)',
        filename='Fig6S_Drug_effects_Vero_mock_supplement',
        output_dir=output_dir,
    )


def plot_c636(df, cond_to_samples, stats_results, output_dir):
    """
    Fig 7 — C636 cell line analysis. 3 panels:
      A: C636 CHIKV vs C636 DMSO             (effect of infection)
      B: C636 Simvastatin+CHIKV vs C636 CHIKV (Simvastatin on infected C636)
      C: C636 Simvastatin vs C636 DMSO        (Simvastatin on mock C636)
    """
    means, stds, _ = compute_stats(df, cond_to_samples, COMPARISONS)

    panels = [
        ('C636 CHIKV',             'C636 DMSO',    'A — CHIKV effect (C636)',
         '#4878CF', '#D65F5F', 'C636 DMSO', 'C636 CHIKV'),
        ('C636 Simvastatin+CHIKV', 'C636 DMSO',    'B — Simvastatin on infected C636',
         '#4878CF', '#6ABF69', 'C636 DMSO', 'C636 Simvastatin+CHIKV'),
        ('C636 Simvastatin',       'C636 DMSO',    'C — Simvastatin on mock C636',
         '#4878CF', '#E8A838', 'C636 DMSO', 'C636 Simvastatin'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(22, 5))

    for ax, (treatment, control, title, col_a, col_b, lab_a, lab_b) in zip(axes, panels):
        key = (treatment, control)
        if key not in stats_results or treatment not in cond_to_samples:
            ax.set_title(f'{title}\n(data not available)', fontsize=9)
            continue

        stars = stats_results[key]['stars']
        _grouped_bar_two_groups(
            ax, df, cond_to_samples,
            group_a=control, group_b=treatment,
            means=means, stds=stds, stars_series=stars,
            label_a=lab_a, label_b=lab_b,
            color_a=col_a, color_b=col_b,
        )
        ax.set_title(title, fontsize=10, fontweight='bold')

    fig.suptitle('C636 cell line — lipid family analysis\n'
                 '(* p<0.05  ** p<0.01  *** p<0.001  — Welch t-test, n=3)',
                 fontsize=12)
    plt.tight_layout()
    save_figure(fig, output_dir, 'Fig7_C636_analysis')


def plot_species_heatmaps(species_stats, output_dir):
    """
    Generate 2 heatmaps per lipid family (Vero + C636) — 30 figures total.
    Rows = individual lipid species, columns = comparisons.
    Same RdBu_r style as _drug_heatmap / Fig5.

    Output: analysis/CHIKV/species/{FAMILY}_Vero.pdf/.jpeg
                                   {FAMILY}_C636.pdf/.jpeg
    """
    species_dir = os.path.join(output_dir, 'species')
    os.makedirs(species_dir, exist_ok=True)

    # Ordered list of families (preserve insertion order from species_stats keys)
    families = list(dict.fromkeys(fam for fam, _, _ in species_stats))

    configs = [
        (SPECIES_COMPARISONS_VERO, 'Vero',
         'Drug effects — Vero cells — species level'),
        (SPECIES_COMPARISONS_C636, 'C636',
         'Drug effects — C636 cells — species level'),
    ]

    total = len(families) * len(configs)
    done  = 0

    for family in families:
        for comparisons_list, suffix, title_prefix in configs:
            # Get ordered species list from first available comparison
            species_list = None
            for (treatment, control), slug, label in comparisons_list:
                key = (family, treatment, control)
                if key in species_stats:
                    species_list = species_stats[key]['log2fc'].index.tolist()
                    break
            if species_list is None:
                done += 1
                continue

            col_labels = [label for (_, _), _, label in comparisons_list]
            fc_matrix   = pd.DataFrame(np.nan, index=species_list, columns=col_labels)
            star_matrix = pd.DataFrame('',     index=species_list, columns=col_labels)

            for (treatment, control), slug, label in comparisons_list:
                key = (family, treatment, control)
                if key not in species_stats:
                    continue
                r = species_stats[key]
                for sp in species_list:
                    fc_matrix.loc[sp, label]   = r['log2fc'].get(sp, np.nan)
                    s = r['stars'].get(sp, '')
                    star_matrix.loc[sp, label] = s if s in ('*', '**', '***') else ''

            # Adaptive figure size: ~0.35" per species (height), ~1.4" per comparison (width)
            fig_height = max(4.0, len(species_list) * 0.35)
            fig_width  = max(8.0, len(col_labels) * 1.4)

            vmax = np.nanpercentile(np.abs(fc_matrix.values.astype(float)), 95)
            vmax = max(vmax, 1.0)

            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            sns.heatmap(
                fc_matrix.astype(float),
                ax=ax,
                cmap='RdBu_r',
                center=0,
                vmin=-vmax, vmax=vmax,
                annot=star_matrix.values,
                fmt='',
                linewidths=0.4,
                linecolor='white',
                cbar_kws={'label': 'log2 Fold Change', 'shrink': 0.6},
                annot_kws={'size': 7},
            )
            ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha='right', fontsize=8)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)
            ax.set_xlabel('')
            ax.set_ylabel('Lipid species', fontsize=10)
            ax.set_title(
                f'{title_prefix}\n{family}\n'
                f'(* p<0.05  ** p<0.01  *** p<0.001  — Welch t-test, n=3)',
                fontsize=10,
            )
            plt.tight_layout()
            save_figure(fig, species_dir, f'{family}_{suffix}')

            done += 1
            if done % 10 == 0 or done == total:
                print(f'    {done}/{total} species heatmaps done')


def main():
    t0 = time.time()
    print('=== CHIKV Lipidomics Analysis ===')
    print(f'Input : {DATA_FILE}')
    print(f'Output: {OUTPUT_DIR}/')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print('\n[1/5] Loading data...')
    df, protein_amounts, cond_to_samples, species_dfs = load_data(DATA_FILE)
    print(f'      {len(df)} lipid families x {len(df.columns)} samples')
    print(f'      {len(cond_to_samples)} condition groups')

    print('[2/5] Computing statistics (family level)...')
    means, stds, stats_results = compute_stats(df, cond_to_samples, COMPARISONS)
    print(f'      {len(stats_results)} comparisons computed')

    print('[3/5] Computing statistics (species level)...')
    species_stats = compute_species_stats(species_dfs, cond_to_samples, SPECIES_COMPARISONS)
    print(f'      {len(species_stats)} (family × comparison) pairs computed')

    print('[4/5] Generating figures...')
    plot_qc(protein_amounts, OUTPUT_DIR)
    plot_composition(df, cond_to_samples, OUTPUT_DIR)
    plot_heatmap(stats_results, OUTPUT_DIR)
    plot_chikv_effect(df, cond_to_samples, stats_results, OUTPUT_DIR)
    plot_drug_effects_infected(df, cond_to_samples, stats_results, OUTPUT_DIR)
    plot_drug_effects_mock(df, cond_to_samples, stats_results, OUTPUT_DIR)
    plot_c636(df, cond_to_samples, stats_results, OUTPUT_DIR)

    print('[5/5] Generating species-level heatmaps (30 figures)...')
    plot_species_heatmaps(species_stats, OUTPUT_DIR)

    elapsed = time.time() - t0
    print(f'\nDone in {elapsed:.1f}s — 10 family figures + 30 species heatmaps'
          f' saved to {OUTPUT_DIR}/')


if __name__ == '__main__':
    main()
