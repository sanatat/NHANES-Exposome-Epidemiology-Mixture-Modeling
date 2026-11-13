#!/usr/bin/env python3
"""Exploratory human-exposome report, never a claim of causal effects."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
PROC = ROOT / 'data' / 'processed'
RESULTS = ROOT / 'results'
RESULTS.mkdir(exist_ok=True)


def main():
    qc = json.loads((PROC / 'qc_report.json').read_text())
    joint = (PROC / 'analysis_joint.csv').exists() and qc.get('joint_weights_verified')
    filename = 'analysis_joint.csv' if joint else 'analysis_pfas.csv'
    df = pd.read_csv(PROC / filename)
    chemicals = qc['joint_analytical_exposures'] if joint else qc['pfas_analytical_exposures']
    print(f'Examining {len(df)} adults, {len(chemicals)} exposures in {filename}')

    summary = []
    weights = df['survey_weight'].to_numpy(float)
    for column in chemicals:
        x = df[column].to_numpy(dtype=float)
        good = np.isfinite(x) & (x > 0) & np.isfinite(weights) & (weights > 0)
        gm = float(np.exp(np.average(np.log(x[good]), weights=weights[good]))) if good.any() else np.nan
        summary.append({'chemical': column, 'specimen': 'serum' if column.startswith('LBX') else 'urine',
                        'units': 'ng/mL', 'n_valid': int(good.sum()),
                        'median_unweighted': float(np.nanmedian(x)),
                        'weighted_geometric_mean': gm})
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(RESULTS / 'exposure_descriptives.csv', index=False)

    logex = np.log2(df[chemicals])
    correlations = logex.corr(method='spearman')
    correlations.to_csv(RESULTS / 'spearman_matrix.csv')
    fig, ax = plt.subplots(figsize=(max(7, len(chemicals) * .9), 6))
    im = ax.imshow(correlations, vmin=-1, vmax=1, cmap='coolwarm', aspect='equal')
    ax.set_xticks(range(len(chemicals)), chemicals, rotation=45, ha='right')
    ax.set_yticks(range(len(chemicals)), chemicals)
    for i in range(len(chemicals)):
        for j in range(len(chemicals)):
            ax.text(j, i, f'{correlations.iloc[i,j]:.2f}', ha='center', va='center', fontsize=9)
    fig.colorbar(im, ax=ax, label='Spearman rho')
    ax.set_title(f'NHANES 2017–2018 chemical co-exposures (n={len(df)})')
    fig.tight_layout()
    fig.savefig(RESULTS / 'exposome_spearman_heatmap.png', dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    y = df['LBXHSCRP'].to_numpy(float)
    ax.hist(y[np.isfinite(y) & (y > 0)], bins=32, color='slategray', edgecolor='white')
    ax.set_xscale('log')
    ax.set_xlabel('hs-CRP (mg/L, logarithmic x-axis)')
    ax.set_ylabel('Participants (unweighted count)')
    ax.set_title('Distribution of inflammation biomarker, hs-CRP')
    fig.tight_layout()
    fig.savefig(RESULTS / 'hscrp_distribution.png', dpi=170)
    plt.close(fig)
    (RESULTS / 'eda_notes.txt').write_text(
        'Source: CDC NHANES 2017–2018\n'
        f'Analysis cohort: {filename}, n={len(df)}\n'
        'Exposures: ' + ', '.join(chemicals) + '\n'
        'Outcome: hs-CRP (inflammation; NOT cancer)\n'
        'Plots descriptive only; prevalence inference needs complex sampling design.\n'
        'Urine analytes are adjusted using log urine creatinine as a model covariate, not divided by it.\n'
        'Non-detects are already provided as LOD/sqrt(2) by CDC.\n'
    )
    print(summary_df.to_string(index=False))
    print(f'Plots and descriptives saved in: {RESULTS}')


if __name__ == '__main__':
    main()
