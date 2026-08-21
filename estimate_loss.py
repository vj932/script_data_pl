#!/usr/bin/env python3
# Viraj Arora

import os
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
import scipy as sp # type: ignore

from get_loss import findPlasmidLossRate


# Hardcoded colony counts: dark colonies (successes) out of total colonies (trials)
COUNTS = {
    'Carbenicillin': {'dark': 42, 'total': 5329},
    'Kanamycin':     {'dark': 9,  'total': 4448},
}

N_ITER = 10000
SEED = 42


def main():
    '''This script estimates plasmid loss rates using a parametric bootstrap approach. Here, the triplicates are pooled over each treatmenet, and each data point is assumed to be i.i.d. '''
    script_dir = os.path.dirname(os.path.abspath(__file__))
    odpcb_file = os.path.join(script_dir, "OD_values_growth_curve_uv_vs_ck_070126.txt")
    viability_file = os.path.join(script_dir, "NORMALIZED_CFU_vs_time.txt")

    rng = np.random.default_rng(SEED)

    # Parametric bootstrap: for each iteration draw k ~ Binomial(total, p_hat) and use p = k / total
    samples = {}
    for treatment, c in COUNTS.items():
        p_hat = c['dark'] / c['total']
        draws = rng.binomial(c['total'], p_hat, size=N_ITER)
        samples[treatment] = draws / c['total']

    # Feed sampled proportions into the same getLoss path as get_loss.py
    prop_list = [(float(p), t) for t, ps in samples.items() for p in ps]

    obj = findPlasmidLossRate(odpcb_file, viability_file, prop_list, treatment='Carbenicillin')
    rates = obj.getLosses()

    stats = {}
    for head, vals in rates.items():
        if len(vals) == 0:
            continue
        stats[head] = [np.mean(vals), np.var(vals, ddof=1)]

    lines = ["Plasmid Loss Rates (per division) [Monte Carlo, parametric bootstrap]:"]
    for head, (m, v) in stats.items():
        lines.append(f"{head}: mean={m:.5f}  var={v:.5g}  (n={len(rates[head])})")
    output_text = "\n".join(lines)
    print(output_text)

    results_path = os.path.join(script_dir, "results_mc.txt")
    with open(results_path, 'w') as f:
        f.write(output_text + "\n")

    # Bar plot in the same style as findPlasmidLossRate.plotLosses, error bars = sqrt(sample variance)
    plt.bar(stats.keys(),
            [v[0] for v in stats.values()],
            yerr=[np.sqrt(v[1]) for v in stats.values()],
            capsize=5, color='gray', edgecolor='black')
    plt.ylabel("Loss Rate (per division)")
    plt.xlabel("Treatment")
    plt.title("Plasmid Loss Rates by Antibiotic (MC)")
    plt.tight_layout()
    plt.savefig(os.path.join(script_dir, "plot_mc.png"), bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    main()
