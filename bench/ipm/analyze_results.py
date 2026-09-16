"""Summarize a corpus validation report produced by validate_corpus.py.

Compares zhighs IPM against the HiGHS references recorded in the same JSON:
status distribution, validated-accuracy count, per-model timing ratios, and
per-phase time breakdown for the zhighs solves.
"""
import argparse
import json
import math
from collections import Counter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=str)
    args = parser.parse_args()
    with open(args.report, encoding='utf-8') as handle:
        report = json.load(handle)

    models = report['models']
    status = Counter(m['zhighs']['status'] for m in models)
    print(f"models: {len(models)}  backend: {report['backend']}")
    print(f"zhighs status distribution: {dict(status)}")
    print(f"accepted (KKT + dimensions): {sum(m['accepted'] for m in models)}/{len(models)}")
    ref_status = Counter(m['reference_simplex']['status'] for m in models)
    print(f"reference (HiGHS simplex) status: {dict(ref_status)}")
    print()

    header = (f"{'model':12} {'rows':>4}x{'cols':<5} {'zh_it':>5} {'hi_it':>4} "
              f"{'zh_ms':>9} {'hiipm_ms':>8} {'hispx_ms':>8} {'x':>5} {'val':>3}")
    print(header)
    ratios = []
    for m in models:
        z, hi, hs = m['zhighs'], m['reference_ipm'], m['reference_simplex']
        ratio = ''
        if z['status'] == 'optimal' and hi['status'] == 'Optimal' and hi['solve_ns']:
            value = z['solve_ns'] / hi['solve_ns']
            ratios.append(value)
            ratio = f'{value:4.2f}'
        validated = m.get('validation', {}).get('accepted', False)
        print(f"{m['model']:12} {m['rows']:>4}x{m['cols']:<5} {z['iterations']:>5} "
              f"{hi['ipm_iterations']:>4} {z['solve_ns'] / 1e6:>9.2f} {hi['solve_ns'] / 1e6:>8.2f} "
              f"{hs['solve_ns'] / 1e6:>8.2f} {ratio:>5} {'ok' if validated else '--':>3}")

    if ratios:
        geo = math.exp(sum(math.log(x) for x in ratios) / len(ratios))
        best, worst = min(ratios), max(ratios)
        wins = sum(x < 1.0 for x in ratios)
        print(f"\ncomparable optimal-optimal: {len(ratios)} models; zhighs/HiGHS-IPM")
        print(f"  geometric mean {geo:.2f}x  best {best:.3f}x  worst {worst:.1f}x  zhighs faster on {wins}")

    # Per-phase breakdown over successful solves.
    totals = Counter()
    for m in models:
        stats = m['zhighs'].get('stats')
        if m['zhighs']['status'] == 'optimal' and stats:
            for key, value in stats.items():
                totals[key] += value
    factor_share = totals['factor_ns'] / max(totals['total_ns'], 1)
    assembly_share = totals['assembly_ns'] / max(totals['total_ns'], 1)
    direction_share = totals['direction_ns'] / max(totals['total_ns'], 1)
    print(f"\nsuccessful-solve time shares: factor {factor_share:.1%}, "
          f"assembly {assembly_share:.1%}, direction {direction_share:.1%}")

    # Stored-normal density of the successful solves: the profile (envelope)
    # factorization stores rows*first-row extents, so density near 1.0 means
    # near-dense fill-in after the structural envelope is respected.
    print(f"\n{'model':12} {'rows':>5} {'nnz':>6} {'stored':>8} {'m*m':>8} {'fill_ratio':>9}")
    for m in sorted(models, key=lambda m: -m['zhighs'].get('stats', {}).get('normal_entries', 0))[:12]:
        stats = m['zhighs'].get('stats')
        if not stats:
            continue
        stored = stats['normal_entries']
        print(f"{m['model']:12} {m['rows']:>5} {m['nnz']:>6} {stored:>8} {m['rows'] ** 2:>8} "
              f"{stored / m['rows'] ** 2:>9.3f}")

    # Failure detail: which optimal-reference models missed acceptance.
    misses = [m for m in models
              if m['reference_simplex']['status'] == 'Optimal'
              and m['zhighs']['status'] == 'optimal'
              and not m['accepted']]
    if misses:
        print(f"\noptimal but not accepted ({len(misses)}):")
        for m in misses:
            failed = m.get('validation', {}).get('failed_checks', ['no-validation'])
            print(f"  {m['model']:12} {', '.join(failed)}")


if __name__ == '__main__':
    main()