import json, sys

d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'bench/ipm/real_corpus_results_final.json'))
ms = d['models']
rows = []
for m in ms:
    name = m['model']
    ref = m.get('reference_ipm') or {}
    refst = ref.get('status')
    zh = m.get('zhighs') or {}
    zhst = zh.get('status')
    zt = zh.get('solve_ns'); rt = ref.get('solve_ns')
    ratio = (zt / rt) if (zt and rt) else None
    rows.append((name, refst, zhst, m.get('accepted'), ratio, m.get('rows'), m.get('cols')))
print(f"{'model':12} {'ref':10} {'zhighs':18} {'acc':5} {'ratio':>8} {'rows':>6} {'cols':>6}")
for r in sorted(rows, key=lambda x: -(x[4] or 0)):
    ratio_s = ('%.2f' % r[4]) if r[4] else '-'
    print(f"{r[0]:12} {str(r[1]):10} {str(r[2]):18} {str(r[3]):5} {ratio_s:>8} {r[5] or 0:>6} {r[6] or 0:>6}")
