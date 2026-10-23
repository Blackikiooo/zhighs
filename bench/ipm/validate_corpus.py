"""Validate every LP in a local corpus with independent HiGHS parsing/KKT checks.

No synthetic models, no discarded failures, no solver fallback for zhighs.
Use a clean Python process and --highspy-path for an isolated reference install.
"""
import argparse
import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path


def validate(lp, answer, reference_objective, np):
    x, y, rc = [np.asarray(answer[key], dtype=float) for key in ('primal', 'dual', 'reduced_cost')]
    if (len(x), len(y), len(rc)) != (lp.num_col_, lp.num_row_, lp.num_col_):
        return {'accepted': False, 'reason': 'dimension_mismatch'}
    if not all(np.isfinite(v).all() for v in (x, y, rc)):
        return {'accepted': False, 'reason': 'nonfinite_solution'}
    c = np.asarray(lp.col_cost_)
    ax = np.zeros(lp.num_row_)
    aty = np.zeros(lp.num_col_)
    aty_abs = np.zeros(lp.num_col_)
    mat = lp.a_matrix_
    for j in range(lp.num_col_):
        start, end = mat.start_[j:j+2]
        rows = np.asarray(mat.index_[start:end], dtype=int)
        vals = np.asarray(mat.value_[start:end])
        np.add.at(ax, rows, vals * x[j])
        aty[j] = vals @ y[rows]
        aty_abs[j] = np.abs(vals) @ np.abs(y[rows])
    objective = float(c @ x + lp.offset_)
    sense = int(lp.sense_)
    abs_p = rel_p = dual_sign = comp = 0.0
    dual_obj = sense * lp.offset_
    for activity, multiplier, lower, upper in (
        (x, sense * rc, lp.col_lower_, lp.col_upper_),
        (ax, sense * y, lp.row_lower_, lp.row_upper_),
    ):
        lower, upper = np.asarray(lower), np.asarray(upper)
        # HiGHS returns +/-infinity for absent bounds.
        for bound, lower_side in ((lower, True), (upper, False)):
            finite = np.isfinite(bound)
            violation = np.maximum((bound - activity) if lower_side else (activity - bound), 0)[finite]
            if len(violation):
                abs_p = max(abs_p, float(np.max(violation)))
                rel_p = max(rel_p, float(np.max(violation / (1 + np.abs(bound[finite])))))
        chosen = np.where(multiplier >= 0, lower, upper)
        valid = np.isfinite(chosen)
        if (~valid).any():
            dual_sign = max(dual_sign, float(np.max(np.abs(multiplier[~valid]))))
        dual_obj += float(multiplier[valid] @ chosen[valid])
        if valid.any():
            comp = max(comp, float(np.max(np.abs(multiplier[valid] * (activity[valid] - chosen[valid])))))
    stationarity = float(np.max(np.abs(c - aty - rc) / (1 + np.abs(c) + aty_abs + np.abs(rc)), initial=0))
    gap = abs(sense * objective - dual_obj) / (1 + abs(objective))
    error = abs(objective - reference_objective) / (1 + abs(reference_objective))
    published_error = abs(objective - answer['objective']) / (1 + abs(objective))
    metrics = dict(primal_abs=abs_p, primal_rel=rel_p, dual_sign_abs=dual_sign,
                   stationarity_rel=stationarity, complementarity_abs=comp,
                   primal_dual_gap_rel=gap, reference_objective_error_rel=error,
                   published_objective_error_rel=published_error, recomputed_objective=objective)
    thresholds = dict(primal_rel=1e-7, dual_sign_abs=1e-7, stationarity_rel=1e-7,
                      primal_dual_gap_rel=1e-7, reference_objective_error_rel=1e-7,
                      published_objective_error_rel=1e-10)
    metrics['failed_checks'] = [k for k, limit in thresholds.items() if not math.isfinite(metrics[k]) or metrics[k] > limit]
    metrics['accepted'] = not metrics['failed_checks']
    metrics['absolute_primal_1e_6_pass'] = abs_p <= 1e-6
    return metrics


def reference(highspy, path, solver, presolve):
    h = highspy.Highs()
    for name, value in dict(output_flag=False, threads=1, parallel='off', solver=solver,
                            presolve=presolve, run_crossover='off', time_limit=10.0,
                            primal_feasibility_tolerance=1e-8, dual_feasibility_tolerance=1e-8,
                            ipm_optimality_tolerance=1e-8, random_seed=0).items():
        if h.setOptionValue(name, value) != highspy.HighsStatus.kOk:
            raise RuntimeError(f'rejected reference option {name}')
    read_status = h.readModel(str(path))
    if read_status == highspy.HighsStatus.kError:
        raise RuntimeError(f'cannot read reference model {path}')
    start = time.perf_counter_ns()
    h.run()
    elapsed = time.perf_counter_ns() - start
    info = h.getInfo()
    return h, dict(status=h.modelStatusToString(h.getModelStatus()), objective=h.getObjectiveValue(), read_status=str(read_status),
                   solve_ns=elapsed, ipm_iterations=info.ipm_iteration_count,
                   simplex_iterations=info.simplex_iteration_count)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--corpus', type=Path, required=True)
    p.add_argument('--runner', type=Path, required=True)
    p.add_argument('--highspy-path', type=Path, action='append', required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--backend', choices=['sparse', 'profile', 'dense'], default='sparse')
    p.add_argument('--refinement-steps', type=int, default=0)
    p.add_argument('--max-iterations', type=int, default=200)
    p.add_argument('--feasibility-tolerance', type=float, default=1e-7)
    p.add_argument('--primal-correction-steps', type=int, default=5)
    p.add_argument('--max-regularization', type=float, default=1e-4)
    p.add_argument('--equilibrate', action='store_true')
    p.add_argument('--ordering', choices=['automatic', 'natural', 'rcm', 'minimum-degree'], default='automatic')
    p.add_argument('--models', nargs='*')
    args = p.parse_args()
    sys.path[:0] = [str(path.resolve()) for path in args.highspy_path]
    import highspy
    import numpy as np
    records = []
    report = dict(highs_version=highspy.Highs().version(), python=sys.version,
                  platform=platform.platform(), backend=args.backend,
                  configuration={'zhighs_iteration_limit': args.max_iterations, 'zhighs_time_limit_s': 5,
                                 'zhighs_refinement_steps': args.refinement_steps,
                                 'zhighs_feasibility_tolerance': args.feasibility_tolerance,
                                 'zhighs_primal_correction_steps': args.primal_correction_steps,
                                 'zhighs_max_regularization': args.max_regularization,
                                 'zhighs_equilibrate': args.equilibrate,
                                 'zhighs_ordering': args.ordering,
                                 'process_timeout_s': 20, 'highs_time_limit_s': 10,
                                 'samples': 1, 'highs_threads': 1, 'highs_crossover': 'off'}, models=records)
    paths = sorted(args.corpus.glob('*.mps'))
    if args.models:
        paths = [args.corpus / (name + '.mps') for name in args.models]
    if not paths:
        raise RuntimeError('empty corpus')
    lock_path = Path(__file__).resolve().parents[1] / 'simplex' / 'end_to_end_corpus.lock.tsv'
    locked = {row.split('\t')[0]: row.split('\t')[1] for row in lock_path.read_text().splitlines() if row and not row.startswith('#')}
    for path in paths:
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        row = dict(model=path.stem, sha256=sha, frozen_corpus_hash_match=locked.get(path.stem) == sha)
        oracle, row['reference_simplex'] = reference(highspy, path, 'simplex', 'on')
        lp = oracle.getLp()
        if lp.integrality_:
            raise RuntimeError(f'{path} is not a continuous LP')
        row.update(rows=lp.num_row_, cols=lp.num_col_, nnz=len(lp.a_matrix_.value_))
        _, row['reference_ipm'] = reference(highspy, path, 'ipm', 'off')
        start = time.perf_counter_ns()
        try:
            run = subprocess.run([str(args.runner.resolve()), str(path.resolve()), args.backend,
                                  str(args.refinement_steps), str(args.max_iterations),
                                  str(args.feasibility_tolerance), str(args.primal_correction_steps),
                                  str(args.max_regularization), 'scaled' if args.equilibrate else 'unscaled',
                                  args.ordering],
                                 capture_output=True, text=True, timeout=20)
            if run.returncode:
                row['zhighs'] = dict(status='process_error', returncode=run.returncode, stderr=run.stderr[-2000:])
            else:
                answer = json.loads(run.stdout)
                if answer['status'] == 'optimal' and row['reference_simplex']['status'] == 'Optimal':
                    row['validation'] = validate(lp, answer, row['reference_simplex']['objective'], np)
                row['zhighs'] = {key: value for key, value in answer.items() if key not in ('primal', 'dual', 'reduced_cost')}
                row['parser_dimensions_match'] = all(answer[key] == row[key] for key in ('rows', 'cols', 'nnz'))
        except subprocess.TimeoutExpired:
            row['zhighs'] = dict(status='process_timeout')
        row['process_ns'] = time.perf_counter_ns() - start
        row['accepted'] = row.get('validation', {}).get('accepted', False) and row.get('parser_dimensions_match', False)
        records.append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
        print(f"{path.stem:12} {row['rows']:5}x{row['cols']:<5} {row['zhighs']['status']:20} ref={row['reference_simplex']['status']:10} accepted={row['accepted']}", flush=True)
    print(f"Validated optimal: {sum(r['accepted'] for r in records)}/{len(records)}; all results: {args.output}")


if __name__ == '__main__':
    main()
