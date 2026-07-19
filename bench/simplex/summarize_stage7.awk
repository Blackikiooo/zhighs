# Summarize the raw output produced by run_stage7_corpus.sh. Timing percentiles
# describe the selected model set; repeated-run benchmark percentiles remain a
# separate measurement and must not be inferred from this single-pass report.

BEGIN { FS = OFS = "\t" }

function absolute(value) { return value < 0 ? -value : value }

function sort_values(values, count,    i, j, temporary) {
  for (i = 2; i <= count; i++) {
    temporary = values[i]
    j = i - 1
    while (j >= 1 && values[j] > temporary) {
      values[j + 1] = values[j]
      j--
    }
    values[j + 1] = temporary
  }
}

$1 == "result" {
  model = $2
  solver = $3
  process = $4
  rss = $6 + 0
  selected[solver SUBSEP model] = 1
  if (rss > peak_rss[solver]) peak_rss[solver] = rss

  # A completed run can contain several zhighs statistics records. Only the
  # native solver record ($7 equals the solver name) carries status and time.
  if (process != "completed") {
    if (!process_seen[solver SUBSEP model]++) process_count[solver SUBSEP process]++
    next
  }
  if ($7 != solver) next

  status = tolower($9)
  status_count[solver SUBSEP status]++
  status_by_model[solver SUBSEP model] = status
  if (solver == "zhighs") {
    objective[solver SUBSEP model] = $11 + 0
    total_ns = $25 + 0
    primal = $13 + 0
    dual = $14 + 0
  } else {
    objective[solver SUBSEP model] = $10 + 0
    total_ns = $16 + 0
    primal = $12 + 0
    dual = $13 + 0
  }
  if (primal > max_primal[solver]) max_primal[solver] = primal
  if (dual > max_dual[solver]) max_dual[solver] = dual
  time_count[solver]++
  time_ns[solver, time_count[solver]] = total_ns
  solvers[solver] = 1
  models[model] = 1
}

END {
  print "section", "solver", "metric", "value", "detail"
  for (solver in solvers) {
    selected_count = 0
    for (key in selected) {
      split(key, part, SUBSEP)
      if (part[1] == solver) selected_count++
    }
    print "coverage", solver, "selected_models", selected_count, ""
    for (key in process_count) {
      split(key, part, SUBSEP)
      if (part[1] == solver) print "process", solver, part[2], process_count[key], ""
    }
    for (key in status_count) {
      split(key, part, SUBSEP)
      if (part[1] == solver) print "status", solver, part[2], status_count[key], ""
    }
    print "residual", solver, "max_primal", max_primal[solver] + 0, ""
    print "residual", solver, "max_dual", max_dual[solver] + 0, ""
    print "memory", solver, "peak_rss_kb", peak_rss[solver] + 0, ""

    n = time_count[solver]
    delete values
    for (i = 1; i <= n; i++) values[i] = time_ns[solver, i]
    sort_values(values, n)
    if (n != 0) {
      median_index = int((n + 1) / 2)
      p95_index = int((95 * n + 99) / 100)
      if (p95_index > n) p95_index = n
      print "time", solver, "total_ns_median", values[median_index], "models=" n
      print "time", solver, "total_ns_p95", values[p95_index], "models=" n
    }
  }

  common = matched = mismatched = 0
  max_gap = 0
  for (model in models) {
    if (status_by_model["zhighs" SUBSEP model] != "optimal" ||
        status_by_model["highs" SUBSEP model] != "optimal") continue
    common++
    z_objective = objective["zhighs" SUBSEP model]
    h_objective = objective["highs" SUBSEP model]
    gap = absolute(z_objective - h_objective) / (1 + absolute(h_objective))
    if (gap > max_gap) max_gap = gap
    if (gap <= 1e-7) matched++
    else mismatched++
  }
  print "comparison", "zhighs-highs", "common_optimal", common, ""
  print "comparison", "zhighs-highs", "objective_matches", matched, "tolerance=1e-7 relative"
  print "comparison", "zhighs-highs", "objective_mismatches", mismatched, ""
  print "comparison", "zhighs-highs", "max_relative_gap", max_gap + 0, ""
}
