# Presolve ablation results

Binary: `zig-out\bin\presolve-bench.exe`; repetitions per sample: 10; samples: 3.
Times are means over the repetitions, then a median over the samples. `off` disables presolve;
`legacy` keeps the pre-existing rules; `implied` and `parallel` enable one new rule each;
`full` enables both. Objective must match `off` whenever the raw solve is optimal.

| model | mode | rows | cols | nnz | implied fix | parallel rows | presolve us | solve us | total us | iters | status |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| forcing_lp.lp | off | 60 | 100 | 120 | 0 | 0 | 0.0 | 731.5 | 731.5 | 101 | optimal |
| forcing_lp.lp | legacy | 40 | 60 | 80 | 0 | 0 | 51.6 | 660.4 | 717.2 | 61 | optimal |
| forcing_lp.lp | implied | 0 | 0 | 0 | 60 | 0 | 27.4 | 10.5 | 40.6 | 0 | optimal |
| forcing_lp.lp | parallel | 40 | 60 | 80 | 0 | 0 | 55.1 | 613.4 | 671.7 | 61 | optimal |
| forcing_lp.lp | full | 0 | 0 | 0 | 60 | 0 | 25.1 | 10.2 | 36.9 | 0 | optimal |
| forcing_lp.lp | speedup vs legacy | | | | | | | | 19.46x | | |
| parallel_lp.lp | off | 40 | 20 | 80 | 0 | 0 | 0.0 | 179.8 | 179.8 | 22 | optimal |
| parallel_lp.lp | legacy | 40 | 13 | 80 | 0 | 0 | 25.4 | 170.9 | 197.7 | 15 | optimal |
| parallel_lp.lp | implied | 40 | 13 | 80 | 0 | 0 | 27.3 | 119.0 | 146.5 | 15 | optimal |
| parallel_lp.lp | parallel | 10 | 13 | 20 | 0 | 30 | 22.0 | 33.4 | 56.8 | 15 | optimal |
| parallel_lp.lp | full | 10 | 13 | 20 | 0 | 30 | 53.5 | 62.4 | 116.9 | 15 | optimal |
| parallel_lp.lp | speedup vs legacy | | | | | | | | 1.69x | | |
| mip_forcing.lp | off | 40 | 120 | 280 | 0 | 0 | 0.0 | 562.5 | 562.5 | 94 | optimal |
| mip_forcing.lp | legacy | 40 | 120 | 280 | 0 | 0 | 28.1 | 531.6 | 560.3 | 94 | optimal |
| mip_forcing.lp | implied | 20 | 40 | 200 | 80 | 0 | 48.5 | 28.2 | 78.0 | 0 | optimal |
| mip_forcing.lp | parallel | 40 | 120 | 280 | 0 | 0 | 76.1 | 1302.8 | 1378.9 | 94 | optimal |
| mip_forcing.lp | full | 20 | 40 | 200 | 80 | 0 | 68.6 | 35.9 | 106.2 | 0 | optimal |
| mip_forcing.lp | speedup vs legacy | | | | | | | | 5.28x | | |
| branching_compare.lp | off | 10 | 16 | 60 | 0 | 0 | 0.0 | 55.4 | 55.4 | 21 | optimal |
| branching_compare.lp | legacy | 10 | 16 | 60 | 0 | 0 | 9.0 | 64.7 | 73.8 | 21 | optimal |
| branching_compare.lp | implied | 10 | 16 | 60 | 0 | 0 | 14.3 | 81.9 | 96.6 | 21 | optimal |
| branching_compare.lp | parallel | 10 | 16 | 60 | 0 | 0 | 8.1 | 58.4 | 66.7 | 21 | optimal |
| branching_compare.lp | full | 10 | 16 | 60 | 0 | 0 | 8.6 | 45.5 | 54.6 | 21 | optimal |
| branching_compare.lp | speedup vs legacy | | | | | | | | 1.35x | | |
| presolve_mip.mps | off | 2 | 3 | 4 | 0 | 0 | 0.0 | 12.7 | 12.7 | 2 | optimal |
| presolve_mip.mps | legacy | 1 | 2 | 2 | 0 | 0 | 2.9 | 7.9 | 11.0 | 2 | optimal |
| presolve_mip.mps | implied | 1 | 2 | 2 | 1 | 0 | 5.9 | 12.1 | 18.3 | 2 | optimal |
| presolve_mip.mps | parallel | 1 | 2 | 2 | 0 | 0 | 3.6 | 9.1 | 12.9 | 2 | optimal |
| presolve_mip.mps | full | 1 | 2 | 2 | 1 | 0 | 4.7 | 8.6 | 14.6 | 2 | optimal |
| presolve_mip.mps | speedup vs legacy | | | | | | | | 0.75x | | |
| warm_restart_mip.lp | off | 14 | 22 | 122 | 0 | 0 | 0.0 | 70.4 | 70.4 | 30 | optimal |
| warm_restart_mip.lp | legacy | 14 | 22 | 122 | 0 | 0 | 6.7 | 69.4 | 76.1 | 30 | optimal |
| warm_restart_mip.lp | implied | 14 | 22 | 122 | 0 | 0 | 10.8 | 70.3 | 81.1 | 30 | optimal |
| warm_restart_mip.lp | parallel | 14 | 22 | 122 | 0 | 0 | 8.3 | 68.2 | 76.5 | 30 | optimal |
| warm_restart_mip.lp | full | 14 | 22 | 122 | 0 | 0 | 13.1 | 74.0 | 87.3 | 30 | optimal |
| warm_restart_mip.lp | speedup vs legacy | | | | | | | | 0.87x | | |
| pipeline_lp.lp | off | 1 | 2 | 2 | 0 | 0 | 0.0 | 11.6 | 11.6 | 2 | optimal |
| pipeline_lp.lp | legacy | 0 | 0 | 0 | 0 | 0 | 4.2 | 8.4 | 13.1 | 0 | optimal |
| pipeline_lp.lp | implied | 0 | 0 | 0 | 1 | 0 | 3.3 | 5.8 | 9.4 | 0 | optimal |
| pipeline_lp.lp | parallel | 0 | 0 | 0 | 0 | 0 | 3.2 | 6.8 | 10.2 | 0 | optimal |
| pipeline_lp.lp | full | 0 | 0 | 0 | 1 | 0 | 4.3 | 7.3 | 11.9 | 0 | optimal |
| pipeline_lp.lp | speedup vs legacy | | | | | | | | 1.10x | | |
