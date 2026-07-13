# isssues

## todolists

- [] There are too many deep copy everywhere, such as compile_model, so I need to figure out the lifetime of
  these variables, it's better to realise zero-copy.
- [] So many struct type may cause waste in memory, and also lower the performance, need to tackle.
