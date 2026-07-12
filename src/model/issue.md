1. [x] objects里面了坐承载了太多的功能，需要vars,models等全部拆分出来，单独实现，避免成为单文件实现；
2. [x] attrs在.zig里面的实现要变得更加zig-like；
3. attrs.fromName函数性能不够高，需要重新优化；
4. model内部的数据太分散需要使用SoA结构重构；