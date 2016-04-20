% Notes on using LINUX
% Jinming LYU
% \today

- Commands to view the contents of library.
  1. For static library, use `ar -t lib*.a`
  2. For dynamic library, we can use `nm -D --defined-only libGreenCPU.so` or `objdump -T libGreenCPU.so | grep text`. The first 
     one will just give all the symbol names defined in this file, the second one will only the symbols in the text section.
