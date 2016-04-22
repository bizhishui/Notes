% Notes on using LINUX
% Jinming LYU
% \today

- Commands to view the contents of library.
  1. For static library, use `ar -t lib*.a`
  2. For dynamic library, we can use `nm -D --defined-only libGreenCPU.so` or `objdump -T libGreenCPU.so | grep text`. The first 
     one will just give all the symbol names defined in this file, the second one will only the symbols in the text section.

- How to create a bootable USB stick on Ubuntu
  1. Insert a USB stick with at least 2GB of free space
  2. Open the dash and search for **Startup Disk Creator**
  3. Select the Startup Disk Creator to launch the app
  4. Click *Other* to choose the downloaded ISO file if it isn’t found automatically, select the file and click *Open*
  5. Select the USB stick in the bottom box and click *Make Startup Disk* and then *Yes*
  6. That's it! When the process completes, you'll be ready to restart your computer and begin installing Ubuntu
