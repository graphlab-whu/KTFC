# (k,t,f)-core Project


## File Description
1. ComputeF.cpp: Source code of t-frequency computation algorithm
2. buildCFI-advanced.cpp: Source code of CF-index Construction algorithm
3. buildCFI-baseline.cpp: Source code of baseline CF-index Construction algorithm
4. query-baseline.cpp: Source code of online query algorithm
5. query-index.cpp: Source code of CF-index query algorithm
6. skyline-baseline.cpp  Source code of baseline skyline query algorithm
7. skyline-idnex.cpp  Source code of index skyline query algorithm
## Documentation

Experiments are conducted on a on a Windows machine equipped with an Intel Core i7 CPU at 2.30GHz and 32GB of RAM. Each algorithm is implemented in C++ and compiled using the g++ compiler with -O3 optimization.


The command to compile the source file and generate an executable file is as follows:

The command to construct the ComputeF is:  ```.\[computeF exe file] [tn] [qt]```, where ```[computeF exe file]``` is the executable file generated from computeF.cpp, ```[tn]``` is the number of timestamps, and  ```[qt]``` is the query parameter of t.

The command to construct the buildCFI-advanced is:  ```.\[buildCFI-advanced exe file] [graph] [choice]```, where ```[buildCFI-advanced exe file]``` is the executable file generated from buildCFI-advanced.cpp, ```[graph]``` is the graph file, and  ```[choice]``` is the choice of how to precalculate to get E_t.

The command to construct the buildCFI-baseline is:  ```.\[buildCFI-baseline exe file] [graph] [choice]```, where ```[buildCFI-baseline exe file]``` is the executable file generated from buildCFI-baseline.cpp, ```[graph]``` is the graph file, and  ```[choice]``` is the choice of how to precalculate to get f.

The command to construct the query-baseline is:  ```.\[query-baseline exe file] [graph] [qk][qt][qf][choice]```, where ```[query-baseline exe file]``` is the executable file generated from query-baseline.cpp, ```[graph]``` is the graph file, ```[qk][qt][qf]``` are three query parameters, and ```[choice]``` is the choice of how to calculate t-frequency.

The command to construct the query-index is:  ```.\[query-index exe file] [graph] [qk][qt][qf]```, where ```[query-index exe file]``` is the executable file generated from query-index.cpp, ```[graph]``` is the graph file, and ```[qk][qt][qf]``` are three query parameters.

The command to construct the skyline-baseline is:  ```.\[skyline-baseline exe file] [graph]```, where ```[skyline-baseline exe file]``` is the executable file generated from skyline-baseline.cpp, ```[graph]``` is the graph file.

The command to construct the  skyline-idnex is:  ```.\[skyline-idnex exe file] [graph]```, where ```[skyline-idnex exe file]``` is the executable file generated from  skyline-idnex.cpp, ```[graph]``` is the graph file.


The dataset used in the experiment can be downloaded from Google Drive:
**link**: https://drive.google.com/drive/folders/1cPxmTVoa7zM2ymGOUA132brQiLzCuxhU?usp=drive_link
The format of the dataset storing the graph is as follows: each line contains three integers u, v, and t, representing a connection between node u and node v at time t.

## Contact
If you have any questions, contact us by sending an email to clock@whu.edu.cn / zhongfandu@whu.edu.cn