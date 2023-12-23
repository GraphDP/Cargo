# Cargo
This is a source code of the following paper:
"CARGO: Crypto-Assisted Differentially Private Triangle Counting without Trusted Servers"

# Directory Structure
* cpp/  C++ code
* data/  four real-world graphs from SNAP
* README.md  This file

# Usage
$ cd cpp/  
$ g++ -o cargo cargo.cpp  
$ ./run_Cargo.sh wiki.txt

# Execution Environment
We use AMD EPYC 7313P 16-Core Processor, 512GB RAM running Ubuntu 20.04.5 LTS with g++ 9.4.0
