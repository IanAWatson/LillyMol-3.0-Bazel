LillyMol Software Enhanced and built with Bazel
===============================================

This was cloned from 

https://github.com/EliLillyCo/LillyMol

This version offers some interesting changes/additions

1. Uses Bazel to build the system.
2. Ths substructure query specification language now uses
   a protocol buffer, so is fully described.
3. The Google test framework is used to create a wide
   variety of tests for the substructure searching code.
4. There are implementations of circular fingerprints,
   linear fingerprints and atom pair fingerprints. These
   are fairly rudimentary, but should be functional enough
   to be useful.

To read a proto query with tsubstructure use

`-q PROTO:file.proto`

To read a file that contains a listing of proto query files try

`-q PROTOFILE:file.txt`

Installation
------------
To get this going with bazel I would recommend the following:

Use the protocol buffer compiler that comes with your distro. I
tried to download a newer version, but it proved very difficult
to make that work.

Similarly, if you can get bazel with your distro, use that.

Download Google Test and Google Mock from 

[GoogleTest](https://github.com/google/googletest)

Edit WORKSPACE here and designate the location where you have installed
Google test/mock on your system. *You will need to change this.*

The code does require c++17, so to build the system try things like

bazel build --cxxopt='-std=c++17' -c opt Molecule_Tools:fileconv
or
bazel build --cxxopt='-std=c++17' -c opt Molecule_Tools:all

A recent debugging session saw this command being useful

/usr/bin/bazel run --compilation_mode=dbg --config asan --cxxopt='-std=c++17' Molecule_Tools:ec_fingerprint -- -J NCEC -D file=$PWD/precedent.txt -D puse -P UST:ACHPY -v -i do=10000 -n -s -l -g all -E autocreate $PWD/bad.smi

Using Bazel outside Google can be quite challenging, it took me several
iterations to get a workable configuration, and I don't claim what I have
is optimal. Particularly problematic was to have multiple versions isntalled.

Much of this new code is less tested than the official release. I am hoping
that the advantages of testing and having a well defined input will ultimately
result in merging.

If Bazel proves too difficult for some, it would not be hard to create Makefiles.
But I am convinced that Bazel is a fundamentally better build system than make.
