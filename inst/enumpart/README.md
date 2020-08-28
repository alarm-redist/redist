Enumpart
=================

# Program that computes the ZDD representing all graph partitions of a given graph using the TdZdd library

This program computes the ZDD representing all graph partitions of a given graph.
This program uses [the TdZdd library](https://github.com/kunisura/TdZdd), mainly developed by Dr. Hiroaki Iwashita.
The file ```enumpart.cpp``` is based on [graphillion_example in TdZdd library](https://github.com/kunisura/TdZdd/tree/master/apps/graphillion) also written by Dr. Hiroaki Iwashita.

## Usage

### Build

```
make
```

### Run

```
./enumpart <graph_file> [<weight_file>] [some options]
```

### example

The following command constructs the ZDD representing all graph partitions of the 2x2 grid graph.

```
./enumpart grid2x2.dat
```

The following command constructs the ZDD representing all graph partitions of the 2x2 grid graph such that the number of components is 2.

```
./enumpart grid2x2.dat -k 2
```

The following command constructs the ZDD representing all graph partitions of the 2x2 grid graph with weight indicated by ```grid2x2_v.dat``` such that the weight of each component is between 3 and 8.

```
./enumpart grid2x2.dat grid2x2_v.dat -k 2 -lower 3 -upper 8
```

The following command constructs the ZDD representing all graph partitions of the 2x2 grid graph with weight indicated by ```grid2x2_v.dat``` such that the number of components is 2 and the ratio of the maximum weight to the minimum weight is at most 2.0.

```
./enumpart grid2x2.dat grid2x2_v.dat -k 2 -ratio 2.0
```

## File format

### graph_file

Example:

```
1 2
1 3
2 3
2 4
```

The format of <graph_file> is a so-called edge list. One line corresponds to an edge
and consists of two integers representing the vertex numbers of the endpoints.
In the above example, there are four vertices 1, 2, 3 and 4, and four edges
connecting 1 with 2, 1 with 3, 2 with 3 and 2 with 4.

### weight_file

Example:

```
1 2 4 8
```

The weight_file specifies weights of vertices.
Each number is separated by a white space.

### output file

If you specify the "-allsols" option, you obtain all the
solutions.

For example,

```
./enumpart grid2x2.dat -k 2 -allsols
```

outputs the following texts:

```
3 4
2 4
2 3
1 4
1 3
1 2
```

A line corresponds to a solution. The numbers in a line represent
edge numbers in the corresponding solution. For example,
the solution corresponding to the first line "3 4" consists of
the 3rd edge "2 3" and the 4th edge "2 4" in "grid2x2.dat".
As a result, the graph includes an isolated component (vertex 1)
and the other component (consists of vertices 2,3 and 4).

This command and its output indicate that the number of ways
how to divide the 2x2 grid graph into two components is six.

The "-comp" option together with "-allsols" outputs solutions in the
"component" format as follows.

Command:

```
./enumpart grid2x2.dat -k 2 -allsols -comp
```

Output:

```
0 1 1 1
0 1 0 0
0 1 0 1
0 0 1 1
0 0 1 0
0 0 0 1
```

A line corresponds to a solution. The i-th number in a line represents
which component the i-th vertex belongs to. For example,
"0 1 1 1" means that vertex 1 belongs to component 0 and
vertices 2,3,4 belong to component 1. If you use "-comp" option,
the input graph has vertex numbers 1,2,...,n (n is the number of vertices).

## Random sampling

If you specify the "-sample" (with an integer) option, random sampling
from all the solutions is carried out and the specified number of solutions
are output.

The following command outputs three random solutions.

```
./enumpart grid2x2.dat -k 2 -sample 3
```

## Storing the constructed ZDD into a file

We can store the constructed ZDD into a file, and then
by using the stored ZDD we can enumerate/sample solutions.
By adding the "-export" option, enumpart outputs the ZDD
to the standard output.

The following command outputs the ZDD to "zdd.txt" file.

```
./enumpart grid2x2.dat -export >zdd.txt
```

The "sample" program reads the ZDD file and enumerate/sample solutions.

The following command reads "zdd.txt" file and output all solutions.

```
./sample grid2x2.dat zdd.txt -allsols
```

The "sample" program supports the "count", "allsols", "solutions <n>",
"comp", "drawsol <n>", and "sample <n>" options. Their usages are
the same as those of the enumpart program.

## Link

* [TdZdd](https://github.com/kunisura/TdZdd/)
* Jun Kawahara, Takashi Horiyama, Keisuke Hotta, and Shin-ichi Minato,
  "Generating All Patterns of Graph Partitions within a Disparity Bound,"
  In Proceedings of the 11th International Conference and Workshops on Algorithms and Computation (WALCOM 2017), vol. 10167, pp. 119--131, 2017, [doi:10.1007/978-3-319-53925-6_10](https://dx.doi.org/10.1007/978-3-319-53925-6_10).
