# Heuristic for Virtual Network Embedding with 1 + 1 Network Protection

## System Requirement
The heuristic is implemented as a multi-threaded program. For better performance
run the heuristic on a machine with atleast 4 to 8 CPU cores. For best results
run the program on a machine with 32 or more CPU cores.

## How to run
```
$ make
$ ./vne_heuristic  --pn_topology_file=<physical_network_topology>\ 
                   --vn_topology_file=<virtual_network_topology>\
                   --location_constraint_file=<location_constraint_file>
```

Two example physical (test_pn.topo) and virtual (test_vn.topo) network topology
files are provided with the distribution. A sample location constraint file is
provided as well (test_location.txt).

For more extensive testing please refer to the inputs in TestSet-0 and TestSet-1 
directories. 

TestSet-0 directory contains a 10 test cases, each within its own 
directory (test0 to test9). A test case directory contains the following files:
  * sn.txt = Specification of physical network
  * vn.txt = Specification of virtual network request
  * vnloc.txt = Location constraint for the virtual network request

TestSet-1 has a bit different directory structure and the directory structure is
as follows (repeated and unnecessary file/directory names are not shown):
```
TestSet-1/
└── case0
   ├── sn2.txt
   └── vnr
       ├── random
       │   ├── test0
       │   │   ├── vnloc.txt
       │   │   ├── vn.txt
       │   ├── test1
       │   └── test2
       ├── ring
       │   ├── test0
       │   ├── test1
       │   └── test2
       └── star
           ├── test0
           ├── test1
           └── test2
```
The difference here is that the virtual network requests are grouped according 
to the type of the request (ring, start, random). The semantics of the files
are as follows:
  * sn2.txt = Specification of physical network
  * vn.txt = Specification of a virtual network request
  * vnloc.txt = Location constraint for the virtual network request

Extensive testing with these test sets can be automated by running the provided
run_experiments.py script. Run the script as follows:
```
$ python run_experiments.py --testcase_root <a_testset_directory>\
                            --executable <name_of_the_executable>
```
For example, the following command will run the heuristic for all the test cases
located under TestSet-1 directory:
```
$ python run_experiments.py --testset_root TestSet-1 --executable vne_heuristic
```
Please refer to the  "Input file format" section for details on the format of 
the input files.

## Input file format

A topology file contains the list of edges. Each line contains a description of
an edge in a comma separated value (CSV) format. Format of a line is as follows:
```
<LinkID>,<SourceNodeId>,<DestinationNodeId>,<PeerID>,<Cost>,<Bandwidth>,<Latency>
```
Where,
  * LinkID = Id of a link. Ignored for both physical and virtual topology.
  * SourceNodeId = 0-based node index of the source of the link*
  * DestinationNodeId = 0-based node index of the destination of the link*
  * PeerID = Ignored
  * Cost = Cost of provisioning unit bandwidth on this link. Cost is ignored for
           virtual links.
  * Bandwidth = Available bandwidth of a physical link. In case of virtual link,
                this is the bandwidth requirement
  * Delay = Latency of a physical link. In case of virtual link, this is the
            maximum delay requirement for the virtual link.

A location constraint file contains as many lines as the number of virtual
nodes. Each line is a comma separated list of values. The first value indicates
the id of a virtual node followed by the ids of physical nodes where this
virtual node can be mapped.

*Nodes are numberded from `0 ... (n - 1)` in a network with `n` nodes.

## Output Files

Currently the program writes the output to different files. Each output file is
prefixed with the virtual topology file name and has the following suffixes
based on the contents:

* .cost = solution cost
* .nmap =primary node mapping
* .emap = primary link mapping
* .snmap = shadow/backup node mapping
* .semap = shadow/backup link mapping
