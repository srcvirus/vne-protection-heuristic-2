FILES = vne_heuristic.cc util.cc 
LIB = -lpthread

all:
	g++ -O2 -std=c++0x $(FILES) $(LIB) -o vne_heuristic

dbg:
	g++ -DDBG -g -std=c++0x $(FILES) $(LIB) -o vne_heuristic

debug:
	g++ -g -std=c++0x $(FILES) $(LIB) -o vne_heuristic
