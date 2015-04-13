.PHONY: all tsp run clean

default: all

tsp:
	g++ TSPConstructor.cpp -o TSPConstructor
	./TSPConstructor

all: tsp
	g++ antal.cpp -o antal
	g++ antalElitist.cpp -o antalElitist
	g++ antalRank.cpp -o antalRank

run:
	python output.py

clean:
	rm -rf TSPConstructor antal antalElitist antalRank *.txt *.o *.pyc
