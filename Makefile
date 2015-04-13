.PHONY: all run clean

default: all

all:
	g++ antal.cpp -o antal
	g++ antalElitist.cpp -o antalElitist
	g++ antalRank.cpp -o antalRank

run:
	#python output.py
	./antal

clean:
	rm -rf antal antalElitist antalRank *.txt *.o *.pyc
