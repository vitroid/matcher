.DELETE_ON_ERROR:
OS=$(shell uname)
ifeq ($(OS), Darwin)
	DEST=~/Library/Application\ Support/GenIce
else
	DEST=~/.genice
endif
LDFLAGS=-lm
EXE=matcher smatcher
%.o: %.c %.h
	$(CC) -std=c99 -c -g -O $< -o $@
#for tests
%: %.c
%: %.o pairlist.o bst.o common.o
	$(CC) -g -O $^ -o $@
all: $(EXE)
prepare:
	cp ../PairList/pairlist.[ch] .
clean:
	-rm *.o *~ $(EXE) @ @@ @@@ @@@@
	-rm */*.o */*.mod 
	-rm *.so
	-rm -rf build *.egg-info dist
	-rm *.gro *.ar3r *.match *.smatch

#python
check:
	python setup.py check
build.:
	-rm *.so
	-rm -rf build
	python setup.py build_ext --inplace
install: matcher.py
	python setup.py install
	install -d $(DEST)
	install -d $(DEST)/formats
	install matcher.py $(DEST)/formats


#end python

7.gro:
	genice 7 -r 6 6 6 -d 1.6> $@
1c.ar3r:
	genice 1c -r 1 1 1 -d 0.8 -f rcom > $@
1c.gro:
	genice 1c -r 1 1 1 -d 0.8 > $@
test: 7.gro
	./smatcher 7.gro 0.8 0.06 > test.smatch
test2: 7.gro
	python ./smatcher.py 7.gro 8.0 0.06 > test2.smatch
test3: 7.gro 1c.ar3r
	./matcher -e 0.03 -v 0.06 7.gro 1c.ar3r > test3.1c.match
test4: 7.gro 1c.gro
	python ./matcher.py 7.gro 1c.gro 0.03 0.06 > test4.1c.match
test5: 7.1c.match
7.1c.match: 1c.gro
	genice 7 -r 6 6 6 --dens 1.6 -f matcher[1c.gro:0.03:0.06:0] > 7.1c.match
%.1c.match.yap: %.match 1c.ar3r
	head $< | ./match2yap.py $*.gro 1c.ar3r > $@

