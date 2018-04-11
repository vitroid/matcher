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
	-rm *.o *~ $(EXE)
	-rm */*.o */*.mod
	-rm *.so
	-rm -rf build

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
test:
	./smatcher 7.gro 0.8 0.06 
test2:
	python ./smatcher.py 7.gro 8.0 0.06
test3:
	./matcher -e 0.03 -v 0.06 7.gro 1c.ar3r
test4:
	python ./matcher.py 7.gro 1c.gro 0.03 0.06
test5:
	genice 7 -r 6 6 6 --dens 1.6 -f matcher[1c.gro:0.03:0.06:0] 
