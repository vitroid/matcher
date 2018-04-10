# CC=c99
LDFLAGS=-lm
EXE=matcher smatcher
%.o: %.c
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


#end python

7.gro:
	genice 7 -r 6 6 6 > $@
test:
	./smatcher 7.gro 0.8 0.6 
test2:
	python ./smatcher.py 7.gro 8.0 0.6
