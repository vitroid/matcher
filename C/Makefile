.DELETE_ON_ERROR:
LDFLAGS=-lm
EXE=smatcher matcher2
%.o: %.c %.h
	$(CC) -std=c99 -c -g -O $< -o $@
#for tests
%: %.c
%: %.o pairlist.o bst.o common.o svd.o neighborlist.o
	$(CC) -g -O $^ -o $@ $(LDFLAGS)
all: $(EXE)
prepare:
	-chmod a+w pairlist.[ch]
	cp ../../PairList/pairlist.[ch] .
	chmod a-w pairlist.[ch]
clean:
	-rm 7.1c.match*
	-rm *.o *~
	-rm $(EXE)

7.gro:
	genice 7 -r 6 6 6 -d 1.6> $@
T2x224.gro:
	genice T2B -r 2 2 4 --add_noise 1 > $@
T2.ar3r:
	genice T2B -f rcom > $@
T2.gro:
	genice T2B > $@
1c.ar3r:
	genice 1c -r 1 1 1 -d 0.8 -f rcom > $@
1c.gro:
	genice 1c -r 1 1 1 -d 0.8 > $@


test: 7.gro
	./smatcher 7.gro 0.8 0.06 > test.smatch
test3: 7.1c.match22
7.1c.match22: 7.gro 1c.gro
	./matcher2 7.gro 1c.gro > $@
test5: 7.1c.match22.yap
%.1c.match22.yap: %.1c.match22 1c.gro
	head $< | ../matcher/matcher2yap.py $*.gro 1c.gro > $@

test6: T2x224.gro T2.gro
	./matcher2 T2x224.gro T2.gro


