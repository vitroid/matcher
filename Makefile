.DELETE_ON_ERROR:
OS=$(shell uname)
ifeq ($(OS), Darwin)
	DEST=~/Library/Application\ Support/GenIce
else
	DEST=~/.genice
endif
clean:
	make -C C clean
	-rm matcher/*.so
	-rm -rf build *.egg-info dist
	-rm *.gro *.ar3r *.match *.smatch
	-rm *~ matcher/*~

#python
check:
	python setup.py check
build.:
	-rm *.so
	-rm -rf build
	python setup.py build_ext --inplace
install: matcher.py
	-make clean
	python setup.py install
	install -d $(DEST)
	install -d $(DEST)/formats
	install smatcher.py matcher.py $(DEST)/formats
pypi: check
	./setup.py sdist bdist_wheel upload


#end python

7.gro:
	genice 7 -r 6 6 6 -d 1.6> $@
T2x224.gro:
	genice T2B -r 2 2 4 --add_noise 1 > $@
T2.ar3r:
	genice T2B -f rcom > $@
T2.gro:
	genice T2B > $@
R.gro:
	genice iceR > $@
T.gro:
	genice iceT > $@
1c.ar3r:
	genice 1c -r 1 1 1 -d 0.8 -f rcom > $@
1c.gro:
	genice 1c -r 1 1 1 -d 0.8 > $@
test2: 7.gro
	python ./smatcher.py 7.gro 8.0 0.06 > test2.smatch
test4: 7.gro 1c.gro
	python ./matcher.py 7.gro 1c.gro 0.03 0.06 > test4.1c.match
test5: 7.1c.match
7.1c.match: 1c.gro
	genice 7 -r 6 6 6 --dens 1.6 -f matcher[1c.gro:0.03:0.06:0] > 7.1c.match

test7: T2x224.gro T2.gro
	python ./matcher.py T2x224.gro T2.gro 0.01 0.4 
test8: T2x224.gro T2.gro  # use with analice
	analice T2x224.gro -f matcher[T2.gro:0.01:0.4]

test10: T2x224.gro T2.gro
	python ./matcher2.py T2x224.gro T2.gro
T2x224.R.match:T2x224.gro R.gro
	python ./matcher2.py T2x224.gro R.gro > $@
%.R.match.yap: %.R.match R.gro
	./matcher2yap.py -v 0.045 $*.gro R.gro R < $< > $@
T2x224.T.match:T2x224.gro T.gro
	python ./matcher2.py T2x224.gro T.gro > $@
%.T.match.yap: %.T.match T.gro
	./matcher2yap.py -v 0.06 $*.gro T.gro T < $< > $@



lattices/T2x22.py:
	genice T2B -f reshape[1,1,0,-1,1,0,0,0,2] --add_noise 1 > $@
T2x22.gro:
	genice T2x22 > $@

