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
	-rm *.gro *.ar3r *.match2 *.smatch *.cmatch2 *.csmatch
	-rm *~ matcher/*~

#python
check:
	python setup.py check
build.:
	-rm *.so
	-rm -rf build
	python setup.py build_ext --inplace
install:
	-make clean
	python setup.py install
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

# Python
# run python module directly
pytest1: 7.smatch 7.1c.match2
%.smatch: %.gro
	python matcher/smatcher.py 7.gro 8.0 0.1 > 7.smatch
%.1c.match2: 7.gro 1c.gro
	python matcher/matcher2.py $*.gro 1c.gro > $@
# run as a genice plugin
pytest2: 7.gsmatch 7.1c.gmatch2
%.gsmatch: %.gro
	genice 7 -r 6 6 6 -f smatcher[8.0:0.1] > 7.smatch
%.1c.gmatch2: 7.gro 1c.gro
	python matcher/matcher2.py $*.gro 1c.gro > $@
# run as an analice plugin

# C and cpython
ctest1: 7.csmatch 7.1c.cmatch2
%.csmatch: %.gro
	C/smatcher 7.gro 8.0 0.1 > $@
%.1c.cmatch2: 7.gro 1c.gro
	C/matcher2 $*.gro 1c.gro > $@
# run as a genice plugin (C)
# run as an analice plugin (C)
