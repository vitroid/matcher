# CC=c99
LDFLAGS=-lm
EXE=matcher slide-matcher
%.o: %.c
	$(CC) -c -g -O $< -o $@
#for tests
%: %.c
%: %.o pairlist.o bst.o
	$(CC) -g -O $^ -o $@
all: $(EXE)
prepare:
	cp ../PairList/pairlist.c .
	sed -e s/double/float/g ../PairList/pairlist.h > pairlist.h
clean:
	-rm *.o *~ $(EXE)
	-rm */*.o */*.mod




%.rect.gro: %.gro
	python3 rectify.py < $< > $@

TPPI.ar3r:
	genice TPPI --rep 1 1 1 --format r > TPPI.ar3r
PPI.ar3r:
	genice PPI --rep 1 1 1 --format r --dens 1.65 > PPI.ar3r
6.ar3r:
	genice 6 --rep 1 1 1 --format r --dens 1.6 > 6.ar3r
7.ar3r:
	genice 7 --rep 1 1 1 --format r --dens 1.55 > 7.ar3r
12.ar3r:
	genice 12 --rep 1 1 1 --format r --dens 1.6 > 12.ar3r
%.TPPI.match: %.gro TPPI.ar3r # Makefile
	./matcher $< TPPI.ar3r 0.03 0.60 > $@
%.TPPI.match.yap: %.TPPI.match TPPI.ar3r
	python3 ./match2yap.py $*.gro TPPI.ar3r T < $< > $@
%.PPI.match: %.gro PPI.ar3r
	./matcher $< PPI.ar3r 0.05 0.30 > $@
%.PPI.match.yap: %.PPI.match PPI.ar3r
	python3 ./match2yap.py $*.gro PPI.ar3r R < $< > $@
%.6.match: %.gro 6.ar3r
	./matcher $< 6.ar3r 0.05 0.40 > $@
%.6.match.yap: %.6.match 6.ar3r
	python3 ./match2yap.py $*.gro 6.ar3r T < $< > $@
%.7.match: %.gro 7.ar3r 
	./matcher $< 7.ar3r 0.05 0.20 > $@
%.7.match.yap: %.7.match 7.ar3r
	python3 ./match2yap.py -e 0 $*.gro 7.ar3r C < $< > $@
%.12.match: %.gro 12.ar3r
	./matcher $< 12.ar3r 0.03 0.40 > $@
%.12.match.yap: %.12.match 12.ar3r
	python3 ./match2yap.py $*.gro 12.ar3r T < $< > $@
%.7+TPPI.yap: %.7.match.yap %.TPPI.match.yap
	perl ../Yaplot/utils/picklayer.pl 3 < $*.7.match.yap > $@
	perl ../Yaplot/utils/picklayer.pl 1 < $*.TPPI.match.yap >> $@
	echo "" >> $@  #newpage
pub:
	for ss in init 1140ps 4740ps; do  echo $$ss.rect.gro;done | xargs make -k -j 8
	for ss in init 1140ps 4740ps; do for ice in 6 7 PPI TPPI; do echo $$ss.rect.$$ice.match;done;done | xargs make -k -j 3
	for ss in init 1140ps 4740ps; do for ice in 6 7 PPI TPPI; do echo $$ss.rect.$$ice.match.yap;done;done | xargs make -k -j 8
	for ss in init 1140ps 4740ps; do  echo $$ss.rect.7+TPPI.yap;done | xargs make -k -j 8
push:
	rsync -av *.c *.h *.ar3r *py 192.168.3.3:/r3/matto/MovieWorkArea/iceT/n_ice7_6e04_300K
	rsync -av *.c *.h *.ar3r *py 192.168.3.3:/r3/matto/MovieWorkArea/iceT/ice6_6e04_300K
	rsync -av *.c *.h *.ar3r *py 192.168.3.3:/r3/matto/MovieWorkArea/iceT/PPI_6e04_300K
	rsync -av *.c *.h *.ar3r *py 192.168.3.3:/r3/matto/MovieWorkArea/iceT/TPPI_6e04_300K
pull:
	rsync -av --include="*/" --include="[^0]*.yap" --exclude="*" 192.168.3.3:/r3/matto/MovieWorkArea/iceT/n_ice7_6e04_300K .
	rsync -av --include="*/" --include="[^0]*.yap" --exclude="*" 192.168.3.3:/r3/matto/MovieWorkArea/iceT/ice6_6e04_300K .
	rsync -av --include="*/" --include="[^0]*.yap" --exclude="*" 192.168.3.3:/r3/matto/MovieWorkArea/iceT/PPI_6e04_300K .
	rsync -av --include="*/" --include="[^0]*.yap" --exclude="*" 192.168.3.3:/r3/matto/MovieWorkArea/iceT/TPPI_6e04_300K .
mov:
	ffmpeg -r 5 -i 'yaplot00_%05d.png' -vf "pad=width=512:height=1024:x=0:y=0:color=white" -pix_fmt yuv420p -crf 24 cat.mp4
