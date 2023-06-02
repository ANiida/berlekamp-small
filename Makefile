#
# Makefile in GoppaDecorder/src/Berlekamp/ by rubato6809
#
OFILES   = Berlekamp.o chash.o debug.o lu.o sha3.o fy.o inv_mat.o vc3000.o
COPTIONS = -mtune=native -march=native -ffast-math -funroll-loops -fopenmp
DEBUGOPT = -Wall -g -pg
OPTIMIZE = -O3

CFLAGS   = $(COPTIONS) $(DEBUGOPT) $(OPTIMIZE)
LDFLAGS  = $(COPTIONS) $(DEBUGOPT)
CC       = gcc # or gcc-12, otherwise clang

all: berlekamp
	touch make.date

# USEDTOBE:
#  gcc-12 -Wall -g -pg -O3
#  -mtune=native -march=native -ffast-math -funroll-loops -fopenmp
#  Berlekamp.c debug.c lu.c sha3.c fy.c inv_mat.c 

berlekamp: $(OFILES)
	$(CC) -o $@ $(LDFLAGS) $(OFILES)

det: det.o debug.o
	$(CC) -o $@ -Wall -g -pg -fopenmp det.o debug.o chash.o

# ヘッダーファイルの依存関係
Berlekamp.o : global.h struct.h gf.h val.h
inv_mat.o   : global.h struct.h gf.h
lu.o        : global.h struct.h
chash.o     : global.h struct.h 
fy.o        : global.h fy.h
sha3.o      : sha3.h
det.o       : global.h struct.h gf.h

# テスト実行
check: berlekamp
	ulimit -s unlimited
	./berlekamp

# 不要なファイルを削除
clean:
	rm -f *~ *.o gmon.out

deisclean:
	rm -f *~ *.o gmon.out make.date berlekamp det

# global.h にパラメータ
# val.h は Berlekamp.c に含めてよいのでは
