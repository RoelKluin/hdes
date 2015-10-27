HOST=$(shell hostname)
#There must be a space between ifeq and (
ifeq	($(HOST),Z)
  CC=		ccache g++
else
  CC=		/opt/rh/devtoolset-3/root/usr/bin/x86_64-redhat-linux-g++
endif
##CCC=		colorgcc
#CC=		clang++-3.5 --analyze
#-fprofile-arcs -ftest-coverage
CFLAGS=		-c -Wall -Wno-unused-function -O3 -std=gnu++11 -fdiagnostics-color=always -I./${EXTERNAL_ZLIB}
#CXXFLAGS += --analyze -Xanalyzer -analyzer-output=text
AR=		ar
VERSION=	0.016
LOBJS=
PROG=		uqct
INCLUDES=	
SUBDIRS=	.
OBJS=		b6.o
EXTERNAL_ZLIB=zlib-1.2.8/
LIBS=		-L. -L./zlib-1.2.8/
DEBUG=		-g #-pg --coverage -rdynamic
OPT=		-O3
SOURCES=	gz.cpp b6.cpp seq.cpp mantra.cpp map.cpp indexio.cpp key_init.cpp \
		fa.cpp fq.cpp main.cpp
DEFINES+=	-DPROGRAM_NAME=\"$(PROG)\" -DPROGRAM_VERSION=\"$(VERSION)\" # -DKEY_LENGTH=11
PREPROCESSED=	$(SOURCES:.cpp=.i)
ASSEMBLIES=	$(SOURCES:.cpp=.s)
OBJECTS=	$(SOURCES:.cpp=.o)
ARCHIVE=	$(PROG)_$(VERSION)
GCOV=		gcov


.PHONY: clean distclean dist
all:$(SOURCES) $(PROG)

prep: $(PREPROCESSED)

asm: $(ASSEMBLIES)

$(PROG):libuqct.a $(OBJECTS)
	$(CC) $(DEFINES) $(DEBUG) $(OBJECTS) ${EXTERNAL_ZLIB}libz.a -o $@ $(LIBS)

.cpp.s:
	$(CC) $(DEFINES) $(DEBUG) $(CFLAGS) $(CXXFLAGS) $< -S -o $@ $(LIBS)

.cpp.o:
	$(CC) $(DEFINES) $(DEBUG) $(CFLAGS) $(CXXFLAGS) $< -o $@ $(LIBS)

.cpp.i:
	$(CC) $(DEFINES) $(DEBUG) $(CFLAGS) $(CXXFLAGS) $< -E -o $@ $(LIBS)

libuqct.a:$(OBJS)
		$(AR) -csr $@ $(OBJS)

coverage:
	$(GCOV) -b $(SOURCES)
clean:
	rm -f $(PREPROCESSED) $(ASSEMBLIES) $(OBJECTS) $(PROG) core vgcore.*

distclean: clean
	rm -rf $(ARCHIVE).tar.gz

dist:
	tar -czf $(ARCHIVE).tar.gz $(SOURCES) faless fqless *.cpp *.c *.h *.pl *.sh .git Makefile doc/hdes.odt documentation.txt PLAN IDEE BUGS dodev

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) $(DEF) -- *.cpp )

gz.o: gz.h

b6.o: b6.h

seq.o: gz.h seq.h b6.h

fq.o: gz.h b6.h seq.h fq.h

fa.o: klib/khash.h gz.h b6.h seq.h fa.h

mantra.o: fa.h

map.o: gz.h b6.h seq.h fa.h

main.o: gz.h b6.h seq.h fa.h fq.h

.PHONY: zlib
external: ${EXTERNAL_ZLIB}libz.a

${EXTERNAL_ZLIB}libz.a:
	cd ${EXTERNAL_ZLIB} && ./configure && ${MAKE} libz.a


