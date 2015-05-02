CC=		g++
#CC=		clang --analyze
CFLAGS=		-c -Wall -Wno-unused-function -O3 -std=gnu++11
AR=		ar
VERSION=	0.010
LOBJS=
PROG=		uqct
INCLUDES=	
SUBDIRS=	.
OBJS=		b6.o
LIBS=		-lz -L.
DEBUG=		-g
OPT=		-O3
SOURCES=	gz.cpp b6.cpp seq.cpp map.cpp indexio.cpp key_init.cpp \
		fa.cpp fq.cpp main.cpp
DEFINES+=	-DPROGRAM_NAME=\"$(PROG)\" -DPROGRAM_VERSION=\"$(VERSION)\" # -DKEY_LENGTH=11
PREPROCESSED=	$(SOURCES:.cpp=.i)
ASSEMBLIES=	$(SOURCES:.cpp=.s)
OBJECTS=	$(SOURCES:.cpp=.o)
ARCHIVE=	$(PROG)_$(VERSION)

.PHONY: clean distclean dist

all:$(SOURCES) $(PROG)

prep: $(PREPROCESSED)

asm: $(ASSEMBLIES)

$(PROG):libuqct.a $(OBJECTS)
	$(CC) $(DEFINES) $(DEBUG) $(OBJECTS) -o $@ $(LIBS)

.cpp.s:
	$(CC) $(DEFINES) $(DEBUG) $(CFLAGS) $< -S -o $@ $(LIBS)

.cpp.o:
	$(CC) $(DEFINES) $(DEBUG) $(CFLAGS) $< -o $@ $(LIBS)

.cpp.i:
	$(CC) $(DEFINES) $(DEBUG) $(CFLAGS) $< -E -o $@ $(LIBS)

libuqct.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)


clean:
	rm -f $(PREPROCESSED) $(ASSEMBLIES) $(OBJECTS) $(PROG) core vgcore.*

distclean: clean
	rm -rf $(ARCHIVE).tar.gz

dist:
	tar -czf $(ARCHIVE).tar.gz $(SOURCES)  fqless *.c *.h *.pl *.sh .git Makefile

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.cpp )

gz.o: gz.h

b6.o: b6.h

seq.o: gz.h seq.h b6.h

fq.o: gz.h b6.h seq.h fq.h

fa.o: klib/khash.h gz.h b6.h seq.h fa.h

map.o: gz.h b6.h seq.h fa.h

main.o: gz.h b6.h seq.h fa.h fq.h


