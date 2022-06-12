
UMFD	= umd2fb.o umd2co.o umd2fa.o umd2f0.o umd2f1.o umd2f2.o \
	  umd2fg.o umd2in.o umd2of.o umd2s2.o umd2sl.o umd2so.o umd2su.o \
	  umd2er.o umd2p1.o umd2p2.o umd2lt.o umd2ut.o \
	  umd2rf.o umd2ra.o umd2r0.o umd2r2.o umd2rg.o

UMFS	= ums2fb.o ums2co.o ums2fa.o ums2f0.o ums2f1.o ums2f2.o \
	  ums2fg.o ums2in.o ums2of.o ums2s2.o ums2sl.o ums2so.o ums2su.o \
	  ums2er.o ums2p1.o ums2p2.o ums2lt.o ums2ut.o \
	  ums2rf.o ums2ra.o ums2r0.o ums2r2.o ums2rg.o

HARWELL	= mc21b.o mc13e.o
# BLASLIB	= libblas.a
BLASLIB	= -lblas -llapack

FFLAGS	= -O -g -ggdb
F77	= gfortran

PREFIX = /usr/local

OBJDMAIN = dmain.o  $(UMFD) $(BLASLIB) $(HARWELL)
OBJSMAIN = smain.o  $(UMFS) $(BLASLIB) $(HARWELL)

PROGRAMS = dmain smain

OUTFILES = dmain.out smain.out

UMFPACKLIB = libumfpack

all:  $(OUTFILES)

dmain.out: dmain in
	./dmain < in > dmain.out

smain.out: smain in
	./smain < in > smain.out

clean:
	- rm *.o *.a

purge: clean
	- rm $(OUTFILES) $(PROGRAMS)

dmain: $(OBJDMAIN)
	$(F77) $(FFLAGS) -o dmain $(OBJDMAIN)

smain: $(OBJSMAIN)
	$(F77) $(FFLAGS) -o smain $(OBJSMAIN)

libumfpack.a: $(UMFD) $(UMFS) $(HARWELL) $(BLASLIB)
	ar ru $@ $^
	ranlib $@

install: $(UMFPACKLIB).a
	install -d -p $(DESTDIR)$(PREFIX)/lib/
	install -m 644 -p $(UMFPACKLIB).a $(DESTDIR)$(PREFIX)/lib/
# 	install -d $(DESTDIR)$(PREFIX)/include/
# 	install -m 644 $(UMFPACKLIB).h $(DESTDIR)$(PREFIX)/include/

uninstall: $(UMFPACKLIB).a
	rm $(DESTDIR)$(PREFIX)/lib/$(UMFPACKLIB).a

.SUFFIXES:
.SUFFIXES: .o .f

.f.o:
	$(F77) $(FFLAGS) -c $<
