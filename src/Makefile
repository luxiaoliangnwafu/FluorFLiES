FC = gfortran

FFLAGS1 = -funroll-all-loops -O2
#FFLAGS1 = -g -fbounds-check -Wall
#FFLAGS1 = -g -fbacktrace -ffpe-trap=zero,overflow,underflow -Wall
TARGET1 = flies

OBJ=main.o iparam.o imath.o math.o rparam.o idivspace.o wdata.o frndi.o frnd.o igtbl.o canort.o mccnp.o mcstm.o mcflr.o r_acos.o elpss.o helpss.o cyls.o cones.o planes.o vegrroulette.o transform.o srfref.o scatvec.o vegtrace.o vegrad.o pf.o ipf.o mc1descape.o arthshift.o vegfeye.o getphs.o prepatm.o mc1doptics.o mc1dtrace.o getext.o artphsfintp.o twovecang.o phsftrunc5.o artsanglut.o nrmlzphsf.o artscat.o getrancircle.o artrroulette.o phsfbiasym2.o i_rvctrssrch1.o ifunc_vctrbinsrch.o nrmlzuvctr.o clai.o fluorweight.o fluorrad.o pretracing.o

#all : clean $(TARGET1)

$(TARGET1) : $(OBJ)
	$(FC) $(FFLAGS1)  $(OBJ) -o $@

.f.o :
	$(FC) -c $(FFLAGS1)  $< -o $@

clean :
	rm -f ./*.o

install :
	mv flies ../

distclean :
	rm -f ./*.o $(TARGET1)
