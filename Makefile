.POSIX:
.SUFFIXES:
.SUFFIXES: .cpp
.SUFFIXES: .o

MPICXX = mpicxx
FLAGS = \
-DLAMMPS_BIGBIG \
-DMPICH_SKIP_MPICXX \
-DOMPI_SKIP_MPICXX \

S = \
atom.o \
atom_vec_atomic.o \
atom_vec.o \
comm_brick.o \
comm.o \
create_atoms.o \
create_box.o \
domain.o \
fix_nve.o \
fix.o \
force.o \
group.o \
input.o \
integrate.o \
main.o \
memory.o \
my_page.o \
nbin.o \
nbin_standard.o \
neighbor.o \
neigh_list.o \
neigh_request.o \
npair_half_bin_atomonly_newton.o \
npair.o \
pair_dpd.o \
pair.o \
procmap.o \
random_mars.o \
random_park.o \
region_block.o \
region.o \
run.o \
universe.o \
update.o \
utils.o \
verlet.o \

main: $(S:.cpp=.o)
	$(MPICXX) -o main $(S:.cpp=.o) $(LDFLAGS)
.cpp.o:
	$(MPICXX) -o $@ -c $< $(FLAGS) $(CXXFLAGS)
clean:
	rm -f main $(S:.cpp=.o)
