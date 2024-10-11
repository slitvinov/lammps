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
arg_info.o \
atom_map.o \
atom.o \
atom_vec_atomic.o \
atom_vec.o \
comm_brick.o \
comm.o \
comm_tiled.o \
compute.o \
compute_pe_atom.o \
compute_pressure.o \
create_atoms.o \
create_box.o \
domain.o \
error.o \
fix_nve.o \
fix.o \
fmtlib_format.o \
fmtlib_os.o \
force.o \
group.o \
hashlittle.o \
input.o \
integrate.o \
label_map.o \
lammps.o \
lattice.o \
main.o \
math_eigen.o \
math_extra.o \
math_special.o \
memory.o \
modify.o \
my_page.o \
my_pool_chunk.o \
nbin.o \
nbin_standard.o \
neighbor.o \
neigh_list.o \
neigh_request.o \
npair_half_bin_atomonly_newton.o \
npair.o \
pair_dpd.o \
pair.o \
platform.o \
procmap.o \
random_mars.o \
random_park.o \
rcb.o \
region_block.o \
region.o \
run.o \
text_file_reader.o \
tokenizer.o \
universe.o \
update.o \
utils.o \
variable.o \
verlet.o \

main: $(S:.cpp=.o)
	$(MPICXX) -o main $(S:.cpp=.o) $(LDFLAGS) -ldl
.cpp.o:
	$(MPICXX) -o $@ -c $< $(FLAGS) $(CXXFLAGS)
clean:
	-rm main $(S:.cpp=.o)
