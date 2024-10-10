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
balance.o \
comm_brick.o \
comm.o \
comm_tiled.o \
compute.o \
compute_pe_atom.o \
compute_pe.o \
compute_pressure.o \
compute_temp.o \
create_atoms.o \
create_box.o \
domain.o \
error.o \
fix_abf.o \
fix_fileforce.o \
fix_nve.o \
fix.o \
fix_pair.o \
fix_rbc.o \
fix_read_restart.o \
fix_sdf_bounceback.o \
fix_sdf_repforce.o \
fmtlib_format.o \
fmtlib_os.o \
force.o \
group.o \
hashlittle.o \
imbalance_group.o \
imbalance_neigh.o \
imbalance.o \
imbalance_store.o \
imbalance_var.o \
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
nbin_multi.o \
nbin.o \
nbin_standard.o \
neighbor.o \
neigh_list.o \
neigh_request.o \
npair_copy.o \
npair_full_bin_atomonly.o \
npair_full_bin.o \
npair_full_multi.o \
npair_half_bin_atomonly_newton.o \
npair_half_bin_newtoff.o \
npair_half_bin_newton.o \
npair_half_bin_newton_tri.o \
npair_halffull_newtoff.o \
npair_halffull_newtoff_trim.o \
npair_halffull_newton.o \
npair_halffull_newton_trim.o \
npair_half_multi_newtoff.o \
npair_half_multi_newton.o \
npair_half_multi_newton_tri.o \
npair_half_size_bin_newtoff.o \
npair_half_size_bin_newton.o \
npair_half_size_bin_newton_tri.o \
npair_half_size_multi_newtoff.o \
npair_half_size_multi_newton.o \
npair_half_size_multi_newton_tri.o \
npair.o \
pair_dpd.o \
pair.o \
platform.o \
procmap.o \
random_mars.o \
random_park.o \
rcb.o \
read_restart.o\
region_block.o \
region.o \
region_prism.o \
run.o \
set.o \
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
