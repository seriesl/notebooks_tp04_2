FC = gfortran
LD = gfortran

FCFLAGS = -w
LDFLAGS = 

LIBS =

OBJDIR = ./obj

EXE = run/mr_nagumo

OBJS = obj/mod_common.o \
       obj/mod_structure.o \
       obj/mod_recherche.o \
       obj/mod_initial.o \
       obj/mod_arbre.o \
       obj/mod_integration.o \
       obj/mod_reaction.o \
       obj/mod_diffusion.o \
       obj/lect_don.o \
       obj/rock4.o \
       obj/rho.o \
       obj/radau5.o \
       obj/util.o \
       obj/lapackc.o \
       obj/lapack.o \
       obj/dc_lapack.o \
       obj/main_mr.o 
       

all: $(OBJDIR) $(EXE)

$(OBJDIR):
	@echo "Creation du repertoire $@"
	@mkdir -p $@
	@echo

$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)

obj/%.o : src/%.f90
	$(FC) $(FCFLAGS) -J ./obj -c $< -o $@

obj/%.o : src/%.f
	$(FC) $(FCFLAGS) -c $< -o $@

obj/mod_structure.o : src/mod_structure.f90 \
                      obj/mod_common.o

obj/mod_recherche.o : src/mod_recherche.f90 \
                      src/recherche_support.f90 \
                      src/recherche_support_diff.f90 \
                      src/cal_indice.f90 \
                      src/creation.f90 \
                      src/recherche.f90 \
                      obj/mod_common.o \
                      obj/mod_structure.o

obj/mod_initial.o : src/mod_initial.f90 \
                    src/chainage.f90 \
                    src/support_racine.f90 \
                    src/init_sol.f90 \
                    src/init_flags.f90 \
                    src/multimesh.f90 \
                    obj/mod_common.o \
                    obj/mod_structure.o \
                    obj/mod_arbre.o

obj/mod_arbre.o : src/mod_arbre.f90 \
                  src/encoding_moy.f90 \
                  src/encoding_det.f90 \
                  src/encoding_moy_diff.f90 \
                  src/seuillage.f90 \
                  src/graduation.f90 \
                  src/graduation_local.f90 \
                  src/elagage.f90 \
                  src/chainage_support.f90 \
                  src/liste_feuille.f90 \
                  src/feuille_fictive.f90 \
                  src/valeur_feuille_fictive.f90 \
                  src/decoding_moy.f90 \
                  src/built_fine_sol.f90 \
                  src/compute_error.f90 \
                  src/liberation.f90 \
                  src/print_leaves.f90 \
                  obj/mod_common.o \
                  obj/mod_structure.o \
                  obj/mod_recherche.o

obj/mod_integration.o : src/mod_integration.f90 \
                        src/integration_reaction.f90 \
                        src/integration_diffusion.f90 \
                        src/construction_phi.f90 \
                        obj/mod_common.o \
                        obj/mod_structure.o \
                        obj/mod_recherche.o \
                        obj/mod_diffusion.o \
                        obj/mod_reaction.o \
                        obj/mod_arbre.o

obj/mod_integration.o : src/mod_reaction.f90 \
                        src/reaceq.f90 \
                        obj/mod_common.o 

obj/mod_diffusion.o : src/mod_reaction.f90 \
                      src/diffeq.f90 \
                      src/flux_diff_1.f90 \
                      src/flux_diff_2.f90 \
                      src/fonction_diffusion.f90

clean:
	rm -f obj/*.o obj/*.mod run/mr_nagumo
