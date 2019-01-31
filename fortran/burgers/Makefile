FC = gfortran
LD = gfortran

FCFLAGS = -w
LDFLAGS = 

LIBS =

OBJDIR = ./obj

EXE = run/mr_burgers

OBJS = obj/mod_common.o \
       obj/mod_structure.o \
       obj/mod_recherche.o \
       obj/mod_fonction.o \
       obj/mod_arbre.o \
       obj/mod_initial.o \
       obj/mod_print.o \
       obj/mod_euler.o \
       obj/cell_avg.o \
       obj/mr_burgers1D.o 
       

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
                      obj/mod_common.o \
                      obj/mod_structure.o

obj/mod_fonction.o : src/mod_fonction.f90 \
                     src/liberation.f90 \
                     src/viscosite.f90 \
                     obj/mod_common.o \
                     obj/mod_structure.o

obj/mod_arbre.o : src/mod_arbre.f90 \
                  src/encoding_moy.f90 \
                  src/encoding_det.f90 \
                  src/seuillage.f90 \
                  src/graduation.f90 \
                  src/graduation_local.f90 \
                  src/elagage.f90 \
                  src/chainage_support.f90 \
                  src/liberation.f90 \
                  src/liste_feuille.f90 \
                  src/feuille_fictive.f90 \
                  src/feuille_fictive_visc.f90 \
                  src/valeur_feuille_fictive.f90 \
                  src/decoding_moy.f90 \
                  obj/mod_common.o \
                  obj/mod_structure.o \
                  obj/mod_fonction.o \
                  obj/mod_recherche.o

obj/mod_initial.o : src/mod_initial.f90 \
                    src/lect_don.f90 \
                    src/chainage.f90 \
                    src/support_racine.f90 \
                    src/init_sol.f90 \
                    src/lect_reprise.f90 \
                    src/recherche.f90 \
                    src/init_flags.f90 \
                    src/multimesh.f90 \
                    obj/mod_common.o \
                    obj/mod_structure.o \
                    obj/mod_arbre.o

obj/mod_print.o : src/mod_print.f90 \
                  src/print_leaves.f90 \
                  src/print_reprise.f90 \
                  obj/mod_common.o \
                  obj/mod_structure.o \
                  obj/mod_recherche.o \
                  obj/mod_fonction.o 

obj/mod_euler.o : src/mod_euler.f90 \
                  src/pas_de_temps.f90 \
                  src/undt_update.f90 \
                  src/u0_update.f90 \
                  src/flux_euler.f90 \
                  src/flux_osmp7.f90 \
                  src/flux_visc.f90 \
                  src/osmp7.f90 \
                  src/integration_euler.f90 \
                  src/integration1_visqueux.f90 \
                  src/integration2_visqueux.f90 \
                  src/mise_a_jour.f90 \
                  src/residu.f90 \
                  obj/mod_common.o \
                  obj/mod_structure.o \
                  obj/mod_recherche.o \
                  obj/mod_fonction.o 

clean:
	rm -f obj/*.o obj/*.mod run/mr_burgers
