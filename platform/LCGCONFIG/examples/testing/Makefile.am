exampledir = $(top_srcdir)/examples/testing

INCLUDES = -I$(prefix)/include

AM_LDFLAGS = -R $(prefix)/lib

LDADD =                          \
$(FLIBS)                         \
-L$(prefix)/lib/ -lPhotospp -lPhotosppHepMC

example_PROGRAMS =

### All examples require HepMC, Pythia8, and MC-TESTER ###
if HAS_HEPMC
if HAS_PYTHIA8
if HAS_MCTESTER

AM_LDFLAGS += -R $(HEPMC_DIR)/lib    \
              -R $(PYTHIA8_DIR)/lib  \
              -R $(MCTESTER_DIR)/lib
              
INCLUDES += -I$(HEPMC_DIR)/include    \
            -I$(PYTHIA8_DIR)/include  \
            -I$(MCTESTER_DIR)/include \
            -I$(ROOTINC)

LDADD    += -L$(HEPMC_DIR)/lib -lHepMC                               \
            -L$(PYTHIA8_DIR)/lib/ -lpythia8                          \
            -L$(MCTESTER_DIR)/lib -lHEPEvent -lHepMCEvent -lMCTester \
            $(ROOTLIBS)

example_PROGRAMS += photos_test.exe

photos_test_exe_SOURCES  = photos_test.cxx

if HAS_TAUOLA

AM_LDFLAGS += -R $(TAUOLA_DIR)/lib
INCLUDES   += -I$(TAUOLA_DIR)/include
LDADD      += -L$(TAUOLA_DIR)/lib -lTauolaFortran -lTauolaCxxInterface

example_PROGRAMS += photos_tauola_test.exe

photos_tauola_test_exe_SOURCES = photos_tauola_test.cxx

endif # TAUOLA

endif # MC-TESTER
endif # PYTHIA8
endif # HEPMC
