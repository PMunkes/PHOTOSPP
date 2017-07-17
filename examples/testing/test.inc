EVENTS = 1000000

VIEWER=gv
COMPARE=../compare.sh

# Programs
PHOTOS_ONLY   =../photos_test.exe
PHOTOS_TAUOLA =../photos_tauola_test.exe

# Pythia modes (list of examples where they are used):
# ee@200GeV    - Wenu, Wmunu, WmunuNLO
# ee@91.187GeV - Zmumu, ZmumuNLO, Ztautau, ScalarNLO, pairs
# pp@14TeV     - Zee, ttbar
# ee@500GeV    - Htautau

# Makefiles of subdirs will choose from the lists of options:

# Special modes: 
TOP_MODE     =1 # Turn on MC-Tester fix for Top decays
NLO_MODE     =2 # Turn on Photos NLO corrections
PAIRS_MODE   =4 # Turn on pairs emission

# Tauola generation modes:
PIPI_MODE    =3 # ee@91.187GeV (for Ztautau)
RHORHO_MODE  =4 # ee@500GeV    (for Htautau)
K-STAR_MODE  =7 # ee@91.187GeV (for ScalarNLO)
