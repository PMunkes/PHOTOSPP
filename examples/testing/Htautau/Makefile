############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.benchmark.Htautau.root


CONF=pythia_H.conf
NAME=Htautau

PROGRAM=$(PHOTOS_TAUOLA)
TAUOLA_MODE=$(RHORHO_MODE)

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) $(EVENTS) $(TAUOLA_MODE)
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root 
