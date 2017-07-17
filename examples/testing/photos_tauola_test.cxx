/**
 * Main program for testing photos C++ interface.
 * Pythia events are generated first, Tauola++ used for tau decays
 * and photos used for FSR.
 *
 * @author Nadia Davidson and Tomasz Przedzinski
 * @date 10 May 2011
 */

//Pythia header files
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

//MC-TESTER header files
#include "Generate.h"
#include "HepMCEvent.H"
#include "Setup.H"

//TAUOLA header files
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"

//PHOTOS header files
#include "Photos/Photos.h"
#include "Photos/PhotosHepMCParticle.h"
#include "Photos/PhotosHepMCEvent.h"
#include "Photos/Log.h"

using namespace std;
using namespace Pythia8;
using namespace Photospp;
using namespace Tauolapp;

unsigned long NumberOfEvents = 10000;
unsigned int EventsToCheck=20;

// elementary test of HepMC typically executed before
// detector simulation based on http://home.fnal.gov/~mrenna/HCPSS/HCPSShepmc.html
// similar test was performed in Fortran
// we perform it before and after Photos (for the first several events)
void checkMomentumConservationInEvent(HepMC::GenEvent *evt)
{
	//cout<<"List of stable particles: "<<endl;

	double px=0.0,py=0.0,pz=0.0,e=0.0;
	
	for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();
	      p != evt->particles_end(); ++p )
	{
		if( (*p)->status() == 1 )
		{
			HepMC::FourVector m = (*p)->momentum();
			px+=m.px();
			py+=m.py();
			pz+=m.pz();
			e +=m.e();
			//(*p)->print();
		}
	}
  cout.precision(6);
  cout.setf(ios_base::floatfield);
	cout<<endl<<"Vector Sum: "<<px<<" "<<py<<" "<<pz<<" "<<e<<endl;
}

int main(int argc,char **argv)
{

	// Program needs at least 3 parameters
	if(argc<4)
	{
		cout<<endl<<"Usage: "<<argv[0]<<" <pythia_conf>  <no_events> <tauola_mode> [ <alpha_order> <ScalarNLO_mode> ]"<<endl;
		cout<<endl<<"   eg. "<<argv[0]<<" pythia_H.conf 10000 4 0 0"<<endl;
		cout<<endl;
		return -1;
	}

	HepMC::Pythia8ToHepMC ToHepMC;

	// Initialization of pythia
	Pythia pythia;
	Event& event = pythia.event;

	/********************************************************
	  Read input parameters from console. List of parameters:
	  1. Pythia configuration filename
	  2. Number of events
	  3. Tauola decay mode (refer to Tauola documentation)
	  4. Photos - use 1-photon mode on/off
	  5. Photos - use ScalarNLO mode on/off

	  Example where all input parameters are used:

	  ./photos_tauola_test.exe pythia_H.conf 100000 4 0 0
	  - use pythia_H.conf
	  - generate 100 000 events
	  - fix TAUOLA decay to channel 4 (RHORHO_MODE)
	  - Photos is not using 1-photon mode (default option)
	  - Photos is not in ScalarNLO mode (default option)
	*********************************************************/

	// 1. Load pythia configuration file (argv[1], from console)
	pythia.readFile(argv[1]);

	// 2. Get number of events (argv[2], from console)
	NumberOfEvents = atoi(argv[2]);

	// 3. Set Tauola decay mode (argv[3], from console)
	// argv[3]=3 (tau => pi nu_tau)    for Ztautau
	// argv[3]=4 (tau => pi pi nu_tau) for Htautau
	Tauola::setSameParticleDecayMode(atoi(argv[3]));
	Tauola::setOppositeParticleDecayMode(atoi(argv[3]));

    pythia.init();
	Tauola::initialize();
	Photos::initialize();

	Photos::setExponentiation(true);
	Photos::setInfraredCutOff(1.e-6);
	Photos::maxWtInterference(3.0);

	// 4. Check if we're using 1-photon mode
	if( argc>4 && atoi(argv[4]) )
	{
		Photos::setDoubleBrem(false);
		Photos::setExponentiation(false);

		// Set infrared cutoff to 10MeV for scale M_Z=91.187GeV or 500 GeV
		if(atoi(argv[2])==1) Photos::setInfraredCutOff(0.01/91.187);
		else                 Photos::setInfraredCutOff(0.01/500.);

		Photos::maxWtInterference(2.0);
	}

	// 5. Check if we're in ScalarNLO mode
	if( argc>5 )
	{
		Tauola::setEtaK0sPi(1,1,0);
    
		// Check if we are using NLO
		if(atoi(argv[5])) Photos::setMeCorrectionWtForScalar(true);
	}

	Log::SummaryAtExit();
	cout.setf(ios::fixed);

	MC_Initialize();

	// Begin event loop
	for(unsigned long iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<iEvent*(100./NumberOfEvents)<<"%)"<<endl;
		if(!pythia.next()) continue;

		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(event, HepMCEvt);

		if(iEvent<EventsToCheck)
		{
			cout<<"                                          "<<endl;
			cout<<"Momentum conservation chceck BEFORE/AFTER Photos"<<endl;
			checkMomentumConservationInEvent(HepMCEvt);
		}

		// Run TAUOLA on the event
		TauolaHepMCEvent * t_event = new TauolaHepMCEvent(HepMCEvt);

		// Since we let Pythia decay taus, we have to undecay them first.
		t_event->undecayTaus();
		t_event->decayTaus();
		delete t_event;

		// Run PHOTOS on the event
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();

		if(iEvent<EventsToCheck)
		{
			checkMomentumConservationInEvent(HepMCEvt);
		}

		// Run MC-TESTER on the event
		HepMCEvent temp_event(*HepMCEvt,false);
		MC_Analyze(&temp_event);

		//clean up
		delete HepMCEvt;
	}
	pythia.stat();
	MC_Finalize();
}
