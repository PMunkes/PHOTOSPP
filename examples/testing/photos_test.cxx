/**
 * Main program for testing photos C++ interface.
 * Pythia events are generated first and photos used for FSR.
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

//PHOTOS header files
#include "Photos/Photos.h"
#include "Photos/PhotosHepMCEvent.h"
#include "Photos/Log.h"

using namespace std;
using namespace Pythia8;
using namespace Photospp;

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

// Finds X Y -> 6 -6 decay and converts it to 100 -> 6 -6, where 100 = X + Y
// Method used only in test for t bar  t pair production
void fixForMctester(HepMC::GenEvent *evt)
{
	for(HepMC::GenEvent::particle_const_iterator p=evt->particles_begin();p!=evt->particles_end(); p++)
	if((*p)->pdg_id()==6)
	{
		HepMC::GenParticle *pt = *p;
		int id=(* pt->production_vertex()->particles_in_const_begin() )->pdg_id();
		if(id!=21 && id!=11 && id>5) continue;

		// Get first mother and add 2x second mother to it
		HepMC::GenParticle *X = (* pt->production_vertex()->particles_in_const_begin());
		HepMC::GenParticle *Y = (* ++(pt->production_vertex()->particles_in_const_begin()) );
		HepMC::FourVector fX = X->momentum();
		HepMC::FourVector fY = Y->momentum();
		HepMC::FourVector fXY(fX.px()+fY.px(),fX.py()+fY.py(),fX.pz()+fY.pz(),fX.e()+fY.e());
		X->set_momentum(fXY);
		// Unique ID for MC-Tester to analyze
		X->set_pdg_id(100);

		// Set 2nd mother as decayed and delete it from production vertex
		Y->set_status(1);
		(* Y->production_vertex()->particles_in_const_begin())->set_status(1);
		pt->production_vertex()->remove_particle(Y);
		return;
	}
}

int main(int argc,char **argv)
{

	// Program needs at least 2 parameters
	if(argc<3)
	{
		cout<<endl<<"Usage: "<<argv[0]<<" <pythia_conf> <no_events> [ <special_mode> ]"<<endl;
		cout<<endl<<"   eg. "<<argv[0]<<" pythia_W.conf 0 10000 4 0"<<endl;
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
	  3. Special mode - default(off), ttbar, NLO
	  4. Photos - use 1-photon mode on/off

	  Example where all input parameters are used:

	  ./photos_test.exe pythia_W.conf 100000 0 0
	    - use pythia_W.conf
	    - generate 100 000 events
	    - default configuration (not using any special mode)
	    - Photos is not using 1-photon mode (default option, except for WmunuNLO and ZmumuNLO)
	*********************************************************/

	// 1. Load pythia configuration file (argv[1], from console)
	pythia.readFile(argv[1]);

	// 2. Get number of events (argv[2], from console)
	NumberOfEvents = atoi(argv[2]);

    pythia.init();
	Photos::initialize();

	Photos::setInfraredCutOff(1.e-6);
	Photos::maxWtInterference(3.0);

	bool topDecays    = false;
	bool pairEmission = false;

	// 3. Check if we're using any special mode
	if(argc>3)
	{
		// Top decays
		if(atoi(argv[3])==1)      topDecays=true;
		// NLO mode
		else if(atoi(argv[3])==2)
		{
			Photos::setMeCorrectionWtForW(true);
			Photos::setMeCorrectionWtForZ(true);
			//Photos::meCorrectionWtForScalar(true);
		}
		// Pairs emission
		else if(atoi(argv[3])==4)
		{
			pairEmission = true;
		}
	}

	// 4. Check if we're using 1-photon mode
	if(argc>4 && atoi(argv[4])==1)
	{
		Photos::setDoubleBrem(false);
		Photos::setExponentiation(false);
		Photos::setInfraredCutOff(0.001);
		Photos::maxWtInterference(2.0);
	}

	// Photon emission is turned on by default
	// Use this flag to turn it off if needed
	//Photos::setPhotonEmission(false);
	Photos::setPairEmission(pairEmission);

	MC_Initialize();

	Photos::iniInfo();
	Log::SummaryAtExit();
	cout.setf(ios::fixed);

	// Begin event loop
	for(unsigned long iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<iEvent*(100./NumberOfEvents)<<"%)"<<endl;

		// --- Event no transmitted inside Photos:: 
		Photos::setEventNo(iEvent);
		// --- may be useful for temporary event specific test prints.

		if (!pythia.next()) continue;

		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(event, HepMCEvt);

		if(iEvent<EventsToCheck)
		{
			cout<<"                                          "<<endl;
			cout<<"Momentum conservation check BEFORE/AFTER Photos"<<endl;
			checkMomentumConservationInEvent(HepMCEvt);
		}

		// Run PHOTOS on the event
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();

		if(iEvent<EventsToCheck)
		{
			checkMomentumConservationInEvent(HepMCEvt);
		}

		// Top decays - we mess with the event so MC-TESTER can work on it as in LC analysis case
		if(topDecays) fixForMctester(HepMCEvt);

		// Run MC-TESTER on the event
		HepMCEvent temp_event(*HepMCEvt,false);
		MC_Analyze(&temp_event);

		// Clean up
		delete HepMCEvt;
	}
	pythia.stat();
	MC_Finalize();

	// Additional printout for pairs
	if( pairEmission ) {
	// PAIR emission
	// Test with formula 11 from UTHEP-93-0301 M. Skrzypek ...
	const double PI=3.141592653589793238462643;     
	const double ALFINV= 137.01;

	double deno=log(91/2/0.000511);
	deno=deno*deno*(deno+log(2.))*4;
        double delta=5;//.125;//0.125; //0.25;
        double L=log(2*delta/0.000511)-5.0/6.0;
        double num=0.99* 4.0/3.0*(L*L*L/3.+ (31./36.- PI*PI/6.)*L+0.5940);
	printf ("    >>>    Soft pairs emission --- probability tests    <<< \n");
	printf (" Delta = %15.8f GeV (set the same in pairs.cxx): \n", delta);
	printf (" Z->ee: \n");
	printf (" Abslolute= %15.8f  Relative to crude= %15.8f \n",num/PI/PI/ALFINV/ALFINV/0.99, num/deno);


       	deno=log(91/2/0.1056);
	deno=deno*deno*(deno+log(2.))*4;
        L=log(2*delta/0.1056)-5.0/6.0;
        num=0.99* 4.0/3.0*(L*L*L/3.+ (31./36.- PI*PI/6.)*L+0.5940);
	printf (" Z->mumu: \n");
	printf (" Abslolute= %15.8f  Relative to crude= %15.8f \n",num/PI/PI/ALFINV/ALFINV/0.99, num/deno);
	}
	
}
