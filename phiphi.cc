#include <boost/program_options.hpp>
#include <string>

#include "Pythia8/Pythia.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"

int main(int argc, char **argv)
{

  int nevents, seed;
  std::string config, output;
  bool verbose, onechannel;
  
  /** process arguments **/
  namespace po = boost::program_options;
  po::options_description desc("Options");
  try {
    desc.add_options()
      ("help", "Print help messages")
      ("nevents,n"  , po::value<int>(&nevents)->default_value(10), "Number of events to be generated")
      ("config,c"   , po::value<std::string>(&config), "Configuration file")
      ("seed,s"     , po::value<int>(&seed)->default_value(0), "Random seed")
      ("output,o"   , po::value<std::string>(&output)->default_value("phiphi.root"), "Output ROOT file")
      ("verbose,V"  , po::bool_switch(&verbose)->default_value(false), "Verbose event listing")
      ("onechannel" , po::bool_switch(&onechannel)->default_value(false), "Force phi(1020) -> K+ K-")
      ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return 1;
    }
  }
  catch(std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }
  
  // pythia, config and init
  Pythia8::Pythia pythia;
  Pythia8::Rndm   rndm;
  // baseline config
  std::vector<std::string> _cmd =
    {
      "Beams:eCM = 14000.",
      "SoftQCD:nonDiffractive = on",
      "ParticleDecays:limitTau0 = on",
      "ParticleDecays:tau0Max = 10",
      "333:m0         = 1.01946",
      "333:mWidth     = 0.00426",
      "333:mMin       = 0.99000",
      "333:mMax       = 1.05000",
      "Random:setSeed = on",
      "Random:seed = " + std::to_string(seed),
      "Next:numberShowEvent = 0"
    };
  for (auto cmd : _cmd)
    pythia.readString(cmd);
  if (onechannel) 
    pythia.readString("333:oneChannel = 1 1. 0 -321  321");
  // config from file
  if (!config.empty() && !pythia.readFile(config)) {
    std::cout << "Error: could not read config file \"" << config << "\"" << std::endl;
    return 1;
  }
  pythia.init();
  rndm.init();

  // histograms and output
  auto fout = TFile::Open(output.c_str(), "RECREATE");
  // K+K-
  short n_ka_plusminus;
  float ka_plusminus_m[1024], ka_plusminus_pt[1024], ka_plusminus_y[1024];
  short ka_plusminus_i[1024], ka_plusminus_j[2014];
  bool  ka_plusminus_iseff[1024], ka_plusminus_isphi[1024];
  auto tout_ka_plusminus = new TTree("ka_plusminus", "Output tree");
  tout_ka_plusminus->Branch("n",     &n_ka_plusminus,     "n/S");
  tout_ka_plusminus->Branch("m",     &ka_plusminus_m,     "m[n]/F");
  tout_ka_plusminus->Branch("pt",    &ka_plusminus_pt,    "pt[n]/F");
  tout_ka_plusminus->Branch("y",     &ka_plusminus_y,     "y[n]/F");
  tout_ka_plusminus->Branch("i",     &ka_plusminus_i,     "i[n]/S");
  tout_ka_plusminus->Branch("j",     &ka_plusminus_j,     "j[n]/S");
  tout_ka_plusminus->Branch("iseff", &ka_plusminus_iseff, "iseff[n]/O");
  tout_ka_plusminus->Branch("isphi", &ka_plusminus_isphi, "isphi[n]/O");
  // K+K+
  short n_ka_plusplus;
  float ka_plusplus_m[1024], ka_plusplus_pt[1024], ka_plusplus_y[1024];
  bool  ka_plusplus_iseff[1024], ka_plusplus_isphi[1024];
  auto tout_ka_plusplus = new TTree("ka_plusplus", "Output tree");
  tout_ka_plusplus->Branch("n",     &n_ka_plusplus,     "n/I");
  tout_ka_plusplus->Branch("m",     &ka_plusplus_m,     "m[n]/F");
  tout_ka_plusplus->Branch("pt",    &ka_plusplus_pt,    "pt[n]/F");
  tout_ka_plusplus->Branch("y",     &ka_plusplus_y,     "y[n]/F");
  tout_ka_plusplus->Branch("iseff", &ka_plusplus_iseff, "iseff[n]/O");
  tout_ka_plusplus->Branch("isphi", &ka_plusplus_isphi, "isphi[n]/O");
  // K-K-
  short n_ka_minusminus;
  float ka_minusminus_m[1024], ka_minusminus_pt[1024], ka_minusminus_y[1024];
  bool  ka_minusminus_iseff[1024], ka_minusminus_isphi[1024];
  auto tout_ka_minusminus = new TTree("ka_minusminus", "Output tree");
  tout_ka_minusminus->Branch("n",     &n_ka_minusminus,     "n/I");
  tout_ka_minusminus->Branch("m",     &ka_minusminus_m,     "m[n]/F");
  tout_ka_minusminus->Branch("pt",    &ka_minusminus_pt,    "pt[n]/F");
  tout_ka_minusminus->Branch("y",     &ka_minusminus_y,     "y[n]/F");
  tout_ka_minusminus->Branch("iseff", &ka_minusminus_iseff, "iseff[n]/O");
  tout_ka_minusminus->Branch("isphi", &ka_minusminus_isphi, "isphi[n]/O");
  
  // event loop
  for (int iev = 0; iev < nevents; ++iev) {

    // generate event
    pythia.next();
    if (verbose) pythia.event.list();

    std::vector<Pythia8::Particle> _ka_plus, _ka_minus;
    
    // particle loop
    for (int ipa = 0; ipa < pythia.event.size(); ++ipa) {

      const auto particle = pythia.event[ipa];
      if (!particle.isFinal()) continue;
      if (particle.id() ==  321) _ka_plus.push_back(particle);
      if (particle.id() == -321) _ka_minus.push_back(particle);
      
    }  // end of particle loop

    // K+ K- pairs
    n_ka_plusminus = 0;
    for (int iplus = 0; iplus < _ka_plus.size(); ++iplus) {
      const auto &ka_plus = _ka_plus.at(iplus);
      for (int iminus = 0; iminus < _ka_minus.size(); ++iminus) {
	const auto &ka_minus = _ka_minus.at(iminus);
	auto P  = ka_plus.p() + ka_minus.p();
	auto m  = P.mCalc();
	auto pt = P.pT();
	auto y  = P.rap();
	auto iseff = fabs(ka_plus.eta()) < 0.8 && fabs(ka_minus.eta()) < 0.8 && ka_plus.pT() > 0.15 && ka_minus.pT() > 0.15;
	auto isphi = ka_plus.mother1() == ka_minus.mother1() && pythia.event[ka_plus.mother1()].id() == 333;
	
	if (m < 0.99 || m > 1.05) continue;
	if (fabs(y) > 1.) continue;
	
	ka_plusminus_m[n_ka_plusminus] = m;
	ka_plusminus_pt[n_ka_plusminus] = pt;
	ka_plusminus_y[n_ka_plusminus] = y;
	ka_plusminus_i[n_ka_plusminus] = iplus;
	ka_plusminus_j[n_ka_plusminus] = iminus;
	ka_plusminus_iseff[n_ka_plusminus] = iseff;
	ka_plusminus_isphi[n_ka_plusminus] = isphi;
	n_ka_plusminus++;

      }	
    }
    tout_ka_plusminus->Fill();
    
    continue;

    // K+ K+ pairs
    n_ka_plusplus = 0;
    for (int iplus = 0; iplus < _ka_plus.size(); ++iplus) {
      const auto &ka_plus_1 = _ka_plus.at(iplus);
      for (int iiplus = iplus + 1; iiplus < _ka_plus.size(); ++iiplus) {
	const auto &ka_plus_2 = _ka_plus.at(iiplus);
	auto P  = ka_plus_1.p() + ka_plus_2.p();
	auto m  = P.mCalc();
	auto pt = P.pT();
	auto y  = P.rap();
	auto iseff = fabs(ka_plus_1.eta()) < 0.8 && fabs(ka_plus_2.eta()) < 0.8 && ka_plus_1.pT() > 0.15 && ka_plus_2.pT() > 0.15;
	auto isphi = ka_plus_1.mother1() == ka_plus_2.mother1() && pythia.event[ka_plus_1.mother1()].id() == 333;
	
	if (m < 0.99 || m > 1.05) continue;
	if (fabs(y) > 1.) continue;
	
	ka_plusplus_m[n_ka_plusplus] = m;
	ka_plusplus_pt[n_ka_plusplus] = pt;
	ka_plusplus_y[n_ka_plusplus] = y;
	ka_plusplus_iseff[n_ka_plusplus] = iseff;
	ka_plusplus_isphi[n_ka_plusplus] = isphi;
	n_ka_plusplus++;

      }	
    }
    tout_ka_plusplus->Fill();

    // K- K- pairs
    n_ka_minusminus = 0;
    for (int iminus = 0; iminus < _ka_minus.size(); ++iminus) {
      const auto &ka_minus_1 = _ka_minus.at(iminus);
      for (int iiminus = iminus + 1; iiminus < _ka_minus.size(); ++iiminus) {
	const auto &ka_minus_2 = _ka_minus.at(iiminus);
	auto P  = ka_minus_1.p() + ka_minus_2.p();
	auto m  = P.mCalc();
	auto pt = P.pT();
	auto y  = P.rap();
	auto iseff = fabs(ka_minus_1.eta()) < 0.8 && fabs(ka_minus_2.eta()) < 0.8 && ka_minus_1.pT() > 0.15 && ka_minus_2.pT() > 0.15;
	auto isphi = ka_minus_1.mother1() == ka_minus_2.mother1() && pythia.event[ka_minus_1.mother1()].id() == 333;
	
	if (m < 0.99 || m > 1.05) continue;
	if (fabs(y) > 1.) continue;
	
	ka_minusminus_m[n_ka_minusminus] = m;
	ka_minusminus_pt[n_ka_minusminus] = pt;
	ka_minusminus_y[n_ka_minusminus] = y;
	ka_minusminus_iseff[n_ka_minusminus] = iseff;
	ka_minusminus_isphi[n_ka_minusminus] = isphi;
	n_ka_minusminus++;

      }	
    }
    tout_ka_minusminus->Fill();

  } // end of event loop

  // write output and close
  fout->cd();
  tout_ka_plusminus->Write();
  //  tout_ka_plusplus->Write();
  //  tout_ka_minusminus->Write();
  fout->Close();

  return 0;
}
