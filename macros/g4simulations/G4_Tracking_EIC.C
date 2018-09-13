#include <vector>

#include "G4_GEM_EIC.C"
#include "G4_Svtx_maps_ladders+intt_ladders+tpc_KalmanPatRec_EIC.C"
//#include "G4_Svtx_maps_ladders+intt_ladders+tpc_KalmanPatRec_EIC2.C"

void TrackingInit(int verbosity = 0)
{
  /* electron-going side detectors */
  EGEM_Init();

  /* hadron-going side detectors */
  FGEM_Init();

  /* central detectors */
  SvtxInit();
}

double Tracking(PHG4Reco* g4Reco, double radius,
            const int absorberactive = 0,
            int verbosity = 0)
{
  /* Place central tracking detectors */
  Svtx(g4Reco, radius);

  /* Place electron-going side tracking detectors */
  EGEMSetup(g4Reco);

  /* Place hadron-going side tracking detectors */
  FGEMSetup(g4Reco);

  return;
}

void Tracking_Reco(int verbosity = 0, TString fileout="auto.root")
{
  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4hough.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer* se = Fun4AllServer::instance();

  //---------------------
  // Kalman Filter
  //---------------------

  PHG4TrackFastSim* kalman = new PHG4TrackFastSim("PHG4TrackFastSim",fileout.Data());
  kalman->Verbosity(0);

  kalman->VertexIn();
  kalman->VertexSmear(50E-4,50E-4); // 50 micros
  std::string phg4hits_names[16] = {"G4HIT_MAPS",   "G4HIT_SVTX",
				    "G4HIT_EGEM_0", "G4HIT_EGEM_1",
				    "G4HIT_EGEM_2O","G4HIT_EGEM_2I",
				    "G4HIT_EGEM_3O","G4HIT_EGEM_3I",
				    "G4HIT_FGEM_0", "G4HIT_FGEM_1",
				    "G4HIT_FGEM_2O","G4HIT_FGEM_2I",
				    "G4HIT_FGEM_3O","G4HIT_FGEM_3I",
				    "G4HIT_FGEM_4O","G4HIT_FGEM_4I"};
  int dettypes[16] = {
    PHG4TrackFastSim::Cylinder, PHG4TrackFastSim::Cylinder, // SVTX TPC | phi, lon
    PHG4TrackFastSim::Vertical_Plane,PHG4TrackFastSim::Vertical_Plane, // EGEM01 | rad, phi
    PHG4TrackFastSim::Vertical_Plane,PHG4TrackFastSim::Vertical_Plane, // EGEM2  | rad, phi
    PHG4TrackFastSim::Vertical_Plane,PHG4TrackFastSim::Vertical_Plane, // EGEM3  | rad, phi
    PHG4TrackFastSim::Vertical_Plane,PHG4TrackFastSim::Vertical_Plane, // FGEM01 | rad, phi
    PHG4TrackFastSim::Vertical_Plane,PHG4TrackFastSim::Vertical_Plane, // FGEM2  | rad, phi
    PHG4TrackFastSim::Vertical_Plane,PHG4TrackFastSim::Vertical_Plane, // FGEM3  | rad, phi
    PHG4TrackFastSim::Vertical_Plane,PHG4TrackFastSim::Vertical_Plane};// FGEM4  | rad, phi
  float rad[16] = {5.0,  10000,10000,10000,10000,10000,10000,10000,
		   10000,10000,10000,10000,10000,10000,10000,10000};
  float phi[16] = {5.0,  150,  50,   50,  100,   50,  100,    50,
		   50,   50,  100,   50,  100,   50,  100,    50};
  float lon[16] = {5.0,  500,  5,    5,    5,    5,    5,    5,
		   5,    5,    5,    5,    5,    5,    5,    5,};
  for(int i=0; i!=16; ++i) { // from um to cm
    rad[i] *= 1e-4;
    phi[i] *= 1e-4;
    lon[i] *= 1e-4;
  }
  
  float eff[16] = {1.0,1.0, 
		   1.0,1.0,1.0,1.0,1.0,1.0,
		   1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  float noi[16] = {0.0,0.0,
		   0.0,0.0,0.0,0.0,0.0,0.0,
		   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  kalman->Configure(phg4hits_names, dettypes, rad, phi, lon, eff, noi, 16);
  se->registerSubsystem(kalman);

  return;
}

void Tracking_Eval(std::string outputfile, int verbosity = 0)
{
  return;
  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libg4hough.so");
  gSystem->Load("libg4eval.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer* se = Fun4AllServer::instance();

  //----------------
  // SVTX evaluation
  //----------------

  SvtxEvaluator* eval;
  eval = new SvtxEvaluator("SVTXEVALUATOR", outputfile.c_str());
  eval->do_cluster_eval(false);
  eval->do_g4hit_eval(false);
  eval->do_hit_eval(false);  // enable to see the hits that includes the chamber physics...
  eval->do_gpoint_eval(false);
  eval->scan_for_embedded(false);  // take all tracks if false - take only embedded tracks if true
  eval->Verbosity(verbosity);
  se->registerSubsystem(eval);

  // MomentumEvaluator* eval = new MomentumEvaluator(outputfile.c_str(),0.2,0.4,Max_si_layer,2,Max_si_layer-4,10.,80.);
  // se->registerSubsystem( eval );

}
void Fast_Tracking_Eval(std::string outputfile, int verbosity = 0)
{
  gSystem->Load("libFastTrackingEval.so");

  Fun4AllServer *se = Fun4AllServer::instance();

  FastTrackingEval *fast_sim_eval = new FastTrackingEval("FastTrackingEval");
  fast_sim_eval->set_filename( outputfile.c_str() );
  se->registerSubsystem( fast_sim_eval );

  return;
}
