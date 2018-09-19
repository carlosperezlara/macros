using namespace std;

void
EGEM_Init()
{

}

void
FGEM_Init()
{

}

void
EGEMSetup(PHG4Reco* g4Reco)
{
  /* Careful with dimensions! If GEM station volumes overlap, e.g. with TPC volume, they will be
   * drawn in event display but will NOT register any hits.
   *
   * Geometric constraints:
   * TPC length = 211 cm --> from z = -105.5 to z = +105.5
   */
  float thickness=3.;
  PHG4SectorSubsystem *gem;

  make_GEM_station("EGEM_0", g4Reco,  -20.5 + thickness, -0.94, -1.95);
  make_GEM_station("EGEM_1", g4Reco,  -69.5 + thickness, -2.07, -3.21);
  gem = make_GEM_station("EGEM_2O", g4Reco, -137.0 + thickness, -1.4, -2.5, 8, 0, true);
  gem->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  gem = make_GEM_station("EGEM_2I", g4Reco, -137.0 + thickness, -2.5, -3.9, 8, 0, true);
  gem->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  gem = make_GEM_station("EGEM_3O", g4Reco, -160.0 + thickness, -1.5, -2.50, 8, 0, true);
  gem->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  gem = make_GEM_station("EGEM_3I", g4Reco, -160.0 + thickness, -2.5, -4.00, 8, 0, true);
  gem->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
}

void
FGEMSetup(PHG4Reco* g4Reco)//
//const double min_eta = 1.245 //
//)
{
  int N_Sector = 8;
  double min_eta = 1.245;
  double tilt = .1;
  PHG4SectorSubsystem *gm;
  
  string name;
  double etamax;
  double etamin;
  double zpos;
  double theta;
  double alpha0, alpha1, alpha2;
  double ta0, ta1, ta2, ct, st, tt;
  PHG4SectorSubsystem *gem;
  st  = TMath::Sin(tilt);
  ct  = TMath::Cos(tilt);
  tt  = TMath::Tan(tilt);

  make_GEM_station("FGEM_0", g4Reco, 17.5, 0.94, 2.73, N_Sector);
  make_GEM_station("FGEM_1", g4Reco, 66.5, 2.07, 4.00, N_Sector);
  
  ////////////////////////////// FGEM 2 ///////////////////////////////////
  alpha0 = 0.5*( +PHG4Sector::Sector_Geometry::eta_to_polar_angle(1.245)
		 +PHG4Sector::Sector_Geometry::eta_to_polar_angle(4.000) ); // 1.245 - 4.0
  alpha1 = 0.5*( +PHG4Sector::Sector_Geometry::eta_to_polar_angle(1.245)
		 +PHG4Sector::Sector_Geometry::eta_to_polar_angle(2.500) ); // 1.245 - 2.5 
  alpha2 = 0.5*( +PHG4Sector::Sector_Geometry::eta_to_polar_angle(2.500)
		 +PHG4Sector::Sector_Geometry::eta_to_polar_angle(4.000) ); // 2.5 - 4.0
  ta0 = TMath::Tan(alpha0);
  ta1 = TMath::Tan(alpha1);
  ta2 = TMath::Tan(alpha2);
  zpos = 134*(1 - tt*(ta1-ta0)/(1+ta1*tt) );
  gem = make_GEM_station("FGEM_2O", g4Reco, zpos, 1.245, 2.50, N_Sector, tilt, true);
  gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm(), 0);
  gem->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  zpos = 134*(1 + tt*(ta0-ta2)/(1+ta2*tt));
  gem = make_GEM_station("FGEM_2I", g4Reco, zpos, 2.500, 4.00, N_Sector, tilt, true);  
  gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm(), 0);
  gem->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());

  ////////////////////////////// FGEM 3 ///////////////////////////////////
  alpha0 = 0.5*( +PHG4Sector::Sector_Geometry::eta_to_polar_angle(2.000)
		 +PHG4Sector::Sector_Geometry::eta_to_polar_angle(4.000) ); // 2.0 - 4.0
  alpha1 = 0.5*( +PHG4Sector::Sector_Geometry::eta_to_polar_angle(2.000)
		 +PHG4Sector::Sector_Geometry::eta_to_polar_angle(2.500) ); // 2.0 - 2.5
  alpha2 = 0.5*( +PHG4Sector::Sector_Geometry::eta_to_polar_angle(2.500)
		 +PHG4Sector::Sector_Geometry::eta_to_polar_angle(4.000) ); // 2.5 - 4.0
  ta0 = TMath::Tan(alpha0);
  ta1 = TMath::Tan(alpha1);
  ta2 = TMath::Tan(alpha2);
  zpos = 157*(1 + tt*(ta0-ta2)/(1+ta2*tt));
  gem = make_GEM_station("FGEM_3I", g4Reco, zpos, 2.500, 4.00, N_Sector, tilt, true);
  //gem = make_GEM_station("FGEM_3I", g4Reco, zpos, 1.250, 4.0);
  gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm(), 0);
  gem->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());

  //==
  if(1) {
    zpos = 157*(1 - tt*(ta1-ta0)/(1+ta1*tt) );
    gem = make_GEM_station("FGEM_3O", g4Reco, zpos, 2.000, 2.50, N_Sector, tilt, true, "FGEM_3OL");
    gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm(), 0);
    //==
    double alphak = PHG4Sector::Sector_Geometry::eta_to_polar_angle(2.000);
    double alpham = PHG4Sector::Sector_Geometry::eta_to_polar_angle(1.245);
    zpos = 157*(1 - (st + ct*tan(alphak-tilt))*st);
    gem = make_GEM_station("FGEM_3O", g4Reco, zpos, 1.245, 2.00, N_Sector, (alpham+alphak)/2, true, "FGEM_3OU");
    gem->get_geometry().set_normal_start(zpos*PHG4Sector::Sector_Geometry::Unit_cm(), alphak);
    gem->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  }
  ////////////////////////////// FGEM 4 ///////////////////////////////////
  zpos = 271*(1 + tt*(ta0-ta2)/(1+ta2*tt));
  gem = make_GEM_station("FGEM_4I", g4Reco, zpos, 2.500, 4.00, N_Sector, tilt, true);
  //gem = make_GEM_station("FGEM_4I", g4Reco, zpos, 1.2500, 4.00);
  gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm(), 0);
  gem->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  //==
  if(1) {
    zpos = 271*(1 - tt*(ta1-ta0)/(1+ta1*tt) );
    gem = make_GEM_station("FGEM_4O", g4Reco, zpos, 2.000, 2.50, N_Sector, tilt, true, "FGEM_4OL");
    gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm(), 0);
    //==
    zpos = 271*(1 - (st + ct*tan(alphak-tilt))*st);
    gem = make_GEM_station("FGEM_4O", g4Reco, zpos, 1.245, 2.00, N_Sector, (alpham+alphak)/2, true, "FGEM_4OU");
    gem->get_geometry().set_normal_start(zpos*PHG4Sector::Sector_Geometry::Unit_cm(), alphak);
    gem->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  }
  ///////////////////////////////////////////////////////////////////////////
}

//! Add drift layers to mini TPC
void
AddLayers_MiniTPCDrift(PHG4SectorSubsystem *gem)
{
  assert(gem);

  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;

  //  const int N_Layers = 70; // used for mini-drift TPC timing digitalization
  const int N_Layers = 1; // simplified setup
  const double thickness = 2 * cm;

  gem->get_geometry().AddLayer("EntranceWindow", "G4_MYLAR", 25 * um, false,
                               100);
  gem->get_geometry().AddLayer("Cathode", "G4_GRAPHITE", 10 * um, false, 100);

  for (int d = 1; d <= N_Layers; d++)
    {
      stringstream s;
      s << "DriftLayer_";
      s << d;

      gem->get_geometry().AddLayer(s.str(), "G4_METHANE", thickness / N_Layers,
                                   true);

    }
}

PHG4SectorSubsystem* make_GEM_station(string name, PHG4Reco* g4Reco, double zpos, double etamin,
				      double etamax, int N_Sector=8, double extratilt=0,
				      bool skip=false, string subname="none") {
  if(subname.compare("none")==0) subname = name;
  double polar_angle = 0;
  if (zpos < 0) {
    zpos = -zpos;
    polar_angle = TMath::Pi();
  }
  if(etamax < etamin) {
    double t = etamax;
    etamax = etamin;
    etamin = t;
  }
  double tanmin = TMath::Tan(2*TMath::ATan(TMath::Exp(-etamax)));
  double tanmax = TMath::Tan(2*TMath::ATan(TMath::Exp(-etamin)));
  cout << " Making " << name.c_str() << " " << subname.c_str() << endl;
  cout << "   z=" << zpos << "  Rmin=" << zpos*tanmin << "  Rmax=" << zpos*tanmax << endl << endl;
  polar_angle += extratilt;
  PHG4SectorSubsystem *gem;
  gem = new PHG4SectorSubsystem(subname.c_str());
  gem->SuperDetector(name);
  gem->get_geometry().set_normal_polar_angle(polar_angle);
  gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
  gem->get_geometry().set_min_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamax));
  gem->get_geometry().set_max_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamin));
  if(!skip) {
    gem->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
    gem->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
  }
  gem->get_geometry().set_N_Sector(N_Sector);
  gem->get_geometry().set_material("G4_METHANE");
  gem->OverlapCheck(overlapcheck);
  AddLayers_MiniTPCDrift(gem);
  gem->get_geometry().AddLayers_HBD_GEM();
  g4Reco->registerSubsystem(gem);
  return gem;
}
