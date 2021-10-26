
// ANCHOR debug output verbosity
//Int_t verbosityERH = 0;

// ANCHOR create histograms globally
TH2F*  h_E_Minv[_active_calo][_active_algo];
TH2F*  h_Erecon_Minv[_active_calo][_active_algo];
TH2F*  h_E_Erecon2photon[_active_calo][_active_algo];
TH2F*  h_E_Ereconsinglecluster[_active_calo][_active_algo];
TH2F*  h_E_M20_photon1[_active_calo][_active_algo];
TH2F*  h_E_M20_photon2[_active_calo][_active_algo];
TH2F*  h_E_M20_singlecluster[_active_calo][_active_algo];
TH2F*  h_E_M02_photon1[_active_calo][_active_algo];
TH2F*  h_E_M02_photon2[_active_calo][_active_algo];
TH2F*  h_E_M02_singlecluster[_active_calo][_active_algo];

TH3F*  h_E_eta_Minv[_active_calo][_active_algo];
TH3F*  h_Erecon_eta_Minv[_active_calo][_active_algo];

TH2F*  h_E_Clusters[_active_calo][_active_algo];

TH2F*  h_Pt_Minv[_active_calo][_active_algo];
TH2F*  h_Pt_Erecon2photon[_active_calo][_active_algo];
TH2F*  h_Pt_Ereconsinglecluster[_active_calo][_active_algo];
TH2F*  h_Pt_M20_photon1[_active_calo][_active_algo];
TH2F*  h_Pt_M20_photon2[_active_calo][_active_algo];
TH2F*  h_Pt_M20_singlecluster[_active_calo][_active_algo];
TH2F*  h_Pt_M02_photon1[_active_calo][_active_algo];
TH2F*  h_Pt_M02_photon2[_active_calo][_active_algo];
TH2F*  h_Pt_M02_singlecluster[_active_calo][_active_algo];
TH3F*  h_Pt_eta_Minv[_active_calo][_active_algo];
TH2F*  h_Pt_Clusters[_active_calo][_active_algo];

TH3F*  h_E_eta_Clusters[_active_calo][_active_algo];
TH3F*  h_Pt_eta_Clusters[_active_calo][_active_algo];

void pi0studies(){

  CheckBarrelCaloVersion();
  
  for(int icalo=0;icalo<_active_calo;icalo++){
  
    //=====Loop only over enabled calorimeters ===========//
    if (!caloEnabled[icalo]) continue;

    for(int ialgo=0;ialgo<_active_algo;ialgo++){

      //if(verbosityERH) cout << "Calo-type: " << icalo << "\t Algorithm: " << ialgo << endl;

      if(!loadClusterizerInput(ialgo,icalo )) continue;
      
      std::sort(_clusters_calo[ialgo][icalo].begin(), _clusters_calo[ialgo][icalo].end(), &acompareCl);

      
      int nbinsx = 400;
      if(icalo==kFHCAL) nbinsx = 200;
      if(icalo==kEHCAL || icalo==kHCALOUT || icalo==kHCALOUT) nbinsx = 100;
      if(icalo==kLFHCAL) nbinsx = 200;
      if(icalo==kEEMC) nbinsx = 200;
      if(icalo==kBECAL) nbinsx = 200;
      if(icalo==kFEMC) nbinsx = 200;
      

      if(!h_E_Minv[icalo][ialgo]) h_E_Minv[icalo][ialgo]  = new TH2F(Form("h_E_Minv_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, nbinsx, 0, 1);
      if(!h_Erecon_Minv[icalo][ialgo]) h_Erecon_Minv[icalo][ialgo]  = new TH2F(Form("h_Erecon_Minv_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, nbinsx, 0, 1);
      if(!h_E_Erecon2photon[icalo][ialgo]) h_E_Erecon2photon[icalo][ialgo]  = new TH2F(Form("h_E_Erecon2photon_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, nbinsx, 0, 5);
      if(!h_E_Ereconsinglecluster[icalo][ialgo]) h_E_Ereconsinglecluster[icalo][ialgo]  = new TH2F(Form("h_E_Ereconsinglecluster_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, nbinsx, 0, 5);
      if(!h_E_M20_photon1[icalo][ialgo]) h_E_M20_photon1[icalo][ialgo]  = new TH2F(Form("h_E_M20_photon1_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);
      if(!h_E_M20_photon2[icalo][ialgo]) h_E_M20_photon2[icalo][ialgo]  = new TH2F(Form("h_E_M20_photon2_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);
      if(!h_E_M20_singlecluster[icalo][ialgo]) h_E_M20_singlecluster[icalo][ialgo]  = new TH2F(Form("h_E_M20_singlecluster_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);
      if(!h_E_M02_photon1[icalo][ialgo]) h_E_M02_photon1[icalo][ialgo]  = new TH2F(Form("h_E_M02_photon1_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);

      if(!h_E_M02_photon2[icalo][ialgo]) h_E_M02_photon2[icalo][ialgo]  = new TH2F(Form("h_E_M02_photon2_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);
      if(!h_E_M02_singlecluster[icalo][ialgo]) h_E_M02_singlecluster[icalo][ialgo]  = new TH2F(Form("h_E_M02_singlecluster_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);

      if(!h_E_eta_Minv[icalo][ialgo]) h_E_eta_Minv[icalo][ialgo]  = new TH3F(Form("h_E_Eta_Minv_singlecluster_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5, 500, 0, 1); 
      if(!h_Erecon_eta_Minv[icalo][ialgo]) h_Erecon_eta_Minv[icalo][ialgo]  = new TH3F(Form("h_Erecon_Etarecon_Minv_singlecluster_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5, 500, 0, 1); 
      if(!h_E_Clusters[icalo][ialgo]) h_E_Clusters[icalo][ialgo]  = new TH2F(Form("h_E_Clusters_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 100, 0, 10);


      if(!h_Pt_Minv[icalo][ialgo]) h_Pt_Minv[icalo][ialgo]  = new TH2F(Form("h_pT_Minv_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, nbinsx, 0, 1);
      if(!h_Pt_Erecon2photon[icalo][ialgo]) h_Pt_Erecon2photon[icalo][ialgo]  = new TH2F(Form("h_pT_Erecon2photon_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, nbinsx, 0, 5);
      if(!h_Pt_Ereconsinglecluster[icalo][ialgo]) h_Pt_Ereconsinglecluster[icalo][ialgo]  = new TH2F(Form("h_pT_Ereconsinglecluster_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, nbinsx, 0, 5);
      if(!h_Pt_M20_photon1[icalo][ialgo]) h_Pt_M20_photon1[icalo][ialgo]  = new TH2F(Form("h_pT_M20_photon1_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);
      if(!h_Pt_M20_photon2[icalo][ialgo]) h_Pt_M20_photon2[icalo][ialgo]  = new TH2F(Form("h_pT_M20_photon2_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);
      if(!h_Pt_M20_singlecluster[icalo][ialgo]) h_Pt_M20_singlecluster[icalo][ialgo]  = new TH2F(Form("h_pT_M20_singlecluster_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);
      if(!h_Pt_M02_photon1[icalo][ialgo]) h_Pt_M02_photon1[icalo][ialgo]  = new TH2F(Form("h_pT_M02_photon1_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);
      if(!h_Pt_M02_photon2[icalo][ialgo]) h_Pt_M02_photon2[icalo][ialgo]  = new TH2F(Form("h_pT_M02_photon2_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);
      if(!h_Pt_M02_singlecluster[icalo][ialgo]) h_Pt_M02_singlecluster[icalo][ialgo]  = new TH2F(Form("h_pT_M02_singlecluster_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5);
      if(!h_Pt_eta_Minv[icalo][ialgo]) h_Pt_eta_Minv[icalo][ialgo]  = new TH3F(Form("h_pT_Eta_Minv_singlecluster_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5, 500, 0, 1);
      if(!h_Pt_Clusters[icalo][ialgo]) h_Pt_Clusters[icalo][ialgo]  = new TH2F(Form("h_pT_Clusters_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 100, 0, 10);


      if(!h_E_eta_Clusters[icalo][ialgo]) h_E_eta_Clusters[icalo][ialgo]  = new TH3F(Form("h_E_eta_Clusters_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5, 500, 0, 10); 
      if(!h_Pt_eta_Clusters[icalo][ialgo]) h_Pt_eta_Clusters[icalo][ialgo]  = new TH3F(Form("hPt_eta_Clusters__%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 250, 500, -5, 5, 500, 0, 10); 
      TLorentzVector pi0Gen(_mcpart_px[0],_mcpart_py[0],_mcpart_pz[0],_mcpart_E[0]);
	
      if(_clusters_calo[ialgo][icalo].size() == 2){

	double Erec = (_clusters_calo[ialgo][icalo].at(0)).cluster_E + (_clusters_calo[ialgo][icalo].at(1)).cluster_E ;
	double E1 = (_clusters_calo[ialgo][icalo].at(0)).cluster_E;
	double E2 = (_clusters_calo[ialgo][icalo].at(1)).cluster_E;
	double theta1 = 2*atan(exp(-(_clusters_calo[ialgo][icalo].at(0)).cluster_Eta));
	double theta2 = 2*atan(exp(-(_clusters_calo[ialgo][icalo].at(1)).cluster_Eta));
	double phi1   = (_clusters_calo[ialgo][icalo].at(0)).cluster_Phi;
	double phi2   = (_clusters_calo[ialgo][icalo].at(1)).cluster_Phi;

	double pX1    = E1*sin(theta1)*cos(phi1);
	double pY1    = E1*sin(theta1)*sin(phi1);
	double pZ1    = E1*cos(theta1);
	TLorentzVector photon1(pX1, pY1, pZ1, E1);
	
	double pX2    = E2*sin(theta2)*cos(phi2);
        double pY2    = E2*sin(theta2)*sin(phi2);
        double pZ2    = E2*cos(theta2);
        TLorentzVector photon2(pX2, pY2, pZ2, E2);

	TLorentzVector pi0 = photon1 + photon2;
	double Minv = pi0.M();
	//if(ialgo==1) cout << _clusters_calo[ialgo][icalo].size() << " " << Minv << endl;
	
	//cout << "====== M ========" << Minv << " " << photon1.M() << " " << photon1.Eta() << " " <<  (_clusters_calo[ialgo][icalo].at(0)).cluster_Eta  << " " << photon1.Phi() << " " << phi1 << endl;	
	//cout << "====== M ========" << Minv << " " << photon2.M() << " " << photon2.Eta() << " " <<  (_clusters_calo[ialgo][icalo].at(1)).cluster_Eta  << " " << photon2.Phi() << " " << phi2 << endl;	
	
	h_Erecon_Minv[icalo][ialgo]->Fill(Erec,Minv);
	h_E_Minv[icalo][ialgo]->Fill(_mcpart_E[0],Minv);
	h_E_Erecon2photon[icalo][ialgo]->Fill(_mcpart_E[0],Erec);
	h_E_M20_photon1[icalo][ialgo]->Fill(_mcpart_E[0],(_clusters_calo[ialgo][icalo].at(0)).cluster_M20);
	h_E_M20_photon2[icalo][ialgo]->Fill(_mcpart_E[0],(_clusters_calo[ialgo][icalo].at(1)).cluster_M20);
	h_E_M02_photon1[icalo][ialgo]->Fill(_mcpart_E[0],(_clusters_calo[ialgo][icalo].at(0)).cluster_M02);
	h_E_M02_photon2[icalo][ialgo]->Fill(_mcpart_E[0],(_clusters_calo[ialgo][icalo].at(1)).cluster_M02);
	h_E_eta_Minv[icalo][ialgo]->Fill(_mcpart_E[0],_mcpart_Eta[0],Minv);
	h_Erecon_eta_Minv[icalo][ialgo]->Fill(pi0.E(),pi0.Eta(),Minv);

	h_Pt_Minv[icalo][ialgo]->Fill(pi0Gen.Pt(),Minv);
 	h_Pt_Erecon2photon[icalo][ialgo]->Fill(pi0Gen.Pt(),Erec);
	h_Pt_M20_photon1[icalo][ialgo]->Fill(pi0Gen.Pt(),(_clusters_calo[ialgo][icalo].at(0)).cluster_M20);;
	h_Pt_M20_photon2[icalo][ialgo]->Fill(pi0Gen.Pt(),(_clusters_calo[ialgo][icalo].at(1)).cluster_M20);
	h_Pt_M02_photon1[icalo][ialgo]->Fill(pi0Gen.Pt(),(_clusters_calo[ialgo][icalo].at(0)).cluster_M02);
	h_Pt_M02_photon2[icalo][ialgo]->Fill(pi0Gen.Pt(),(_clusters_calo[ialgo][icalo].at(1)).cluster_M02);
	h_Pt_eta_Minv[icalo][ialgo]->Fill(pi0Gen.Pt(),_mcpart_Eta[0],Minv);
	
      }

      if( _clusters_calo[ialgo][icalo].size() == 1){
	h_E_Ereconsinglecluster[icalo][ialgo]->Fill(_mcpart_E[0],(_clusters_calo[ialgo][icalo].at(0)).cluster_E);
	h_E_M20_singlecluster[icalo][ialgo]->Fill(_mcpart_E[0],(_clusters_calo[ialgo][icalo].at(0)).cluster_M20);
	h_E_M02_singlecluster[icalo][ialgo]->Fill(_mcpart_E[0],(_clusters_calo[ialgo][icalo].at(0)).cluster_M02);
       
	h_Pt_Ereconsinglecluster[icalo][ialgo]->Fill(pi0Gen.Pt(),(_clusters_calo[ialgo][icalo].at(0)).cluster_E);
	h_Pt_M20_singlecluster[icalo][ialgo]->Fill(pi0Gen.Pt(),(_clusters_calo[ialgo][icalo].at(0)).cluster_M20);
	h_Pt_M02_singlecluster[icalo][ialgo]->Fill(pi0Gen.Pt(),(_clusters_calo[ialgo][icalo].at(0)).cluster_M02);
      }
      h_E_Clusters[icalo][ialgo]->Fill(_mcpart_E[0],_clusters_calo[ialgo][icalo].size());
      h_Pt_Clusters[icalo][ialgo]->Fill(pi0Gen.Pt(),_clusters_calo[ialgo][icalo].size());	
      h_E_eta_Clusters[icalo][ialgo]->Fill(pi0Gen.E(),_mcpart_Eta[0],_clusters_calo[ialgo][icalo].size());
      h_Pt_eta_Clusters[icalo][ialgo]->Fill(pi0Gen.Pt(),_mcpart_Eta[0],_clusters_calo[ialgo][icalo].size());
   
    }
  }
}

// ANCHOR save function after event loop

void pi0studiesSave(){
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_pi0.root",outputDir.Data()),"RECREATE");

  // write histograms
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      //====== Eta =====================================================/
      if(h_E_Minv[icalo][ialgo])  h_E_Minv[icalo][ialgo]->Write();
      if(h_Erecon_Minv[icalo][ialgo]) h_Erecon_Minv[icalo][ialgo]->Write();
      if(h_E_Erecon2photon[icalo][ialgo]) h_E_Erecon2photon[icalo][ialgo]->Write();
      if(h_E_Ereconsinglecluster[icalo][ialgo]) h_E_Ereconsinglecluster[icalo][ialgo]->Write();
      if(h_E_M20_photon1[icalo][ialgo]) h_E_M20_photon1[icalo][ialgo]->Write();
      if(h_E_M20_photon2[icalo][ialgo]) h_E_M20_photon2[icalo][ialgo]->Write();
      if(h_E_M20_singlecluster[icalo][ialgo]) h_E_M20_singlecluster[icalo][ialgo]->Write();
      if(h_E_M02_photon1[icalo][ialgo]) h_E_M02_photon1[icalo][ialgo]->Write();
      if(h_E_M02_photon2[icalo][ialgo]) h_E_M02_photon2[icalo][ialgo]->Write();
      if(h_E_M02_singlecluster[icalo][ialgo])h_E_M02_singlecluster[icalo][ialgo]->Write();
      
      if(h_E_eta_Minv[icalo][ialgo]) h_E_eta_Minv[icalo][ialgo]->Write();
      if(h_Erecon_eta_Minv[icalo][ialgo]) h_Erecon_eta_Minv[icalo][ialgo]->Write();
      if(h_E_Clusters[icalo][ialgo]) h_E_Clusters[icalo][ialgo]->Write(); 
      if(h_Pt_Minv[icalo][ialgo]) h_Pt_Minv[icalo][ialgo]->Write();
      if(h_Pt_Erecon2photon[icalo][ialgo]) h_Pt_Erecon2photon[icalo][ialgo]->Write();
      if(h_Pt_Ereconsinglecluster[icalo][ialgo]) h_Pt_Ereconsinglecluster[icalo][ialgo]->Write();
      if(h_Pt_M20_photon1[icalo][ialgo]) h_Pt_M20_photon1[icalo][ialgo]->Write();
      if(h_Pt_M20_photon2[icalo][ialgo]) h_Pt_M20_photon2[icalo][ialgo]->Write();
      if(h_Pt_M20_singlecluster[icalo][ialgo]) h_Pt_M20_singlecluster[icalo][ialgo]->Write();
      if(h_Pt_M02_photon1[icalo][ialgo]) h_Pt_M02_photon1[icalo][ialgo]->Write();
      if(h_Pt_M02_photon2[icalo][ialgo]) h_Pt_M02_photon2[icalo][ialgo]->Write();
      if(h_Pt_M02_singlecluster[icalo][ialgo]) h_Pt_M02_singlecluster[icalo][ialgo]->Write();
      if(h_Pt_eta_Minv[icalo][ialgo]) h_Pt_eta_Minv[icalo][ialgo]->Write();
      if(h_Pt_Clusters[icalo][ialgo]) h_Pt_Clusters[icalo][ialgo]->Write();

      if(h_E_eta_Clusters[icalo][ialgo]) h_E_eta_Clusters[icalo][ialgo]->Write();
      if(h_Pt_eta_Clusters[icalo][ialgo]) h_Pt_eta_Clusters[icalo][ialgo]->Write();


    }	
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
