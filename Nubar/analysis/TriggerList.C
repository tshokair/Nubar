/*
 *  TriggerList.h
 *  
 *
 *  Created by Tim Shokair on 1/9/12.
 *  Copyright 2012 University of Pennsylvania. All rights reserved.
 *
 */

{ 
	
	gROOT->SetStyle("Plain");
	Float_t nhits;
	Float_t ene;
	Float_t eneu;
	Float_t jdy;
	Float_t ut1;
	Float_t ut2;
	Float_t ut3;
	Float_t jdy_i;
	Float_t ut1_i;
	Float_t ut2_i;
	Float_t ut3_i;	
	Float_t itr;
	Float_t thetaij;
	Float_t beta1;
	Float_t beta4;
	Float_t beta14;
	Float_t xi;
	Float_t yi;
	Float_t zi;
	Float_t ri;
	Float_t xfit;
	Float_t yfit;
	Float_t zfit;
	Float_t rfit;
	Float_t rf;
	Float_t dt1;
	Float_t dt2;
	Float_t dt3;
	Float_t rdamn1;
	Float_t rdamn2;
	Float_t evid;
	
	//files
	Int_t numfiles;
	TString files[1000];
	
	//counters
	Int_t numTrig=0;
	
	//parameters
	Float_t nhitCut=20;
	Float_t nhitMax=160;
	Float_t fiducialVolume=600;
	Float_t timeWindow=0.050;
	Float_t distanceWindow=200;
	Float_t startWindow=100;
	Float_t betamin=-0.12;
	Float_t betamax=0.95;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=3.5;
	Float_t ethresh=2.7;
	Float_t outerVolume=650;
	Float_t innerVolume=550;
	Float_t sourceX=0.00;
	Float_t sourceY=-21.59;
	Float_t sourceZ=0.25;
	Float_t sourceBox=30.0;
	
	//readin the file names
	ifstream runList;
	runList.open("testcf.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	ofstream trigList;
	trigList.open("trigList.txt");

	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k)//loop over ntuples/root files
    {	
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("h360");
		nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		//get the necessary information from ntuple
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Enekp",&ene);
		Td->SetBranchAddress("Enerau",&eneu);
		Td->SetBranchAddress("Itrp",&itr);
		Td->SetBranchAddress("Ev_tid",&evid);	
		Td->SetBranchAddress("Time_jdy",&jdy); //julian day
		Td->SetBranchAddress("Time_s",&ut1); //time in seconds
		Td->SetBranchAddress("Time_us",&ut2); //time in nanoseconds
		Td->SetBranchAddress("Time_ns",&ut3); //time in nanoseconds
		Td->SetBranchAddress("Xfp",&xfit);
		Td->SetBranchAddress("Yfp",&yfit);
		Td->SetBranchAddress("Zfp",&zfit);
		Td->SetBranchAddress("Rfp",&rf);
		Td->SetBranchAddress("Thetaijp",&thetaij);
		Td->SetBranchAddress("Beta1p",&beta1);
		Td->SetBranchAddress("Beta4p",&beta4);
		Td->SetBranchAddress("Rdamn1",&rdamn1);
		Td->SetBranchAddress("Rdamn2",&rdamn2);
		for(Int_t j=0; j<nd; ++j){//loop over events
			Td->GetEntry(j);
			beta14=beta1+4*beta4;
			if (TMath::Abs(xfit-sourceX)<sourceBox &&TMath::Abs(yfit-sourceY)<sourceBox&&TMath::Abs(zfit-sourceZ)<sourceBox ) {
				//cerr<<"Triggered event :"<<j<<endl;
				numTrig++;
				trigList<<j<<endl;
			}
		
		}//end loop over events
	}//end loop over files
		
		cerr<<numTrig<<" source events of "<<nd<<" total events in file"<<endl;
		trigList.close();
}//end code