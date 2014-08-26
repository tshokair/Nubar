/*
 *  d2o_in.C
 *  
 *
 *  Created by Tim Shokair on 11/8/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *  Code will be used for a mock analysis, rewrite of the look back search.
 *  Looks at the light water and acrylic events 
 *  12/9 fixed problem with duplicates two fold events (problem is with the sorted root file,
 *  but the analysis code was fixed to make sure no coincidences occured with different event numbers and the same gtid)
 *  1/31 Fixed problem with smaller window yeilding more events, it was a problem with duplicates. 
 *  2/1 Version verified against forward search. Ready for analysis.  
 *	4/20 Updated to account for total vs. kinetic energy and input dr dt plot
 *  5/15 fixed typo that was overcounting two-fold
 *  6/10 write out more information to the files
 */

{ 
	
	gROOT->SetStyle("Plain");
	//snoman variables and temp snoman variables
	Float_t nhits;
	Float_t ene;
	Float_t eneu;
	Float_t enu;
	Float_t egen;
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
	Float_t xe;
	Float_t ye;
	Float_t ze;
	Float_t ri;
	Float_t xfit;
	Float_t yfit;
	Float_t zfit;
	Float_t rfit;
	Float_t rf;
	Float_t r0,r1,r2;
	Float_t dt1;
	Float_t dt2;
	Float_t dt3;
	Float_t rdamn1;
	Float_t rdamn2;
	Float_t evid;
	//counts;
	Int_t numCoin=0;
	Int_t numPred[500000]={0.0};
	Int_t coin_events[500000]={0.0};
	Int_t numTwoFold=0;
	Int_t numThreeFold=0;
	Int_t ct=0;
	Int_t doubleCounts=0;
	Int_t nd;
	Int_t numAcc;
	Int_t numCut;
	Int_t numPi=0;
	Int_t numOnePred=0;
	Int_t numTwoPred=0;
	Int_t twoFoldEn[15];
	Int_t threeFoldEn[15];
	//flags
	Int_t damnFlag1=0;
	Int_t damnFlag2=0;
	Int_t cleanFlag1=0;
	Int_t cleanFlag2=0;
	Int_t checkR=0;
	Int_t check2[500000]={0.0};
	Int_t twoFlag=0;
	Int_t threeFlag=0;
	Int_t duplicateFlag=0;
	Int_t duplicateFlag2=0;
	
	//times
	double totalTime;
	Int_t day_diff;
	double timediff;
	Double_t timeOfLast[500000][4];
	Double_t timeOfFirst[500000][4];
	//temps
	Int_t lastClean;
	Int_t tempID[500000];
	Int_t tempNhits;
	double tempTime[500000];
	Float_t tempEn[500000];
	double tempDt[500000];
	Int_t tID;
	Float_t tEn;
	double tTime;
	Float_t tBeta;
	Float_t tITR;
	Int_t l=0;
	Float_t pEnu;
	Float_t predX[100];
	Float_t predY[100];
	Float_t predZ[100];
	Int_t predID[100];
	Float_t predEn[100];
	double predDt[100];
	Float_t predNhit[100];
	Float_t predITR[100];
	Float_t predBeta[100];
	double evTime[100];
	Float_t deltaR[3];
	Float_t tempR;
	Float_t piEn;
	Float_t pe1En;
	Float_t pe2En;
	double piTime;
	double pe1Time;
	double pe2Time;
	Int_t piID;
	Int_t pe1ID;
	Int_t pe2ID;
	Int_t piNhit;
	Int_t pe1Nhit;
	Int_t pe2Nhit;
	Int_t lastCoinID1=0;
	Int_t lastCoinID2=0;
	
	//files
	Int_t numfiles;
	TString files[1000];
	
	//parameters
	Float_t nhitCut=20;
	Float_t nhitMax=150;
	Float_t fiducialVolume=600;
	Float_t timeWindow=0.150;
	Float_t distanceWindow=500;
	Float_t startWindow=100;
	Float_t betamin=-0.12;
	Float_t betamax=0.95;
	Float_t thetamin=0.75;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=4.5+.511;
	Float_t ethresh=3.9+.511;
	Float_t outerVolume=650;
	Float_t innerVolume=550;
	
	//readin the file names
	
	ifstream runList;
	//runList.open("d2o_data.runlist");
	//runList.open("d2o_partial.runlist");
	//runList.open("flatD2O.txt");
	runList.open("testNCD.txt");
	//runList.open("reactorInD2O.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	TH1F * h1= new TH1F("h1","Gen neutrino Energy",150,0,15);
	TH1F * h2= new TH1F("h2","4 Mev neutrino Energy",150,0,15);
	TH1F * h3= new TH1F("h3","5 Mev neutrino Energy",150,0,15);
	TH1F * h4= new TH1F("h4","6 Mev neutrino Energy",150,0,15);
	TH1F * h5= new TH1F("h5","7 Mev neutrino Energy",150,0,15);
	TH1F * h6= new TH1F("h6","8 Mev neutrino Energy",150,0,15);
	TH1F * h7= new TH1F("h7","9 Mev neutrino Energy",150,0,15);
	TH1F * h8= new TH1F("h8","10 Mev neutrino Energy",150,0,15);
	TH1F * h9= new TH1F("h9","11 Mev neutrino Energy",150,0,15);
	TH1F * h10= new TH1F("h10","12 Mev neutrino Energy",150,0,15);
	TH1F * h11= new TH1F("h11","13 Mev neutrino Energy",150,0,15);
	TH1F * h12= new TH1F("h12","14 Mev neutrino Energy",150,0,15);
	

	//graph points
	Double_t dt[100000];
	Double_t dr[100000];
	Int_t nGr=0;
	
	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k)//loop over ntuples/root files
    {	
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("t370");
		nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		//get the necessary information from ntuple
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Enekp",&ene);
		Td->SetBranchAddress("Enerau",&eneu);
		Td->SetBranchAddress("Enu",&enu);	
		Td->SetBranchAddress("Egen",&egen);
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
		Td->SetBranchAddress("Xe",&xe);
		Td->SetBranchAddress("Ye",&ye);
		Td->SetBranchAddress("Ze",&ze);
		Td->SetBranchAddress("Thetaijp",&thetaij);
		Td->SetBranchAddress("Beta1p",&beta1);
		Td->SetBranchAddress("Beta4p",&beta4);
		Td->SetBranchAddress("Rdamn1",&rdamn1);
		Td->SetBranchAddress("Rdamn2",&rdamn2);
		lastClean=0;
		
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			//cerr<<egen<<" "<<enu<<" "<<ene<<endl;
			h1->Fill(enu);
				
				if(Int_t(enu)==4){
					h2->Fill(ene);
				}
				else if(Int_t(enu)==5){
					h3->Fill(ene);
				}
				else if(Int_t(enu)==6){
					h4->Fill(ene);
				}
				else if(Int_t(enu)==7){
					h5->Fill(ene);
				}
				else if(Int_t(enu)==8){
					h6->Fill(ene);
				}
				else if(Int_t(enu)==9){
					h7->Fill(ene);
				}
				else if(Int_t(enu)==10){
					h8->Fill(ene);
				}
				else if(Int_t(enu)==11){
					h9->Fill(ene);
				}else if(Int_t(enu)==12){
					h10->Fill(ene);
				}else if(Int_t(enu)==13){
					h11->Fill(ene);
				}else if(Int_t(enu)==14){
					h12->Fill(ene);
				}

			//lastID=tID;
		}//end loop over events 
		delete Td;
		delete f1;
		cerr<<"next ntuple"<<endl;
		}//end loop over ntuples
		
		
		
		
		}//end code
