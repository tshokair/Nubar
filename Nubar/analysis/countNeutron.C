/*
 *  lookBack_salt_out.h
 *  
 *
 *  Created by Tim Shokair on 11/11/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *  Code will count the number of single neutron background events from MC
 *  
 */

{ 

	gROOT->SetStyle("Plain");
	//snoman variables and temp snoman variables
	Float_t nhits;
	Float_t ene;
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
	//counts;
	Int_t numNeutrons=0;
	//flags
	Int_t damnFlag1=0;
	Int_t damnFlag2=0;
	Int_t cleanFlag1=0;
	Int_t cleanFlag2=0;
	Int_t checkR=0;
	Int_t check2[50000]={0.0};
	Int_t twoFlag=0;
	Int_t threeFlag=0;
	//times
	double liveTime=0;
	double runTime=0;
	double timeOfLast=0;
	double startTime=0;
	double endTime=0;
	Int_t startDay=0;
	Int_t endDay;
	Int_t totalDays;
	//files
	Int_t numfiles;
	TString files[2000];
	
	//parameters
	Float_t nhitCut=20;
	Float_t nhitMax=160;
	Float_t fiducialVolume=550;
	Float_t timeWindow=0.050;
	Float_t distanceWindow=200;
	Float_t startWindow=100;
	Float_t thetamin=0.63;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=3.5;
	Float_t ethresh=2.7;
	Float_t emax=100;
	Float_t outerVolume=650;
	Float_t innerVolume=550;

	//readin the file names
	ifstream runList;
	runList.open("test.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	TH1F * h1= new TH1F("h1","Time Differences Salt",300,0,.1);
	TH1F * h2= new TH1F("h2","Primary Energy",100,0,20);
	TH1F * h3= new TH1F("h3","Predicessor Energy",100,0,20);
	TH1F * h4= new TH1F("h4","Coin ID's",30000,0,15000);

	h2->SetLineColor(7);
	h3->SetLineColor(6);


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
		Td->SetBranchAddress("Rdamn1",&rdamn1);
		Td->SetBranchAddress("Rdamn2",&rdamn2);
		Td->GetEntry(0);
		lastClean = evid;
		startTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
		timeOfLast=startTime;
		startDay=jdy;
		runTime=0;
		cerr<<" First event at time "<<startTime<<" on day "<<jdy<<endl;
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			//establish a clean flag
			cerr<<"ID: "<<evid<<endl;
			if (rf<fiducialVolume &&ene>emin &&nhits<nhitMax &&ene<emax) {
				cleanFlag1=1;
			}
			else {
				cleanFlag1=0;
			}
			//if event is clean start a lookback
			if (cleanFlag1==1){
				numNeutrons++;
			}
			endDay=jdy;
			endTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;

		}
		//get the run and livetime
		totalDays=endDay-startDay;
		if(totalDays==0){
			runTime=endTime-startTime;
		}
		else if( totalDays>0){
			runTime=(864000.0-startTime)+(totalDays-1)*86400.0+endTime;
		}
		else {runTime=0;}
		liveTime+=runTime;
		delete Td;
		delete f1;
	}

		cerr<<" Last event at time "<<endTime<<" on day "<<endDay<<endl;
		cerr<<numNeutrons<<" neutrons in time "<<liveTime<<endl;
		cerr<<"Neutron rate is "<<numNeutrons/liveTime<<endl;
}//end code