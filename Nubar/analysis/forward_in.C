/*
 *  search_in.C
 *  
 *
 *  Created by Tim Shokair on 11/8/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *  Code will be used for a mock analysis, rewrite of the look back search.
 *  Looks at the light water and acrylic events 
 *  12/9 fixed problem with duplicates two fold events (problem is with the sorted root file,
 *  but the analysis code was fixed to make sure no coincidences occured with different event numbers and the same gtid)
 *  1/31 Fixed problem with smaller window yeilding more events, it was a problem with duplicates. 
 */

{ 
	
	gROOT->SetStyle("Plain");
	//snoman variables and temp snoman variables
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
	Float_t xe;
	Float_t ye;
	Float_t ze;
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
	Int_t mult=0;
	Int_t lastMult=0;
	Int_t numPi=0;
	
	
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
	double lastTimediff;
	double timeOfTwo;
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
	Int_t tID2;
	Float_t tEn;
	double tTime;
	Int_t l=0;
	Float_t predX[100];
	Float_t predY[100];
	Float_t predZ[100];
	Int_t predID[100];
	Float_t predEn[100];
	double predDt[100];
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
	Float_t nhitMax=160;
	Float_t fiducialVolume=600;
	Float_t timeWindow=0.150;
	Float_t distanceWindow=500;
	Float_t startWindow=100;
	Float_t betamin=-0.12;
	Float_t betamax=0.95;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=2.7+.511;
	Float_t ethresh=2.7+.511;
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
	TH1F * h1= new TH1F("h1","Time Differences D2O",500,0,.5);
	TH1F * h2= new TH1F("h2","Primary Energy",100,0,20);
	TH1F * h3= new TH1F("h3","Predicessor Energy",100,0,20);
	TH1F * h4= new TH1F("h4","Coin IDs",30000,0,30000);
	TH1F * h5= new TH1F("h5","Coin Ev nums",30000,0,30000);
	TH1F * h6= new TH1F("h6","Delta R",1300,0,1300);
	TH1F * h7= new TH1F("h7","Reconstructed radius",1300,0,1300);
	TH1F * h8= new TH1F("h8","Generated radius",1300,0,1300);
	
	h2->SetLineColor(7);
	h3->SetLineColor(6);
	
	
	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k){//loop over ntuples/root files
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("t360");
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
		Td->SetBranchAddress("Xe",&xe);
		Td->SetBranchAddress("Ye",&ye);
		Td->SetBranchAddress("Ze",&ze);
		Td->SetBranchAddress("Thetaijp",&thetaij);
		Td->SetBranchAddress("Beta1p",&beta1);
		Td->SetBranchAddress("Beta4p",&beta4);
		Td->SetBranchAddress("Rdamn1",&rdamn1);
		Td->SetBranchAddress("Rdamn2",&rdamn2);
		lastClean=0;
		Td->GetEntry(0);
		ut1_i=ut1;
		ut2_i=ut2;
		ut3_i=ut3;
		jdy_i=jdy;
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			beta14=beta1+4*beta4;
			//establish a clean flag
			if (rf<fiducialVolume && ene>ethresh) {
				cleanFlag1=1;
			}
			else {
				cleanFlag1=0;
			}
			//if event is clean start a lookback
			if (cleanFlag1==1 ){
				numPred[numCoin]=0;
				day_diff=jdy-jdy_i;
				beta14=beta1+4*beta4;
				if(Int_t(day_diff)==-1){
					dt1=(86400-ut1_i)+ut1;
					if(Int_t(dt1)==0){
						dt2=ut2/1.E6-ut2_i/1.E6;
						dt3=ut3/1.E9-ut3_i/1.E9;
						timediff=dt1+dt2 +dt3;
					}
					else if(Int_t(dt1)==1){
						dt2=(1-ut2/1.E6)+ut2_i/1.E6;
						dt3=(1-ut3/1.E9)+ut3_i/1.E9;
						timediff=dt2+dt3;
					}
					else {timediff=10;}
				}
				//daydiff=0
				else if(Int_t(day_diff)==0){
					dt1=ut1_i-ut1;
					if(Int_t(dt1)==0 &&ut2>ut2_i){
						dt2=ut2/1.E6-ut2_i/1.E6;
						dt3=ut3/1.E9-ut3_i/1.E9;
						timediff=dt1+dt2+dt3;
					}
					else if(Int_t(dt1)==0 &&ut2<ut2_i){
						dt2=ut2_i/1.E6-ut2/1.E6;
						dt3=ut3_i/1.E9-ut3/1.E9;
						timediff=dt1+dt2+dt3;
					}
					else if(Int_t(dt1)==1){
						dt2=(1-ut2/1.E6)+ut2_i/1.E6;
						dt3=(1-ut3/1.E9)+ut3_i/1.E9;
						timediff=dt2+dt3;
					}
					else if(Int_t(dt1)==-1){
						dt2=(1-ut2_i/1.E6)+ut2/1.E6;
						dt3=(1-ut3_i/1.E9)+ut3/1.E9;
						timediff=dt2+dt3;
					}
					else timediff =10;
					
				}
				//daydiff =1
				
				else if(Int_t(day_diff)==1){
					dt1=(86400-ut1)+ut1_i;
					dt2=ut2;
					dt3=ut3;
					
					if(Int_t(dt1)==0){
						dt2=ut2_i/1.E6-ut2/1.E6;
						dt3=ut3_i/1.E9-ut3/1.E9;
						timediff=dt1+dt2 +dt3;
					}
					
					else if(Int_t(dt1)==1){
						dt2=(1-ut2/1.E6)+ut2_i/1.E6;
						dt3=(1-ut3/1.E9)+ut3_i/1.E9;
						timediff=dt2+dt3;
					}
					else timediff=10;
					
				}
				
				//cerr<<timediff<<endl;
				if (timediff<timeWindow) {
					//cerr<<evid<<" occured "<<timediff<< " s after last event"<<endl;
					mult++;
					//cerr<<mult<<endl;
				}
				else{
					//cerr<<timediff<<endl;
					mult=0;
					numPi++;
				}
				if (mult!=0) {
					//cerr<<mult+1<<"-fold event for "<<tID<<" at time "<<tTime<<endl;
					h4->Fill(evid);
					if (mult==1) {
						h4->Fill(tID);
					}
				}
				if (mult==1) {
					tID2=evid;
					lastTimediff=timediff;
					timeOfTwo=totalTime;
				}
				if (lastMult==1) {
					if(mult==2) {
						cerr<<"3-fold "<<endl;
						cerr<<tID<<" at time "<<tTime<<endl;
						cerr<<tID2<<" "<<lastTimediff<<" s later at time "<<timeOfTwo<<endl;
						cerr<<evid<<" "<<timediff<<" s later at time "<<totalTime<<endl;
						numThreeFold++;
						//cerr<<"last mult "<<lastMult<<endl;
						//h4->Fill(evid);
						//h4->Fill(tID2);
						//h4->Fill(tID);
					}
					else if (mult==0){
						cerr<<"2-fold"<<endl;
						cerr<<tID<<" at time "<<tTime<<endl;
						cerr<<tID2<<" "<<lastTimediff<<" s later at time "<<timeOfTwo<<endl;
						numTwoFold++;
						//h4->Fill(tID);
						//h4->Fill(evid);
					}
				}
				if (mult>2) {
					cerr<<mult<<"-fold event"<<endl;
				}
				lastMult=mult;
				if(mult==0){
					tTime=totalTime;
					tID=evid;
					ut1_i=ut1;
					ut2_i=ut2;
					ut3_i=ut3;
					jdy_i=jdy;
				}
				//nextlastID=lastID;
				//lastID=tID;
			}//end if statement on clean
		}//end loop over events 
		delete Td;
		delete f1;
		cerr<<"next ntuple"<<endl;
		//end loop over ntuples
		}
		cerr<<numTwoFold<<" two fold events"<<endl;
		cerr<<numThreeFold<<" three fold events"<<endl;
		//h3->Draw();
		//h2->Draw("same");
		h4->Draw();
		cerr<<numPi<<" primaries"<<endl;
		
	}//end code