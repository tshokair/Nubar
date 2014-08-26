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
	Float_t rut1;
	Float_t rut2;
	Float_t rut3;
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
	Float_t egen;
	Float_t rgen;
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
	Int_t numPi=0;
	Int_t numOnePred=0;
	Int_t numTwoPred=0;
	
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
	Float_t tR;
	double tTime;
	Int_t l=0;
	Float_t predX[100];
	Float_t predY[100];
	Float_t predZ[100];
	Float_t predR[100];
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
	double rtrigTime=0;

	//files
	Int_t numfiles;
	TString files[1000];
	
	//parameters
	Float_t nhitCut=20;
	Float_t nhitMax=160;
	Float_t fiducialVolume=590;
	Float_t timeWindow=.1;
	Float_t distanceWindow=1000;
	Float_t startWindow=100;
	Float_t betamin=-0.12;
	Float_t betamax=0.95;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=0.0;
	Float_t ethresh=0.0;
	Float_t outerVolume=650;
	Float_t innerVolume=550;
	Double_t retrigTime=0.00005;

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
	TH2F * h= new TH2F("h"," Primary Energy vs Predecessor Enegry",100,0,10,100,0,10);
	h->SetXTitle("Primary Energy [MeV]");
	h->SetYTitle("Predecessor Energy [MeV]");
	h2->SetLineColor(7);
	h3->SetLineColor(6);


	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k)//loop over ntuples/root files
    {	
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("t360");
		nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		//get the necessary information from ntuple
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Enekp",&ene);
		Td->SetBranchAddress("Enerau",&eneu);
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
			totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			beta14=beta1+4*beta4;
			h7->Fill(rf);
			h8->Fill(TMath::Sqrt(xe*xe+ye*ye+ze*ze));
			//set the damn flag for this event (UNIDOC)
			
			if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)& 0x6FE1)==0){
				damnFlag1 =1;  
			}
			else {
				damnFlag1=1;
			}
			//establish a clean flag
			if (rf<fiducialVolume  && damnFlag1==1 &&ene>emin&& itr<itrmax && itr>itrmin &&beta14<betamax && beta14>betamin) {
				cleanFlag1=1;
			}
			else {
				cleanFlag1=0;
			}
			//if event is clean start a lookback
			if (cleanFlag1==1 ){
				//cerr<<ene<<endl;
				//cerr<<evid<<" "<<lastID<<" "<<endl;
				//lastID=evid;
				numPi++;
				tTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
				numPred[numCoin]=0;
				ut1_i=ut1;
				ut2_i=ut2;
				ut3_i=ut3;
				jdy_i=jdy;
				tID=evid;
				tempNhits=nhits;
				tEn=ene;
				tR=TMath::Sqrt(xe*xe+ye*ye+ze*ze);
				xi=xfit;
				yi=yfit;
				zi=zfit;
				ri=rfit;
				l=lastClean+1;
				checkR=0;
				twoFlag=0;
				threeFlag=0;
				
				int lbEnd=0;
				if (j+20<nd) {
					lbEnd=j+20;
				}
				else {
					lbEnd=nd-1;
				}
				rut1=ut1;
				rut2=ut2;
				rut3=ut3;
				//start a lookback
				while (l<=lbEnd) {
					Td->GetEntry(l);
					//get the time difference
					day_diff=jdy_i-jdy;
					beta14=beta1+4*beta4;
					if (TMath::Abs(Int_t(ut1)-Int_t(rut1))==0) {
						rtrigTime=TMath::Abs(ut2-rut2)+TMath::Abs((ut3-rut3)/1.E3);
						/*
						if(rtrigTime<5){
							cerr<<"retrigger time "<<rtrigTime<<" "<<ut3<<" "<<rut3<<endl;
						}
						*/
					}
					else if (TMath::Abs(Int_t(ut1)-Int_t(rut1))==1) {
						//rtrigTime=TMath::Abs(ut2+(1-rut2))+TMath::Abs((ut3+1-rut3)/1.E3);
						if (ut1-rut1==1) {
							rtrigTime=ut2+(1.E6-rut2)+(ut3+1.E9-rut3)/1.E3;
							
						}
						else if (Int_t(ut1)-Int_t(rut1)==-1) {
							rtrigTime=rut2+(1.E6-ut2)+(rut3+1.E9-ut3)/1.E3;
							
						}
						else {
							rtrigTime=999;
						}

						/*
						if(rtrigTime<5 ){
							cerr<<"retrigger time "<<rtrigTime<<" "<<ut3<<" "<<rut3<<endl;
						}
						*/
					}
					else {
						rtrigTime=999;
					}

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
						else timediff=10;
						
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
					//if the time difference is within the window, check to see if it is clean
					if (TMath::Abs(timediff)<timeWindow && timediff>.000005 && rtrigTime>5) {
						//get the second damnFlag
						if(rtrigTime<10){
							cerr<<rtrigTime<<endl;
						}
						if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)& 0x6FE1)==0){
							damnFlag2 =1;  
						}
						else {
							damnFlag2=1;
						}
						//check if event is clean
						//&& itr<itrmax && itr>itrmin &&beta14<betamax && beta14>betamin && damnFlag2==1 
						if ( ene>ethresh && rf<fiducialVolume&& itr<itrmax && itr>itrmin &&beta14<betamax && beta14>betamin && damnFlag2==1 ) {
							//cerr<<"pred found for ev "<<tempID<<" it was ev "<<evid<<" at time "<<(ut1)+(ut2)/1.E6+(ut3)/1.E9<<" time diff is "<<timediff<<endl;
							if (timediff>.150) {
								cerr<<evid<<" "<<tID<<" "<<(ut1_i)+(ut2_i)/1.E6+(ut3_i)/1.E9<<" "<<(ut1)+(ut2)/1.E6+(ut3)/1.E9<<" "<<tR<<" "<<TMath::Sqrt(xe*xe+ye*ye+ze*ze)<<endl;
							}
							if (numPred[numCoin]==0) {
								
								evTime[numPred[numCoin]]=(ut1_i)+(ut2_i)/1.E6+(ut3_i)/1.E9;
								predX[numPred[numCoin]]=xi;
								predY[numPred[numCoin]]=yi;
								predZ[numPred[numCoin]]=zi;
								predR[numPred[numCoin]]=tR;
								predID[numPred[numCoin]]=tID;
								predEn[numPred[numCoin]]=tEn;
								predDt[numPred[numCoin]]=0;
								tempID[ct]=tID;
								tempEn[ct]=tEn;
								tempDt[ct]=0;
								tempTime[ct]=tTime;
								tempID[ct+1]=evid;
								tempEn[ct+1]=ene;
								tempDt[ct+1]=timediff;
								tempTime[ct+1]=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
								ct+=2;
								timeOfLast[numCoin][0]=jdy;
								timeOfLast[numCoin][1]=ut1;
								timeOfLast[numCoin][2]=ut2;
								timeOfLast[numCoin][3]=ut3;
								//h2->Fill(tempEn);
								
							}
							if (rtrigTime<10) {
								cerr<<rtrigTime<<" "<<evid<<endl;
							}
							evTime[numPred[numCoin]+1]=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
							predX[numPred[numCoin]+1]=xfit;
							predY[numPred[numCoin]+1]=yfit;
							predZ[numPred[numCoin]+1]=zfit;
							predR[numPred[numCoin]+1]=TMath::Sqrt(xe*xe+ye*ye+ze*ze);
							predID[numPred[numCoin]+1]=evid;
							predEn[numPred[numCoin]+1]=ene;
							predDt[numPred[numCoin]+1]=timediff;
							checkR=1;
							numPred[numCoin]++;
							lastClean=l;
							h1->Fill(timediff);
							timeOfFirst[numCoin][0]=jdy;
							timeOfFirst[numCoin][1]=ut1;
							timeOfFirst[numCoin][2]=ut2;
							timeOfFirst[numCoin][3]=ut3;
							//h3->Fill(ene);
							
						}
					
					}
					l++;
					rut1=ut1;
					rut2=ut2;
					rut3=ut3;
				}//end lookback
				//lbEnd=l;
				//lastClean=lbEnd-1;
				//last check to see if the coincidences make the rcut
				if(checkR==1){
					twoFlag=0;
					threeFlag=0;
					if (numPred[numCoin]==0) {
						twoFlag=0;
						threeFlag=0;
						
					}
					else if	(numPred[numCoin]==1){
						numOnePred++;
						deltaR[0]=TMath::Sqrt((predX[1]-xi)*(predX[1]-xi)+(predY[1]-yi)*(predY[1]-yi)+(predZ[1]-zi)*(predZ[1]-zi));
						deltaR[1]=0;
						deltaR[2]=0;
						tempR=deltaR[0];
						piEn=predEn[0];
						piTime=evTime[0];
						pe1En=predEn[1];
						pe1Time=evTime[1];
						piID=predID[0];
						pe1ID=predID[1];
						twoFlag=1;
						if (deltaR[0]<distanceWindow) {
							check2[numCoin]=1;
							twoFlag=1;	
							if (predDt[0] >.075 || predDt[1]>.0750 || predDt[2]>.0750) {
							cerr<<"2-Fold"<<endl;
							cerr<<predID[0]<<" "<<predEn[0]<<" "<<predR[0]<<" "<<evTime[0]<<endl;
							cerr<<predID[1]<<" "<<predEn[1]<<" "<<predR[1]<<" "<<evTime[1]<<endl;
							}
							h4->Fill(predID[0]);
							h4->Fill(predID[1]);
							h->Fill(predEn[0],predEn[1]);
						}
						h6->Fill(deltaR[0]);
						//cerr<<deltaR[0]<<" "<<xi<<" "<<predX[0]<<endl;

					}
					else if(numPred[numCoin]==2){
						numTwoPred++;
						threeFlag=0;
						twoFlag=0;
						deltaR[0]=TMath::Sqrt((predX[1]-xi)*(predX[1]-xi)+(predY[1]-yi)*(predY[1]-yi)+(predZ[1]-zi)*(predZ[1]-zi));
						deltaR[1]=TMath::Sqrt((predX[2]-xi)*(predX[2]-xi)+(predY[2]-yi)*(predY[2]-yi)+(predZ[2]-zi)*(predZ[2]-zi));
						deltaR[2]=TMath::Sqrt((predX[1]-predX[2])*(predX[1]-predX[2])+(predY[1]-predY[2])*(predY[1]-predY[2])+(predZ[1]-predZ[2])*(predZ[1]-predZ[2])); 
						h6->Fill(deltaR[0]);
						h6->Fill(deltaR[1]);
						h6->Fill(deltaR[2]);
						//cerr<<deltaR[0]<<" "<<deltaR[1]<<" "<<deltaR[2]<<endl;
						
						//case a
						if(deltaR[0]<distanceWindow &&deltaR[1]<distanceWindow &&deltaR[2]<distanceWindow){
							coin_events[numCoin]=numPred[numCoin];
							threeFlag=1;
							twoFlag=0;
							if (predDt[0] >.0750 || predDt[1]>.0750 || predDt[2]>.0750) {
								cerr<<"3-Fold"<<endl;
								cerr<<predID[0]<<" "<<predEn[0]<<" "<<predR[0]<<" "<<evTime[0]<<endl;
								cerr<<predID[1]<<" "<<predEn[1]<<" "<<predR[1]<<" "<<evTime[1]<<endl;
								cerr<<predID[2]<<" "<<predEn[2]<<" "<<predR[2]<<" "<<evTime[2]<<endl;
							}
							h4->Fill(predID[0]);
							h4->Fill(predID[1]);
							h4->Fill(predID[2]);
							h->Fill(predEn[0],predEn[1]);
							h->Fill(predEn[0],predEn[2]);

						}
						//case e
						else if (deltaR[0]<distanceWindow &&deltaR[1]<distanceWindow &&deltaR[2]>distanceWindow ){
							coin_events[numCoin]=numPred[numCoin];
							threeFlag=1;
							twoFlag=0;
							//cerr<<"3-Fold"<<endl;
							//cerr<<predID[0]<<" "<<predEn[0]<<" "<<evTime[0]<<endl;
							//cerr<<predID[1]<<" "<<predEn[1]<<" "<<evTime[1]<<endl;
							//cerr<<predID[2]<<" "<<predEn[2]<<" "<<evTime[2]<<endl;
							h4->Fill(predID[0]);
							h4->Fill(predID[1]);
							h4->Fill(predID[2]);
							h->Fill(predEn[0],predEn[1]);
							h->Fill(predEn[0],predEn[2]);
						}
						//case f
						else if (deltaR[0]<distanceWindow &&deltaR[1]>distanceWindow &&deltaR[2]>distanceWindow ){
							coin_events[numCoin]=numPred[numCoin]-1;
							twoFlag=1;
							threeFlag=0;
							if (predDt[0] >.0750 || predDt[1]>.0750 || predDt[2]>.0750) {
								cerr<<"2-Fold"<<endl;
								cerr<<predID[0]<<" "<<predEn[0]<<" "<<predR[0]<<" "<<evTime[0]<<endl;
								cerr<<predID[1]<<" "<<predEn[1]<<" "<<predR[1]<<" "<<evTime[1]<<endl;ce
							}

							h4->Fill(predID[0]);
							h4->Fill(predID[1]);
							h->Fill(predEn[0],predEn[1]);
						}
						//case h
						else if (deltaR[0]>distanceWindow &&deltaR[1]>distanceWindow &&deltaR[2]<distanceWindow ){
							coin_events[numCoin]=numPred[numCoin]-2;
							twoFlag=1;
							threeFlag=0;
							if (predDt[0] >.0750 || predDt[1]>.0750 || predDt[2]>.0750) {
								cerr<<"2-Fold"<<endl;
								cerr<<predID[1]<<" "<<predEn[1]<<" "<<predR[0]<<" "<<evTime[1]<<endl;
								cerr<<predID[2]<<" "<<predEn[2]<<" "<<predR[1]<<" "<<evTime[2]<<endl;
							}
							h4->Fill(predID[2]);
							h4->Fill(predID[1]);
							h->Fill(predEn[1],predEn[2]);
						}
						//case b
						else if (deltaR[0]>distanceWindow &&deltaR[1]>distanceWindow &&deltaR[2]>distanceWindow ){
							coin_events[numCoin]=numPred[numCoin]-2;
							twoFlag=0;
							threeFlag=0;
						}
						//case c
						else if (deltaR[0]<distanceWindow &&deltaR[1]>distanceWindow &&deltaR[2]<distanceWindow){
							coin_events[numCoin]=numPred[numCoin]-1;
							twoFlag=1;
							threeFlag=0;
							if (predDt[0] >.0750 || predDt[1]>.0750 || predDt[2]>.0750) {
								cerr<<"2-Fold"<<endl;
								cerr<<predID[0]<<" "<<predEn[0]<<" "<<predR[0]<<" "<<evTime[0]<<endl;
								cerr<<predID[1]<<" "<<predEn[1]<<" "<<predR[1]<<" "<<evTime[1]<<endl;
							}
							h4->Fill(predID[0]);
							h4->Fill(predID[1]);
							h->Fill(predEn[0],predEn[1]);
						}
						//case g
						else if (deltaR[0]<distanceWindow &&deltaR[1]>distanceWindow &&deltaR[2]<distanceWindow ){
							coin_events[numCoin]=numPred[numCoin]-1;
							twoFlag=1;
							threeFlag=0;
							if (predDt[0] >.0750 || predDt[1]>.0750 || predDt[2]>.0750) {
								cerr<<"2-Fold"<<endl;
								cerr<<predID[0]<<" "<<predEn[0]<<" "<<predR[0]<<" "<<evTime[0]<<endl;
								cerr<<predID[2]<<" "<<predEn[2]<<" "<<predR[2]<<" "<<evTime[2]<<endl;
							}
							h4->Fill(predID[0]);
							h4->Fill(predID[2]);
							h->Fill(predEn[0],predEn[2]);
						}
						//case d
						else if (deltaR[0]>distanceWindow &&deltaR[1]<distanceWindow &&deltaR[2]<distanceWindow ){
							coin_events[numCoin]=numPred[numCoin]-1;
							twoFlag=1;
							threeFlag=0;
							if (predDt[0] >.0750 || predDt[1]>.0750 || predDt[2]>.0750) {
							cerr<<"2-Fold"<<endl;
							cerr<<predID[0]<<" "<<predEn[0]<<" "<<predR[0]<<" "<<evTime[0]<<endl;
							cerr<<predID[2]<<" "<<predEn[2]<<" "<<predR[2]<<" "<<evTime[2]<<endl;
							}
							h4->Fill(predID[0]);
							h4->Fill(predID[2]);
							h->Fill(predEn[0],predEn[2]);
						}
					}
					else if(numPred[numCoin]>2){
						cerr<<numPred[numCoin]<<" predecessor event"<<endl;
						for (int pr=0; pr<=numPred[numCoin]; pr++) {
							cerr<<predID[pr]<<" "<<predR[pr]<<" "<<predEn[pr]<<" "<<predDt[pr]<<endl;
						}
					}




				}
				if (twoFlag==1) {
					numTwoFold++;
				}
				else if (threeFlag==1) {
					numThreeFold++;
				}

				numCoin++;
				if (lastClean>j) {
					j=lastClean;
				}
				else {
					lastClean=j;
				}

				
				//nextlastID=lastID;
				//lastID=tID;
			}//end if statement that checks to lookback
			//lastID=tID;
		}//end loop over events 
		delete Td;
		delete f1;
		cerr<<"next ntuple"<<endl;
	}//end loop over ntuples

	cerr<<numTwoFold<<" two fold events "<<numOnePred<<endl;
	cerr<<numThreeFold<<" three fold events "<<numTwoPred<<endl;
	//h3->Draw();
	//h2->Draw("same");
	h4->Draw();
	h.Draw("colz");

	cerr<<numPi<<" primaries"<<endl;

}//end code