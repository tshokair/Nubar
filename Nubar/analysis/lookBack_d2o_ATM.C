/*
 *  lookBack_salt_out.h
 *  
 *
 *  Created by Tim Shokair on 11/8/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *  Code will be used for a mock analysis, rewrite of the look back search.
 *  Looks at the light water and acrylic events 
 *  12/9 fixed problem with duplicates two fold events (problem is with the sorted root file,
 *  but the analysis code was fixed to make sure no coincidences occured with different event numbers and the same gtid)
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
	Float_t ri;
	Float_t xfit;
	Float_t yfit;
	Float_t zfit;
	Float_t rfit;
	Float_t rf;
	Float_t dt;
	Float_t dt0;
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
	//temps
	Int_t lastClean;
	Int_t tempID[500000];
	Int_t tempNhits;
	double tempTime[500000];
	Float_t tempEn[500000];
	Int_t tID;
	Float_t tEn;
	double tTime;
	Int_t l=0;
	Float_t predX[100];
	Float_t predY[100];
	Float_t predZ[100];
	Int_t predID[100];
	Float_t predEn[100];
	double evTime[100];
	Float_t deltaR[3];
	Float_t tempR;
	Float_t piEn;
	Float_t pe1En;
	Float_t pe2En;
	double piTime;
	double pe1Time;
	double pe2Time;
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
	Float_t emin=3.5;
	Float_t ethresh=2.7;
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
	TH1F * h1= new TH1F("h1","Time Differences D2O",300,0,.1);
	TH1F * h2= new TH1F("h2","Primary Energy",100,0,20);
	TH1F * h3= new TH1F("h3","Predicessor Energy",100,0,20);
	TH1F * h4= new TH1F("h4","Coin IDs",30000,0,30000);
	TH1F * h5= new TH1F("h5","Coin Ev nums",30000,0,30000);
	TH1F * h6= new TH1F("h6","Delta R",1300,0,1300);
	h2->SetLineColor(7);
	h3->SetLineColor(6);


	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k)//loop over ntuples/root files
    {	
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("t102");
		nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		//get the necessary information from ntuple
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Enep",&ene);
		Td->SetBranchAddress("Itrp",&itr);
		Td->SetBranchAddress("Ev_tid",&evid);	
		Td->SetBranchAddress("T0_dif",&dt0); //julian day
		Td->SetBranchAddress("T_dif",&dt); //time in seconds

		Td->SetBranchAddress("Rfp",&rf);
		Td->SetBranchAddress("Thetaijp",&thetaij);

		Td->SetBranchAddress("Beta_14",&beta14);
		Td->SetBranchAddress("Rdamn1",&rdamn1);
		Td->SetBranchAddress("Rdamn2",&rdamn2);
		lastClean=0;
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			//totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			//beta14=beta1+4*beta4;
			//set the damn flag for this event (UNIDOC)
			
			if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
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
				//cerr<<evid<<" "<<lastID<<" "<<endl;
				//lastID=evid;
				numPred[numCoin]=0;

				tID=evid;
				tempNhits=nhits;
				tEn=ene;

				l=lastClean;
				checkR=0;
				twoFlag=0;
				threeFlag=0;

				//start a lookback
				while (l<=j) {
					Td->GetEntry(l);
					cerr<<dt0<<endl;
					//get the time difference
					beta14=beta1+4*beta4;
										//if the time difference is within the window, check to see if it is clean
					if (dt0<timeWindow && dt0!=0) {
						//get the second damnFlag
						if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
							damnFlag2 =1;  
						}
						else {
							damnFlag2=1;
						}
						//check if event is clean
						if (damnFlag2==1 && ene>ethresh && rf<fiducialVolume && itr<itrmax && itr>itrmin &&beta14<betamax && beta14>betamin) {
							//cerr<<"pred found for ev "<<tempID<<" it was ev "<<evid<<" at time "<<(ut1)+(ut2)/1.E6+(ut3)/1.E9<<" time diff is "<<timediff<<endl;
							if (numPred[numCoin]==0) {
								

								predID[numPred[numCoin]]=tID;
								predEn[numPred[numCoin]]=tEn;
								tempID[ct]=tID;
								tempEn[ct]=tEn;
								tempID[ct+1]=evid;
								tempEn[ct+1]=ene;
								ct+=2;
								//h2->Fill(tempEn);
								
							}
							predX[numPred[numCoin]+1]=xfit;
							predY[numPred[numCoin]+1]=yfit;
							predZ[numPred[numCoin]+1]=zfit;
							predID[numPred[numCoin]+1]=evid;
							predEn[numPred[numCoin]+1]=ene;
							checkR=1;
							numPred[numCoin]++;
							//lastClean=l;
							//h1->Fill(timediff);
							//h3->Fill(ene);
							
						}
					
					}
					l++;
				}//end lookback
				lastClean=j;
				//last check to see if the coincidences make the rcut
				if(checkR==1){
					if (numPred[numCoin]==0) {
						twoFlag=0;
						threeFlag=0;
					
					}
					else if	(numPred[numCoin]==1){
						piEn=predEn[0];
						pe1En=predEn[1];
						twoFlag=1;
						check2[numCoin]=1;
						
						//cerr<<deltaR[0]<<" "<<xi<<" "<<predX[0]<<endl;
					}
					else if(numPred[numCoin]==2){
						threeFlag=1;
						twoFlag=0;

					}
					else if(numPred[numCoin]>2){
						cerr<<numPred[numCoin]<<" predecessor event"<<endl;
					}




				}
				if (numCoin>0) {
					duplicateFlag=0;
					duplicateFlag2=0;

					int dupN;

					if ((ct-10)<0){
						dupN=0;
					}
					else {
						dupN=ct-10;
					}

					for (int n=dupN; n<ct-2; n++) {
						if(evid==tempID[n]){
							duplicateFlag=1;
							//cerr<<ct<<" -coin , id= "<<evid<<" is duplicate with j= "<< n<<" id = "<<tempID[n]<<endl;
						}
					}
					for (int n=dupN; n<ct; n++) {
						if(evid==tempID[n]){
							duplicateFlag2=1;
							//cerr<<ct<<" -coin , id= "<<evid<<" is duplicate with j= "<< n<<" id = "<<tempID[n]<<endl;
						}
					}
					if(threeFlag==1 && twoFlag==0 && duplicateFlag==0 ){
						numThreeFold++;
						h2->Fill(predEn[0]);
						h3->Fill(predEn[1]);
						h3->Fill(predEn[2]);

						cerr<<"Three Fold "<<endl;
						//cerr<<"dr0="<<deltaR[0]<<" "<<" dr1= "<<deltaR[1]<<" dr2= "<<deltaR[2]<<endl;
						for (int m=0; m<=numPred[numCoin]; m++) {
							h4->Fill(predID[m]);
							cerr<<predID[m]<<" "<<predEn[m]<<endl;
						}
						
						
					}
					else if (threeFlag==0 && check2[numCoin-1]==1 && duplicateFlag2==0  ) {
						numTwoFold++;
						h2->Fill(tempEn[ct-1]);
						h3->Fill(tempEn[ct-2]);
						h4->Fill(tempID[ct-1]);
						h4->Fill(tempID[ct-2]);
						//cerr<<tID<<" "<<lastID<<" "<<nextlastID<<endl;
						/*
						cerr<<"Two fold"<<endl;
						cerr<<" j: "<<j-1<<" "<<tempID[ct-1]<<" "<<tempTime[ct-1]<<" "<<tempEn[ct-1]<<endl;
						cerr<<" j: "<<j<<" "<<tempID[ct-2]<<" "<<tempTime[ct-2]<<" "<<tempEn[ct-2]<<endl;
						*/
						
					}
					

				}
				numCoin++;
				//nextlastID=lastID;
				//lastID=tID;
			}//end if statement that checks to lookback
			//lastID=tID;
		}//end loop over events 
		delete Td;
		delete f1;
		cerr<<"next ntuple"<<endl;
	}//end loop over ntuples

	cerr<<numTwoFold<<" two fold events"<<endl;
	cerr<<numThreeFold<<" three fold events"<<endl;
	//h3->Draw();
	//h2->Draw("same");
		h6->Draw();

}//end code