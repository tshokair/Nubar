/*
 *  lb_lw_salt.C
 *  
 *
 *  Created by Tim Shokair on 9/7/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *
 */
//start code
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
	Float_t xi;
	Float_t yi;
	Float_t zi;
	Float_t ri;
	Float_t xfit;
	Float_t yfit;
	Float_t zfit;
	Float_t rfit;
	Float_t rf;
	Float_t xgen;
	Float_t ygen;
	Float_t zgen;
	Float_t rgen;
	Int_t nd;
	Float_t evid;
	Float_t rdamn1;
	Float_t rdamn2;
	Float_t itr;
	Float_t thetaij;
	//variables for time differences
	Float_t dt1;
	Float_t dt2;
	Float_t dt3;
	Int_t day_diff;	
	Int_t totalDays;
	Double_t totalTime;
	Double_t runTime;
	Double_t liveTime;
	Double_t startTime;
	Double_t endTime;
	Float_t startDay;
	Float_t endDay;
	Double_t timediff;
	
	//flags
	Int_t damnflag;
	Int_t damnflag2;
	Int_t cutFlag;
	int neutronFlag;
	
	//lookback information
	Float_t predX[100]={0.0};
	Float_t predY[100]={0.0};
	Float_t predZ[100]={0.0};
	Float_t predR[100]={0.0};
	Float_t deltaR[3];
	Double_t timeOfLast=0;
	Int_t evNumLast=0;
	int preds;
	
	//counted things
	Int_t numRawEvents=0;
	int totalClean=0;
	int numN=0;
	int numPred[100000]={0.0};
	int numPos=0;
	int coins=0;
	int totalCoins=0;
	//temps
	int tempEventID=0;
	int neutronEn;
	Double_t deltaT;	
	//parameters used in cuts
	Float_t nhitCut=20;
	Float_t nhitMax=160;
	Float_t fiducialVolume=600;
	Float_t timeWindow=0.250;
	Float_t distanceWindow=1000;
	Float_t startWindow=100;
	Float_t thetamin=0.63;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=3.5;
	Float_t ethresh=2.7;
	Float_t outerVolume=650;
	Float_t avIn=600;
	Float_t avOut=606;
	Int_t rstep=25;
	Int_t numstep=8;
	//files info
	Int_t numfiles=0;
	TString files[1000];

	
	//readin the file names
	ifstream runList;
	runList.open("test.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	
	
	//make histograms
	TH1F * h1= new TH1F("h1","DeltaR two fold",100,0,500);
	TH1F * h2= new TH1F("h2","DeltaR two fold (1 in 1 out)",100,0,500);	
	TH1F * h3= new TH1F("h3","Neutron Energy",100,0,20);
	TH1F * h5= new TH1F("h5","Positron Energy",100,0,20);
	TH1F * h4= new TH1F("h4","Energy",100,0,20);


	h2->SetLineColor(6);	
	
	
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
		Td->SetBranchAddress("Itrp",&itr);
		Td->SetBranchAddress("Ev_tid",&evid);	
		Td->SetBranchAddress("Time_jdy",&jdy); //julian day
		Td->SetBranchAddress("Time_s",&ut1); //time in seconds
		Td->SetBranchAddress("Time_us",&ut2); //time in nanoseconds
		Td->SetBranchAddress("Time_ns",&ut3); //time in nanoseconds
		Td->SetBranchAddress("Xfp",&xfit);
		Td->SetBranchAddress("Yfp",&yfit);
		Td->SetBranchAddress("Zfp",&zfit);
		Td->SetBranchAddress("Xe",&xgen);
		Td->SetBranchAddress("Ye",&ygen);
		Td->SetBranchAddress("Ze",&zgen);
		Td->SetBranchAddress("Rfp",&rf);
		Td->SetBranchAddress("Thetaijp",&thetaij);
		Td->SetBranchAddress("Rdamn1",&rdamn1);
		Td->SetBranchAddress("Rdamn2",&rdamn2);
		Td->GetEntry(0);
		//initialize for first event
		jdy_i=jdy;
		ut1_i=0;
		ut2_i=0;
		ut3_i=0;
		xi=0;
		yi=0;
		zi=0;
		timeOfLast=0;
		evNumLast=0;
		
		//loop over events in file
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			numRawEvents++;
			//tid=evid;
			//get the vector between this event and the last one
			rfit=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi));
			//get the day differce to calculate time difference correctly
			day_diff=Int_t(jdy)-Int_t(jdy_i);
			totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			neutronFlag=0;	
			if (ene!=-9999) {
				//cerr<<"e= "<<ene<<endl;
				h4->Fill(ene);
			}
			
			//set the damn flag for this event (UNIDOC)
	
			if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
				damnflag =1;  
			}
			else {
				damnflag=0;
			}
			if (rf<fiducialVolume  && damnflag==1 &&ene>emin) {
				cutFlag=1;
				
				h3->Fill(ene);
			}
			else if(rf<outerVolume && damnflag==1 && ene>ethresh){
				cutFlag=0;
				totalClean++;
				neutronFlag=1;
				//cerr<<"energy: "<<ene<<endl;

			}
			else {
				cutFlag=0;
			}
			//if event makes the cuts start a lookback
			if (cutFlag==1){
				
				numPred[numN]=0;
				tempEventID=evid;
				ut1_i=ut1;
				ut2_i=ut2;
				ut3_i=ut3;
				jdy_i=jdy;
				xi=xfit;
				yi=yfit;
				zi=zfit;
				ri=rf;
				//lookback
				int l=evNumLast;
				//cerr<<"starting a lookback"<<endl;
				while (l<j) {
					Td->GetEntry(l);
					if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
						damnflag2 =1;  
					}
					else {
						damnflag2=0;
					}
					day_diff=jdy_i-jdy;
					rfit=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi));
					//if cleanget the time diff
					if (damnflag2==1) {
						
						
						//if dayfiff=-1
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
							else timediff=999;
							
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
							else timediff =999;
							
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
						
					}
					if (timediff<timeWindow && rf<outerVolume && rf>fiducialVolume && ene>ethresh) {
						predX[numPred[numN]]=xfit;
						predY[numPred[numN]]=yfit;
						predZ[numPred[numN]]=zfit;
						predR[numPred[numN]]=rf;
						numPred[numN]++;
						deltaT=timediff;
						
					}
					l++;
				
				}//end lookback
				if (numPred[numN]>0) {
					coins++;
					cerr<<"Particle "<<tempEventID<<" had "<<numPred[numN]<<" predecessors"<<endl;
					if(numPred[numN]==1){
						rgen=TMath::Sqrt(xgen*xgen+ygen*ygen+zgen*zgen);
						deltaR[0]=TMath::Sqrt((predX[0]-xi)*(predX[0]-xi)+(predY[0]-yi)*(predY[0]-yi)+(predZ[0]-zi)*(predZ[0]-zi));
						cerr<<"delta r= "<<deltaR[0]<<" cm from ri ="<<rgen<<" to rf= "<<ri<<" in a time of "<<deltaT<<endl;
						h1->Fill(deltaR[0]);
					}
				}
				numN++;
				
			}
			if (neutronFlag==1){
				preds=0;
				//lookback
				int l=evNumLast;
				//cerr<<"starting a lookback"<<endl;
				while (l<j) {
					Td->GetEntry(l);
					if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
						damnflag2 =1;  
					}
					else {
						damnflag2=0;
					}
					day_diff=jdy_i-jdy;

					//if cleanget the time diff
					if (damnflag2==1) {
						
						
						//if dayfiff=-1
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
							else timediff=999;
							
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
							else timediff =999;
							
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
						
					}
					if (timediff<timeWindow && ene>ethresh) {
						rfit=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi));
						preds++;	
						h2->Fill(rfit);
						h5->Fill(ene);
					}
					l++;
					
				}//end lookback
				if (preds>0) {
					totalCoins++;
					//cerr<<"In lw, ev "<<j<<" had "<<preds <<" predecessors -"<<totalCoins<<endl;
				}
				evNumLast=j;
			}
			
		}//end loop over events
		
		
	}//end loop over ntuples
		cerr<<coins<<" coincidences found inside "<< totalCoins<<" found total from "<<totalClean<<" events"<< endl;
		//h2->Draw();
		//h1->Draw("same");
		h4->Draw();
}//end code