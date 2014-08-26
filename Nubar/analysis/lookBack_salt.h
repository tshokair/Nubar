/*
 *  lookBack_salt.h
 *  
 *
 *  Created by Tim Shokair on 11/8/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *  Code will be used for a mock analysis, rewrite of the look back search.
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
	//counts;
	Int_t numCoin=0;
	Int_t numPred[10000]={0.0};
	Int_t coin_events[10000]={0.0};
	//flags
	Int_t damnFlag1=0;
	Int_t damnFlag2=0;
	Int_t cleanFlag1=0;
	Int_t cleanFlag2=0;
	Int_t checkR=0;
	Int_t twoFlag=0;
	Int_t threeFlag=0;
	//times
	Int_t totalTime;
	//temps
	Int_t lastClean;
	Int_t tempID;
	Int_t tempNhits;
	Float_t tempEn;
	Int_t l=0;
	Float_t PredX[100];
	Float_t PredY[100];
	Float_t PredZ[100];
	Int_t PredID[100];
	double evTime[100];
	Float_t deltaR[3];

	//parameters
	Float_t nhitCut=20;
	Float_t nhitMax=160;
	Float_t fiducialVolume=600;
	Float_t timeWindow=0.200;
	Float_t distanceWindow=300;
	Float_t startWindow=100;
	Float_t thetamin=0.63;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=3.5;
	Float_t ethresh=2.7;
	Float_t outerVolume=600;
	
	//readin the file names
	ifstream runList;
	runList.open("test.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	
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
		lastClean=0;
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			//set the damn flag for this event (UNIDOC)
			
			if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
				damnflag1 =1;  
			}
			else {
				damnflag1=0;
			}
			//establish a clean flag
			if (rf<fiducialVolume  && damnflag==1 &&ene>emin) {
				cleanFlag1=1;
			}
			else {
				cleanFlag1=0;
			}
			//if event is clean start a lookback
			if (cleanFlag1==1){
				numPred[numCoin]=0;
				ut1_i=ut1;
				ut2_i=ut2;
				ut3_i=ut3;
				jdy_i=jdy;
				tempID=evid;
				tempNhits=nhits;
				tempEn=ene;
				xi=xfit;
				yi=yfit;
				zi=zfit;
				ri=rfit;
				l=lastClean;
				checkR=0;
				//start a lookback
				while (l<=j) {
					Td->GetEntry(l);
					//get the time difference
					day_diff=jdy_i-jdy;
					
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
					if (timediff<timeWindow && timediff!=0) {
						
						//get the second damnFlag
						if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
							damnflag2 =1;  
						}
						else {
							damnflag2=0;
						}
						//check if event is clean
						if (damnflag2==1 && ene>ethresh && rfit<outerVolume) {
							evTime[numPred[numCoin]]=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
							PredX[numPred[numCoin]]=xfit;
							PredY[numPred[numCoin]]=yfit;
							PredZ[numPred[numCoin]]=zfit;
							PredID[numPred[numCoin]]=evid;
							checkR=1;
							numPred[numCoin]++;
							lastClean=l;
						}
					}
					l++;
				}//end lookback
				//last check to see if the coincidences make the rcut
				if(checkR==1){
					if (numPred[numCoin]==0) {
						twoFlag=0;
						threeFlag=0;
					}
					else if	(numPred[numCoin]==1){
						deltaR[0]=TMath::Sqrt((predX[0]-xi)*(predX[0]-xi)+(predY[0]-yi)*(predY[0]-yi)+(predZ[0]-zi)*(predZ[0]-zi));
						deltaR[1]=0;
						deltaR[2]=0;
						twoflag=1;
					}
					else if(numPred[numCoin]==2){
						threeFlag=1;
						twoflag=0;
						deltaR[0]=TMath::Sqrt((predX[0]-xi)*(predX[0]-xi)+(predY[0]-yi)*(predY[0]-yi)+(predZ[0]-zi)*(predZ[0]-zi));
						deltaR[1]=TMath::Sqrt((predX[1]-xi)*(predX[1]-xi)+(predY[1]-yi)*(predY[1]-yi)+(predZ[1]-zi)*(predZ[1]-zi));
						deltaR[2]=TMath::Sqrt((predX[0]-predX[1])*(predX[0]-predX[1])+(predY[0]-predY[1])*(predY[0]-predY[1])+(predZ[0]-predZ[1])*(predZ[0]-predZ[1])); 
					}
					
					if(deltaR[0]<distanceWindow &&deltaR[1]<distanceWindow &&deltaR[2]<distanceWindow){
						coin_events[numCoin]=numPred[numCoin];
					}
					else if (deltaR[0]<distanceWindow &&deltaR[1]<distanceWindow &&deltaR[2]>distanceWindow ){
						coin_events[numCoin]=numPred[numCoin];
					}
					else if (deltaR[0]<distanceWindow &&deltaR[1]>distanceWindow &&deltaR[2]>distanceWindow ){
						coin_events[numCoin]=numPred[numCoin]-1;
					}
					else if (deltaR[0]>distanceWindow &&deltaR[1]>distanceWindow &&deltaR[2]<distanceWindow ){
						coin_events[numCoin]=numPred[numCoin]-2;
					}
					else if (deltaR[0]>distanceWindow &&deltaR[1]>distanceWindow &&deltaR[2]>distanceWindow ){
						coin_events[numCoin]=numPred[numCoin]-2;
					}
					else if (deltaR[0]<distanceWindow &&deltaR[1]>distanceWindow &&deltaR[2]>distanceWindow){
						coin_events[numCoin]=numPred[numCoin]-1;
					}
					else if (deltaR[0]<distanceWindow &&deltaR[1]>distanceWindow &&deltaR[2]<distanceWindow ){
						coin_events[numCoin]=numPred[numCoin]-1;
					}
					else if (deltaR[0]>distanceWindow &&deltaR[1]<distanceWindow &&deltaR[2]<distanceWindow ){
						coin_events[numCoin]=numPred[numCoin]-1;
					}
				}
				if (numPred[numCoin]>0) {
					for (int m=0; m<=numPred[numCoin]; m++) {
						cerr<<PredID[m]<<" "<<evTime[m]<<endl;
					}
				}
				numCoin++;
			}//end if statement that checks to lookback
			
		}//end loop over events 
		
		
	}//end loop over ntuples
	

}//end code