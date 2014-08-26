/*
 *  lookback at Cf data to determine see how the coinicidence rates change with the r cut.
 *  
 *
 *  Created by Tim Shokair on 3/29/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *
 *  Modified 6/2/11
 Modified 6/8/11 to include a weighted mean of distance for rcut
 Modified 7/29/11 to look at SNOMAN Data
 8/18/11 changed for Salt 
 
 */
{
	//declare variables
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
	Int_t numev[10];
	Int_t nump;
	Int_t numg1;
	Int_t numg2;
	Int_t is_pos;
	Int_t is_gamma;
	Float_t xi;
	Float_t yi;
	Float_t zi;
	Float_t ri;
	Float_t xfit;
	Float_t yfit;
	Float_t zfit;
	Float_t rfit;
	Float_t rf;
	Float_t rmag;
	Float_t dt1;
	Float_t dt2;
	Float_t dt3;
	Float_t temp_pos;
	Float_t epos;
	Float_t neutronEn;
	Float_t thetaij;
	Int_t numfiles;
	TString files[1000];
	Int_t nd;
	Int_t day_diff;
	Int_t eventID[500000][100]={0.0};
	Int_t eventTime[500000][100]={0.0};
	Float_t evid;
	Int_t partNum;
	Int_t tempEventID;
	Int_t isFirst;
	Int_t isSecond;
	Int_t isThird;
	Double_t tempTime;
	Int_t numCoincidence;
	Double_t totalTime;	
	Double_t totalTime_i;
	Double_t timediff;
	Float_t rdamn1;
	Float_t rdamn2;
	Int_t damnflag;
	Int_t damnflag2;
	Int_t numRawEvents;
	Int_t cutFlag;
	Int_t evNumLast;
	Int_t coin_events[1000000][50]={0.0};
	Float_t timeOfLast;
	Int_t l;
	Int_t lLast;
	Int_t neutronFlag;
	Float_t neutronTheta;
	Int_t totalDays;
	Double_t runTime;
	Double_t liveTime;
	Double_t startTime;
	Double_t endTime;
	Float_t startDay;
	Float_t endDay;
	Int_t num2fold[50];
	Int_t num3fold[50];
	Int_t num4fold[50];
	Int_t numPred[1000000]={0.0};
	Int_t checkR;
	Float_t predX[100]={0.0};
	Float_t predY[100]={0.0};
	Float_t predZ[100]={0.0};
	Float_t deltaR[3];
	Int_t weight;
	Int_t rFlag;
	Int_t rstep;
	Int_t total2fold=0;
	Int_t total3fold=0;
	Float_t acc;
	Float_t acc2;
	Float_t acc3;
	Float_t maxR;
	Int_t tid;
	Int_t counted=0;
	//flags and temps
	Int_t twoflag=0;
	Int_t threeflag=0;
	Float_t rtemp=0;
	Float_t tempEnPi=0;
	Float_t tempEnPd=0;
	Float_t tempTheta=0;
	Float_t tempNhits=0;
	Float_t posEn;
	Float_t nEn;
	Int_t check2[1000000]={0.0};
	
	//parameters used in cuts
	Float_t nhitCut=20;
	Float_t nhitMax=160;
	Float_t fiducialVolume=600;
	Float_t timeWindow=0.050;
	Float_t distanceWindow=300;
	Float_t startWindow=100;
	Float_t thetamin=0.63;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=3.5;
	Float_t ethresh=2.7;
	Float_t outerVolume=620;
	Float_t avIn=600;
	Float_t avOut=606;
	Int_t rstep=25;
	Int_t numstep=8;
	//readin the file names
	ifstream runList;
	runList.open("mcsalt.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	
	//initialize the counted variables
	//numev[=0;
	nump=0;
	numg1=0;
	numg2=0;
	numCoincidence=0;
	numRawEvents=0;
	liveTime=0;
	runTime=0;
	Int_t real3fold=0;
	Int_t real2fold=0;
	//make histograms
	TH1F * h1= new TH1F("h1","Positron Energy",100,0,20);
	TH1F * h2= new TH1F("h2","time",300,0,3);
	TH1F * h3= new TH1F("h3","Thetaij Neutron",100,.3,1.8);
	TH1F * h4= new TH1F("h4","Nhit Neutron",150,0,150);
	TH1F * h5= new TH1F("h5","Neutron energy",100,0,20);
	TH1F * h6= new TH1F("h6","DeltaR two fold",2400,0,1200);
	TH1F * h7= new TH1F("h7","DeltaR three fold",2400,0,1200);
	TH1F * h8= new TH1F("h8","",400,0,300);
	TH1F * h9= new TH1F("h9","Acceptance Ratios ",400,0,300);
	TH1F * h10= new TH1F("h10","Acceptance Ratios ",400,0,300);
	TH1F * h11= new TH1F("h11","time",300,0,3);
	TH1F * h12= new TH1F("h12","DeltaR of Any Coincidence",2400,0,1200);
	TH1F * h13= new TH1F("h13","Neutron Time Differences",300,0,.3);

	h6->SetLineColor(5);
	h7->SetLineColor(6);
	h8->SetLineColor(9);
	h8->SetLineWidth(3);
	h9->SetLineColor(2);
	h9->SetLineWidth(3);
	h10->SetLineColor(8);
	h10->SetLineWidth(3);
	h12->SetLineColor(4);
	h8->GetXaxis()->SetTitle("R-Cut [cm]");
	h8->GetYaxis()->SetTitle("coincidences tagged/coincidences generated");
	nump=0;
	liveTime=0;
	runTime=0;
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
		Td->GetEntry(0);
		
		//initialize variables for first event
		isFirst=0;
		isSecond=0;
		isThird=0;
		jdy_i=jdy;
		ut1_i=0;
		ut2_i=0;
		ut3_i=0;
		xi=0;
		yi=0;
		zi=0;
		timeOfLast=0;
		evNumLast=0;
		l=0;
		startTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
		timeOfLast=startTime;
		startDay=jdy;
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			numRawEvents++;
			tid=evid;
			//get the vector between this event and the last one
			rfit=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi));
			//get the day differce to calculate time difference correctly
			day_diff=Int_t(jdy)-Int_t(jdy_i);
			totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;

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
			h11->Fill(timediff);
			if(timediff!=10){
				counted++;
			}
			jdy_i=jdy;
			ut1_i=ut1;
			ut2_i=ut2;
			ut3_i=ut3;
			//cerr<<"event "<<evid<<" started at "<<ut1<<"+"<<ut2<<endl;
			//set the damn flag for this event (UNIDOC)
			
			if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
				damnflag =1;  
			}
			else {
				damnflag=0;
			}
			
			/*
			 //Orrell DAMN mask
			 if((Int_t(rdamn1)& 0x02b1)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
			 damnflag =1;  
			 }
			 else {
			 damnflag=0;
			 }
			 */
			
			//Melissa DAMN mask
			/*
			 if((Int_t(rdamn1)& 0x047B)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
			 damnflag =1;  
			 }
			 else {
			 damnflag=0;
			 }
			 */
			//damnflag=1;
			
			
			//establish a flag for the cuts
			if (nhits>nhitCut && nhits<nhitMax && rf<fiducialVolume  && damnflag==1 &&ene>emin) {
				cutFlag=1;
			}
			else {
				cutFlag=0;
			}
			//look for a vertex if cuts are passes and vertex flag has not been set
			if (cutFlag==1) {
				numPred[nump]=0;
				coin_events[nump][0]=0;
				coin_events[nump][1]=0;
				coin_events[nump][2]=0;
				coin_events[nump][3]=0;
				coin_events[nump][4]=0;
				coin_events[nump][5]=0;
				tempEventID=evid;
				//eventID[nump][coin_events]=evid;
				//eventTime[nump][coin_events]=totalTime;
				temp_pos=nhits;
				neutronEn=ene;
				neutronTheta=thetaij;
				ut1_i=ut1;
				ut2_i=ut2;
				ut3_i=ut3;
				jdy_i=jdy;
				xi=xfit;
				yi=yfit;
				zi=zfit;
				ri=rf;
				isFirst=1;
				timediff=0;
				l=evNumLast;
				neutronFlag=0;
				checkR=0;
				threeflag=0;
				//cerr<<"starting lookback for event "<<evid<<" at time "<<ut1<<" + "<<ut2<<endl;
				
				while(l<=j) { //lookback at events
					//find the time difference from last event if:
					
					
					Td->GetEntry(l);
					if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
						damnflag2 =1;  
					}
					else {
						damnflag2=0;
					}
					day_diff=jdy_i-jdy;
					rfit=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi));
					//cerr<<day_diff<<" "<<endl;
					//daydiff=-1
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
					if(timediff!=10){
						h2->Fill(timediff);
					}
					if (numPred[nump]==1) {
						h13->Fill(timediff);
					}
					if (timediff<timeWindow &&damnflag2==1 && ene>ethresh &&rf<outerVolume) {
						if (numPred[nump]==0) {
							posEn=ene;
						}
						else if (numPred[nump]==1){
							nEn=ene;
						}
						
						
						eventID[nump][numPred[nump]]=evid;
						//eventTime[nump][coin_events+1]=timediff;
						predX[numPred[nump]]=xfit;
						predY[numPred[nump]]=yfit;
						predZ[numPred[nump]]=zfit;
						numPred[nump]++;
						//numev[nump]++;
						checkR=1;
						evNumLast=l;
					}
					
					l++;
				}//end lookback loop
				if(checkR==1){
					if (numPred[nump]==0) {
						twoflag=0;
						threeflag=0;
					}
					else if(numPred[nump]==1){
						deltaR[0]==TMath::Sqrt((predX[0]-xi)*(predX[0]-xi)+(predY[0]-yi)*(predY[0]-yi)+(predZ[0]-zi)*(predZ[0]-zi));
						deltaR[1]=0;
						deltaR[2]=0;
						twoflag=1;
						check2[nump]=1;
						rtemp=deltaR[0];
						tempTheta=neutronTheta;
						tempEnPi=neutronEn;
						tempEnPd=posEn;
						tempNhits=temp_pos;
						h6->Fill(deltaR[0]);
						h12->Fill(deltaR[0]);
						
					}
					else if (numPred[nump]==2){
						threeflag=1;
						twoflag=0;
						deltaR[0]=TMath::Sqrt((predX[0]-xi)*(predX[0]-xi)+(predY[0]-yi)*(predY[0]-yi)+(predZ[0]-zi)*(predZ[0]-zi));
						deltaR[1]=TMath::Sqrt((predX[1]-xi)*(predX[1]-xi)+(predY[1]-yi)*(predY[1]-yi)+(predZ[1]-zi)*(predZ[1]-zi));
						deltaR[2]=TMath::Sqrt((predX[0]-predX[1])*(predX[0]-predX[1])+(predY[0]-predY[1])*(predY[0]-predY[1])+(predZ[0]-predZ[1])*(predZ[0]-predZ[1])); 
						h7->Fill(deltaR[0]);
						h7->Fill(deltaR[1]);
						h7->Fill(deltaR[2]);
						h12->Fill(deltaR[0]);
						h12->Fill(deltaR[1]);
						h12->Fill(deltaR[2]);
						for(int p=0;p<=numstep;p++){
							if(deltaR[0]<distanceWindow-rstep*p &&deltaR[1]<distanceWindow-rstep*p &&deltaR[2]<distanceWindow-rstep*p ){
								coin_events[nump][p]=numPred[nump];
							}
							else if (deltaR[0]<distanceWindow-rstep*p &&deltaR[1]<distanceWindow-rstep*p &&deltaR[2]>distanceWindow-rstep*p ){
								coin_events[nump][p]=numPred[nump];
							}
							else if (deltaR[0]<distanceWindow-rstep*p &&deltaR[1]>distanceWindow-rstep*p &&deltaR[2]>distanceWindow-rstep*p ){
								coin_events[nump][p]=numPred[nump]-1;
							}
							else if (deltaR[0]>distanceWindow-rstep*p &&deltaR[1]>distanceWindow-rstep*p &&deltaR[2]<distanceWindow-rstep*p ){
								coin_events[nump][p]=numPred[nump]-2;
							}
							else if (deltaR[0]>distanceWindow-rstep*p &&deltaR[1]>distanceWindow-rstep*p &&deltaR[2]>distanceWindow-rstep*p ){
								coin_events[nump][p]=numPred[nump]-2;
							}
							else if (deltaR[0]<distanceWindow-rstep*p &&deltaR[1]>distanceWindow-rstep*p &&deltaR[2]>distanceWindow-rstep*p ){
								coin_events[nump][p]=numPred[nump]-1;
							}
							else if (deltaR[0]<distanceWindow-rstep*p &&deltaR[1]>distanceWindow-rstep*p &&deltaR[2]<distanceWindow-rstep*p ){
								coin_events[nump][p]=numPred[nump]-1;
							}
							else if (deltaR[0]>distanceWindow-rstep*p &&deltaR[1]<distanceWindow-rstep*p &&deltaR[2]<distanceWindow-rstep*p ){
								coin_events[nump][p]=numPred[nump]-1;
							}
						}
						total3fold++;
						neutronFlag=1;
					}
					
				}
				else if (checkR!=1){
					neutronFlag=0;
				}
				maxR=0;
				Int_t q=0;
				//if really a two fold fill in the info for the two fold
				if (nump>0) {
					if(check2[nump-1]==1 &&threeflag==0){
						//cerr<<"doing something for 2fold- "<<real2fold<<endl;
						real2fold++;
						total2fold++;
						neutronFlag=1;

						h4->Fill(tempNhits);
						h5->Fill(tempEnPi);
						h1->Fill(tempEnPd);
						h3->Fill(neutronTheta);
						for(int p=0;p<=numstep;p++){
							if(rtemp<distanceWindow-rstep*p){
								coin_events[nump][p]=1;
							}
						}	
					}
				}
 				
				if (neutronFlag==1 && threeflag==1 &&twoflag==0) {
					//cerr<<"doing something for 3fold- "<<total3fold<<endl;
					
					h5->Fill(nEn);
					h4->Fill(temp_pos);
					h5->Fill(neutronEn);
					h3->Fill(neutronTheta);
					h1->Fill(posEn);
					
				}
				//evNumLast=j;
				
				timeOfLast=totalTime;
			}
			/*
			 cerr<<"event "<<tid<<" had " <<coin_events[nump][0]<<" predecessors";
			 if(coin_events[nump][0]==1){
			 cerr<<": "<<eventID[nump][0]<<endl;
			 }
			 else if (coin_events[nump][0]==2){
			 real3fold++;
			 cerr<<": "<<eventID[nump][0]<<" and "<<eventID[nump][1]<<" for "<<real3fold<<" 3fold coincedences "<<endl;
			 
			 }
			 else {
			 cerr<<endl;
			 }
			 */
			
			endTime=totalTime;
			endDay=jdy;
			nump++;
			
		}//end loop over events
		delete Td;
		delete f1;
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
		}//end loop over ntuples
		cerr<<total3fold<<" real 3fold events and "<<real2fold<< " 2fold evens"<<endl;
		for(int s=0;s<numstep;s++){
			num2fold[s]=0;
			num3fold[s]=0;
			num4fold[s]=0;
			for (Int_t p=0; p<=nump; p++) {
				
				//cerr<<"coincidece " <<p<<" has "<<coin_events[p][s]<<" predecessors in the time window "<<endl;
				if(coin_events[p][s]==1){
					num2fold[s]++;
				}
				if(coin_events[p][s]==2){
					num3fold[s]++;
				}
				if(coin_events[p][s]>2){
					num4fold[s]++;
				}
			}
			
			acc2= (Float_t)num2fold[s]/total2fold; 
			acc3= (Float_t)num3fold[s]/total3fold;
			acc= (Float_t)(num2fold[s]+num3fold[s])/(total2fold +total3fold);
			h8->Fill(distanceWindow-rstep*s-2,acc2);
			h9->Fill(distanceWindow-rstep*s,acc3);
			h10->Fill(distanceWindow-rstep*s+2,acc);
			cerr<<"the rate of tagged events is at rcut of "<<distanceWindow-rstep*s<<" cm is "<<nump/liveTime<<" Hz"<<endl;
			cerr<<"acceptance of 2 fold coincidences at rcut of  "<<distanceWindow-rstep*s<<" cm is "<<num2fold[s]<<" / "<<total2fold<<" = "<<acc2<<" "<<endl;
			cerr<<"acceptance of 3 fold coincidences at rcut of  "<<distanceWindow-rstep*s<<" cm is "<<num3fold[s]<<"/"<<total3fold<<" = "<<acc3<<endl;
			cerr<<"the acceptance coincidences at "<<distanceWindow-rstep*s<<" cm is "<<acc<<endl;
			cerr<<endl;
		}
		cerr<<nump<<" total clean events"<<endl;
		h8->Draw();
		h9->Draw("same");
		h10->Draw("same");
		std::stringstream cap;
		cap<<"Acceptance Ratios for Coincidences in D2O ["<<total2fold+total3fold<<" MC coincidences total]";
		TString st=cap.str();
		leg= new TLegend(0.6,0.7,0.89,0.89);
		//win= new TLegend(0.2,0.2,0.89,0.1);
		win = new TLegend(0.01,0.9401613,0.71,0.995,"blNDC");
		win->SetHeader(st);
		win->Draw();
		leg->AddEntry(h8,"2-Fold", "l");
		leg->AddEntry(h9,"3-Fold", "l");
		leg->AddEntry(h10,"All", "l");
		leg->SetHeader("Coincidence Type");
		leg->Draw();
		cerr<<counted<<" of "<<numRawEvents<<" made the time cut"<<endl;
}//end code
		
