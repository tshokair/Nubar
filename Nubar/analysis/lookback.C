/*
 *  coin.C
 *  
 *
 *  Created by Tim Shokair on 3/29/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *
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
	Int_t numev;
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
	Float_t thetaij;
	Int_t numfiles;
	TString files[1000];
	Int_t nd;
	Int_t day_diff;
	Int_t eventID[100000][50]={0.0};
	Int_t eventTime[100000][50]={0.0};
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
	Int_t numRawEvents;
	Int_t cutFlag;
	Int_t evNumLast;
	Int_t coin_events;
	Float_t timeOfLast;
	Int_t l;
	
	//parameters used in cuts
	Float_t nhitCut=27.5;
	Float_t nhitMax=160;
	Float_t fiducialVolume=550;
	Float_t timeWindow=0.150;
	Float_t distanceWindow=600;
	Float_t thetamin=0.63;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=4.5;
	
	//readin the file names
	ifstream runList;
	runList.open("runlist.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	
	//initialize the counted variables
	numev=0;
	nump=0;
	numg1=0;
	numg2=0;
	numCoincidence=0;
	numRawEvents=0;
	
	//make histograms
	TH1F * h1= new TH1F("h1","Nhits",100,0,200);
	TH1F * h2= new TH1F("h2","time",100,0,100000);
	
	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k)//loop over ntuples/root files
    {	
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("h330");
		nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		
		//get the necessary information from ntuple
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Eneu",&ene);
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
		
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			numRawEvents++;
			//get the vector between this event and the last one
			rfit=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi));
			//get the day differce to calculate time difference correctly
			day_diff=Int_t(jdy)-Int_t(jdy_i);
			totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			//set the damn flag for this event (UNIDOC)
			/*
			if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
				damnflag =1;  
			}
			else {
				damnflag=0;
			}
			*/
			
			//Orrell DAMN mask
			if((Int_t(rdamn1)& 0x02b1)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
				damnflag =1;  
			}
			else {
				damnflag=0;
			}
			
			
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
			if (nhits>nhitCut && nhits<nhitMax && rf<fiducialVolume && thetaij>thetamin && thetaij<thetamax && itr>itrmin && itr<itrmax && damnflag==1 &&ene>emin) {
				cutFlag=1;
			}
			else {
				cutFlag=0;
			}
			//look for a vertex if cuts are passes and vertex flag has not been set
			if (cutFlag==1) {
				coin_events=0;
				eventID[nump][coin_events]=evid;
				eventTime[nump][coin_events]=totalTime;
				temp_pos=nhits;
				epos=ene;
				ut1_i=ut1;
				ut2_i=ut2;
				ut3_i=ut3;
				jdy_i=jdy;
				xi=xfit;
				yi=yfit;
				zi=zfit;
				isFirst=1;
				timediff=-.1;
				l=evNumLast;
				while(timediff<timeWindow) { //lookback at events
					//find the time difference from last event if:
					
					
					Td->GetEntry(l);
					day_diff=jdy_i-jdy;
					//daydiff=-1
					if(Int_t(day_diff)==-1){
						dt1=(86400-ut1)+ut1_i;
						if(Int_t(dt1)==0){
							dt2=ut2_i/1.E6-ut2/1.E6;
							dt3=ut3_i/1.E9-ut3/1.E9;
							timediff=dt1+dt2 +dt3;
							
						}
						else if(Int_t(dt1)==-1){
							dt2=(1-ut2/1.E6)+ut2_i/1.E6;
							dt3=(1-ut3/1.E9)+ut3_i/1.E9;
							timediff=dt2+dt3;
						}
						else timediff=1;
					}
					//daydiff=0
					else if(Int_t(day_diff)==0){
						dt1=ut1-ut1_i;
						if(Int_t(dt1)==0){
							dt2=ut2/1.E6-ut2_i/1.E6;
							dt3=ut3/1.E9-ut3_i/1.E9;
							timediff=dt1+dt2+dt3;
						}
						else if(Int_t(dt1)==-1){
							dt2=(1-ut2/1.E6)+ut2_i/1.E6;
							dt3=(1-ut3/1.E9)+ut3_i/1.E9;
							timediff=dt2+dt3;
						}
						else timediff =1;
					}
					//daydiff =1
					else if(Int_t(day_diff)==1){
						dt1=(86400-ut1_i)+ut1;
						if(Int_t(dt1)==0){
							dt2=ut2/1.E6-ut2_i/1.E6;
							dt3=ut3/1.E9-ut3_i/1.E9;
							timediff=dt1+dt2 +dt3;
							
						}
						else if(Int_t(dt1)==1){
							dt2=(1-ut2_i/1.E6)+ut2/1.E6;
							dt3=(1-ut3_i/1.E9)+ut3/1.E9;
							timediff=dt2+dt3;
						}
						else timediff=1;
						
						
					}
					cerr<<timediff<<" "<<l<<" "<<j<<endl;
					h1->Fill(nhits);
					h2->Fill(ene);
					//cerr<<timediff<<" "<<ut1<<" "<<ut1_i<<endl;
					eventID[nump][coin_events]=evid;
					eventTime[nump][coin_events]=timediff;
					coin_events++;
					numev++;
					l++;
				}//end lookback loop
				evNumLast=j;
				nump++;
				timeOfLast=totalTime;
			}
		
		}//end loop over events

	}//end loop over ntuples
		Int_t n=0;
	for (Int_t p=0; p<nump; p++) {
		n=0;
		while (eventID[p][n+1]!=0) {
			n++;
		}
		cerr<<"tagged event id is " <<eventID[p][0]<<" and has "<<n<<" predecessors in the time window "<<endl;
	}
		
		
}//end code

