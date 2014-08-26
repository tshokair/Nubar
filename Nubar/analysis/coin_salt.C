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
	Int_t eventID[100000][3]={0.0};
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
	Float_t clock[10000000];
	Float_t totaldays=0;
	Float_t clockedtime;
	Float_t rate;
	
	//parameters used in cuts
	Float_t nhitCut=20;
	Float_t nhitMax=160;
	Float_t fiducialVolume=550;
	Float_t timeWindow=5;
	Float_t distanceWindow=1100;
	Float_t thetamin=0.63;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=5.0;
	
	//readin the file names
	ifstream runList;
	runList.open("runlist_salt.txt");
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
	timecounter=0;
	
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
		
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			numRawEvents++;
			
			//get the vector between this event and the last one
			rfit=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi));
			//get the day differce to calculate time difference correctly
			day_diff=Int_t(jdy)-Int_t(jdy_i);
			totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			if(Int_t(day_diff)==1){ totaldays++;}

			clock[j]=ut1+totaldays*84000;
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
			
			//find the time difference from last event if:
			//daydiff =1
			if(Int_t(day_diff)==1){
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
				else timediff=100;
				

			}
			//daydiff=0
			else if(Int_t(day_diff)==0){
				dt1=ut1-ut1_i;
				if(Int_t(dt1)==0){
					dt2=ut2/1.E6-ut2_i/1.E6;
					dt3=ut3/1.E9-ut3_i/1.E9;
					timediff=dt1+dt2+dt3;
				}
				else if(Int_t(dt1)==1){
					dt2=(1-ut2_i/1.E6)+ut2/1.E6;
					dt3=(1-ut3_i/1.E9)+ut3/1.E9;
					timediff=dt2+dt3;
				}
				else timediff =100;
			}
			//establish a flag for the cuts
			if (nhits>nhitCut && nhits<nhitMax && rf<fiducialVolume /*&& thetaij>thetamin && thetaij<thetamax && itr>itrmin && itr<itrmax */ && damnflag==1) {
				cutFlag=1;
			}
			else {
				cutFlag=0;
			}
			//look for a vertex if cuts are passes and vertex flag has not been set
			if (cutFlag==1 /*&& (isFirst==0||timediff>timeWindow ||isThird==1)*/) {
				tempEventID=evid;
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
				nump++;
				
			}
			else {
				jdy_i=jdy;
			}

			//look for a 2 or 3 fold coincidence if vertex flag is set and cuts are passed
			/*
			else if (cutFlag==1 && isFirst==1 && timediff<timeWindow) {
				//check to see if first coincidence is flagged if not flag it
				if (isSecond==0) {
					eventID[numCoincidence][0]=tempEventID;
					eventID[numCoincidence][1]=evid;
					numCoincidence++;
					isSecond=1;
					numg1++;
					jdy_i=jdy;
					//cerr<<Int_t(tempEventID)<<" with energy "<<epos<<" at time "<<ut1_i<<" s+ "<<ut2_i<<" us"<<endl;
					//cerr<<Int_t(evid)<<" with energy "<<ene<<" at time "<<ut1<<" s+ "<<ut2<<" us"<<endl;

				}
				else if (isSecond==1) {
					eventID[numCoincidence-1][2]=evid;
					isThird=1;
					numg2++;
					//cerr<<Int_t(evid)<<" with energy "<<ene<<" at time "<<ut1<<" s+ "<<ut2<<" us"<<endl;
					//reset and start search over ***could be changed to look for bursts
					jdy_i=jdy;
					isFirst=0;
					isSecond=0;
					ut1_i=ut1;
					ut2_i=ut2;
					ut3_i=ut3;
					xi=xfit;
					yi=yfit;
					zi=zfit;
				}


			}
			//reset search for vertex if time window has closed and only cut events occured
			else if (cutFlag==0 && timediff>timeWindow) {
				isFirst=0;
				isSecond=0;
				isThird=0;
				xi=0;
				yi=0;
				zi=0;
			}
			
			if (evid==1406278) {
				cerr<<rdamn1<<" "<<rdamn2<<endl;
			}
			 */
			
			
		}//end loop over events

	}//end loop over ntuples
	cerr<<numCoincidence<<" candidate antinuetrino(s)"<<endl;
		clockedtime = clock[nd-1]-clock[0];
		rate= nump/clockedtime;
	for (Int_t n=0; n<numCoincidence; n++) {
		if(eventID[n][2]==0){
			cerr<<"2Fold event with event ids "<<eventID[n][0]<<" "<<eventID[n][1]<<" "<<endl;
		}
		if(eventID[n][2]!=0){
			cerr<<"3Fold event with event ids "<<eventID[n][0]<<" "<<eventID[n][1]<<" "<<eventID[n][2]<<endl;
		}
	}
		cerr<<nump<<" events in  "<<clockedtime<<" sec over "<<totaldays<<" days for rate of "<<rate<<endl;
		
}//end code

