/*
 *  coin.C
 *  
 *
 *  Created by Tim Shokair on 3/29/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *
 */
Int_t coin_backwards()
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
	Double_t timeSinceLast;
	Double_t vertextime;
	Double_t vertext1;
	Double_t vertext2;
	Double_t vertext3;
	Double_t vertexday;
	Int_t evNumOfLast;
	
	//parameters used in cuts
	Float_t nhitCut=27.5;
	Float_t nhitMax=160;
	Float_t fiducialVolume=550;
	Float_t timeWindow=0.5;
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
	TH1F * h2= new TH1F("h2","Nhits",100,0,100000);

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
		timediff=ut1+ut2/1.E6+ut3/1.E9;
		evNumOfLast=0;
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
			
			//find the time difference from last event if:
			//daydiff =1

			//establish a flag for the cuts
			if (nhits>nhitCut && nhits<nhitMax && rf<fiducialVolume && thetaij>thetamin && thetaij<thetamax && itr>itrmin && itr<itrmax && damnflag==1 &&ene>emin) {
				cutFlag=1;
				//cerr<<evid<<endl;
				
				vertextime=totalTime;
				ut1_i=ut1;
				ut2_i=ut2;
				ut3_i=ut3;
				jdy_i=jdy;
				
			
			}
			else {
				cutFlag=0;
			}
			//if an event is stored, look back to see if other events occured in the time window
			if (cutFlag==1) {
				//cerr<<evid<<endl;
				//Td->GetEntry(evNumOfLast);
				//cerr<<ut1<<" "<<ut1_i<<" "<<timeDifference(ut1,ut1_i,ut2,ut2_i,ut3,ut3_i,jdy_i,jdy)<<endl;
				for (Int_t l=evNumOfLast; l<=j; ++l) {
					Td->GetEntry(evNumOfLast);
					timediff=timeDifference(ut1,ut1_i,ut2,ut2_i,ut3,ut3_i,jdy,jdy_i);
					cerr<<evid<<" "<<timediff<<endl;
					if (timediff<timeWindow) {
						h1->Fill(nhits);
						h2->Fill(ut1);
						cerr<<timediff<<" "<<ut1<<" "<<ut1_i<<endl;
						numev++;
					}
				}
			
				numCoincidence++;
				evNumOfLast=j;
			}
		
			
		}//end loop over events

	}//end loop over ntuples
	cerr<<numCoincidence<<" candidate antinuetrino(s)"<<endl;
	cerr<<numev<<" total particles found"<<endl;

	h1->Draw();
	return 0;	
}//end code

		
Float_t timeDifference(Float_t sec_i,Float_t sec,Float_t usi,Float_t us, Float_t nsi, Float_t ns, Float_t day_i, Float_t day)
{
	Float_t deltat;
	Float_t dts;
	Float_t dtus;
	Float_t dtns;
	if(Int_t(day-day_i)==1){
		dts=(86400-sec_i)+sec;
		if(Int_t(dts)==0){
			dtus=us/1.E6-usi/1.E6;
			dtns=ns/1.E9-nsi/1.E9;
			deltat=dts+dtus +dtns;
			
		}
		else if(Int_t(dts)==1){
			dtus=(1-usi/1.E6)+us/1.E6;
			dtns=(1-nsi/1.E9)+ns/1.E9;
			deltat=dtus+dtns;
		}
		else deltat=1;
		
		
	}
	//daydiff=0
	else if(Int_t(day-day_i)==0){
		dts=sec-sec_i;
		if(Int_t(dts)==0){
			dtus=us/1.E6-usi/1.E6;
			dtns=ns/1.E9-nsi/1.E9;
			deltat=dts+dtus+dtns;
		}
		else if(Int_t(dts)==1){
			dtus=(1-usi/1.E6)+us/1.E6;
			dtns=(1-nsi/1.E9)+ns/1.E9;
			deltat=dtus+dtns;
		}
		else deltat =1;
	}
	return deltat;
	
}
