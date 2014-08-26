/*
 8/16/2011 avn.C
 Code to count the neutrons that are able to travel into the AV from the light water
 
 
 
 */
{
	gROOT->SetStyle("Plain");

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
	Float_t evid;
	Int_t partNum;
	Int_t tempEventID;
	Int_t isFirst;
	Int_t isSecond[100]={0.0};
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
	Int_t tid;
	Int_t trf;
	Float_t neutron;
	Float_t xg;
	Float_t yg;
	Float_t zg;
	Float_t rg;
	Float_t rg_i;
	Float_t pnhits;
	//counted variables
	Int_t coins[100]={0.0};
	Int_t penNeutrons[100]={0.0};
	Int_t numIn=0;
	Int_t totalEvents=0;
	
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
	Float_t nhitCut=5;
	Float_t nhitMax=160;
	Float_t fiducialVolume=600;
	Float_t timeWindow=1.000;
	Float_t startWindow=100;
	Float_t thetamin=0.63;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=5.1;
	Float_t ethresh=4.3;
	Float_t outerVolume=670;
	Float_t fV=0;
	Int_t rstep=10;
	Int_t numstep=7;
	//readin the file names
	ifstream runList;
	runList.open("mcruns.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	
	//initialize the counted variables
	//numev[=0;
	nump=0;

	numCoincidence=0;
	numRawEvents=0;
	liveTime=0;
	runTime=0;
	Int_t real3fold=0;
	Int_t real2fold=0;
	//make histograms
	TH1F * h1= new TH1F("h1","Particle Pos as function of event time",50000,75000,80000);
	TH1F * h2= new TH1F("h2","Time Difference",10000,0,2000);
	TH1F * h3= new TH1F("h3","R Dist",800,0,800);
	TH1F * h4= new TH1F("h4","R vs. EVID",500,0,500);
	TH1F * h5= new TH1F("h5","Particle Pos as function of event time",50000,75000,80000);
	TH1F * h6= new TH1F("h6","event time [ut1]",50000,75000,80000);
	TH1F * h7= new TH1F("h7","Neutron energy",500,0,20);
	TH1F * h8= new TH1F("h8","Pos energy",500,0,20);
	TH1F * h9= new TH1F("h9","Neutron NHit in AV",500,0,100);	
	TH1F * h10= new TH1F("h10","Neutron NHit",500,0,100);	
	TH1F * h11= new TH1F("h11","Positron NHit",500,0,100);	
	h7->SetLineColor(8);
	h8->SetLineColor(7);
	h9->SetLineColor(5);
	h10->SetLineColor(6);
	h11->SetLineColor(4);
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
		totalEvents+=nd;
		cerr<<nd<<" events in file "<<k<<" and "<<totalEvents<<" events total"<<endl;
		
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
		Td->SetBranchAddress("Xe",&xg);
		Td->SetBranchAddress("Ye",&yg);
		Td->SetBranchAddress("Ze",&zg);
		//Td->SetBranchAddress("Neutron",&neutron);
		Td->GetEntry(0);
		
		//initialize variables for first event
		isFirst=0;
		//isSecond=0;
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
			//tid=evid;
			//get the vector between this event and the last one
			rfit=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi));
			rg=TMath::Sqrt(xg*xg+yg*yg+zg*zg);
			//get the day differce to calculate time difference correctly
			day_diff=Int_t(jdy)-Int_t(jdy_i);
			totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			if(Int_t(day_diff)==-1){
				dt1=(86400-ut1)+ut1_i;
				if(Int_t(dt1)==0){
					dt2=ut2/1.E6-ut2_i/1.E6;
					dt3=ut3/1.E9-ut3_i/1.E9;
					timediff=dt1+dt2 +dt3;
					
				}
				else if(Int_t(dt1)==-1){
					dt2=(1-ut2/1.E6)+ut2_i/1.E6;
					dt3=(1-ut3/1.E9)+ut3_i/1.E9;
					timediff=dt2+dt3;
				}
				else timediff=dt1;
				
			}
			//daydiff=0
			else if(Int_t(day_diff)==0){
				dt1=ut1-ut1_i;
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
				else {
					timediff =dt1;
					//cerr<<"dt1: "<<dt1<<" t ="<<totalTime<<endl;
				}
				
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
				else {
					timediff=dt1;
					//cerr<<"dt1: "<<dt1<<endl;
				}
			}
			else {
				timediff=-999;
				cerr<<"day diff:" <<day_diff<<endl;
			}
			if(timediff!=-999){
				h2->Fill(timediff);
				//cerr<<timediff<<endl;
			}

			
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
			if (nhits>nhitCut && nhits<nhitMax && rf>0&& damnflag==1 &&ene>ethresh) {
				cutFlag=1;
				h3->Fill(rf);
				h4->Fill(evid,rf);
				h5->Fill(totalTime,rf);
				h6->Fill(ut1);
				if(rf<620){
					numIn++;
				}
			}
			else {
				cutFlag=0;
				//cerr<<evid<<" did not pass the cuts, energy: "<<ene<<" nhits " <<nhits<<endl;

			}
			
			if( timediff<timeWindow){
				isSecond[0]=1;
				isFirst=0;
				coins[0]++;
			}
			else {
				isSecond[0]=0;
				isFirst=1;
			}


			
			nump++;
			trf=rf;
			tempTime=totalTime;
			tid=evid;
			jdy_i=jdy;
			ut1_i=ut1;
			ut2_i=ut2;
			ut3_i=ut3;
			rg_i=rg;
			epos=ene;
			pnhits=nhits;
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
		cerr<<numIn<<" events inside 620 out of "<<totalEvents<<" events total"<<endl;
		for (Int_t j=0; j<numstep; j++) {
			cerr<<coins[j]<<" coincidences found for an outer volume of "<<outerVolume-rstep*j<<endl;
			cerr<<penNeutrons[j]<<" of those coincidences had neutrons that made it into the AV"<<endl;
		}
		h10->Draw();
		h9->Draw("same");
		h11->Draw("same");
		//gROOT->SetStyle("Plain");
		
}//end code
		
