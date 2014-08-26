/*
 *  ncd_lb.C
 *  2/8 First pass at a lookback for the ncd phase.   
 *  2/13 Version 1.0 attemps to put in distance cut
 *  5/8 Version 1.1 makes sure coincidences are from the same generated particle for timing purposes.
 *  5/9 1.2 put in the analysis cuts and high level cuts
 *  6/25 1.3 adds the damn cuts
 */

{ 
	
	gROOT->SetStyle("Plain");
	//snoman variables and temp snoman variables
	Float_t nhits;
	Float_t ene;
	Float_t eneu;
	Float_t eneuo;
	Float_t egen;
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
	Float_t xeI;
	Float_t yeI;
	Float_t zeI;
	Float_t ri;
	Float_t x2;
	Float_t y2;
	Float_t xfit;
	Float_t yfit;
	Float_t zfit;
	Float_t rfit;
	Float_t rf;
	Float_t dt1;
	Float_t dt2;
	Float_t dt3;
	Float_t dr;
	Float_t dr2;
	Float_t rdamn1;
	Float_t rdamn2;
	Float_t rdamn3;
	Float_t rdamn4;
	Float_t evid;
	Float_t ncdstat;
	Float_t nes;
	Float_t nes_str;
	Float_t mux_str;
	//counts;
	Int_t numNCD=0;
	Int_t numPMT=0;
	Int_t numPi=0;
	Int_t numPred[5000000]={0.0};
	Int_t coins=0;
	Int_t numFN[5000000]={0.0};
	Int_t num3fold=0;
	Int_t num2fold=0;
	Int_t numNN=0;
	Int_t numNe=0;
	//Flags
	Int_t ncdFlag=0;
	Int_t ncdFlag2=0;
	Int_t pmtFlag=0;
	Int_t damnFlag=0;
	Int_t damnFlag0=0;
	Int_t damnFlag1=0;
	Int_t damnFlag2=0;
	Int_t damnFlag3=0;
	Int_t damnFlag4=0;
	//files
	Int_t numfiles;
	TString files[1000];
	//times
	double totalTime;
	Int_t day_diff;
	double startTime;
	double endTime;
	Float_t startDay;
	Float_t endDay;
	double timediff;
	double piTime;
	double peTime;
	double fnTime;
	Int_t predID[100]={0.0};
	Int_t predSt[100]={0.0};
	Float_t predEn[100]={0.0};
	Double_t predDt[100]={0.0};
	Float_t predDr[100]={0.0};
	Float_t predX[100]={0.0};
	Float_t predY[100]={0.0};
	Float_t predZ[100]={0.0};
	
	//temps
	Int_t tID;
	Float_t tEn;
	Int_t tempSt;
	Int_t l;
	Int_t lastN=0;
	Int_t fnIn=0;
	Int_t fnID=0;
	Int_t peID=0;
	Int_t lastpe=0;
	//parameters
	Float_t nhitCut=20;
	Float_t nhitMax=150;
	Float_t fiducialVolume=600;
	Float_t timeWindow=.300;
	Float_t distanceWindow=1275;
	Float_t dwNCD=1275;
	Float_t startWindow=100;
	Float_t betamin=0.2;
	Float_t betamax=0.8;
	Float_t thetamin=0.9;
	Float_t thetamax=1.4;
	Float_t itrmin= 0.65;
	Float_t itrmax=0.90;
	Float_t emin=4.0+.511;
	Float_t emax=50.0;
	Float_t nesmin=0.4;
	Float_t nesmax=0.9;
	Float_t outerVolume=650;
	Float_t innerVolume=550;
	Float_t stringX[40];
	Float_t stringY[40];
	
	stringX[0]=50.0;
	stringY[0]=50.0;
	stringX[1]=50.0;
	stringY[1]=150.0;
	stringX[2]=50.0;
	stringY[2]=250.0;
	stringX[3]=50.0;
	stringY[3]=350.0;
	stringX[4]=150.0;
	stringY[4]=250.0;
	stringX[5]=150.0;
	stringY[5]=150.0;
	stringX[6]=250.0;
	stringY[6]=150.0;
	stringX[7]=150.0;
	stringY[7]=50.0;
	stringX[8]=250.0;
	stringY[8]=50.0;
	stringX[9]=350.0;
	stringY[9]=50.0;
	stringX[10]=350.0;
	stringY[10]=-50.0;
	stringX[11]=250.0;
	stringY[11]=-50.0;
	stringX[12]=150.0;
	stringY[12]=-50.0;
	stringX[13]=250.0;
	stringY[13]=-150.0;
	stringX[14]=50.0;
	stringY[14]=-50.0;
	stringX[15]=150.0;
	stringY[15]=-150.0;
	stringX[16]=150.0;
	stringY[16]=-250.0;
	stringX[17]=50.0;
	stringY[17]=-150.0;	
	stringX[18]=50.0;
	stringY[18]=-250.0;
	stringX[19]=50.0;
	stringY[19]=-350.0;
	stringX[20]=-50.0;
	stringY[20]=-350.0;
	stringX[21]=-50.0;
	stringY[21]=-250.0;
	stringX[22]=-50.0;
	stringY[22]=-150.0;
	stringX[23]=-150.0;
	stringY[23]=-250.0;
	stringX[24]=-150.0;
	stringY[24]=-150.0;
	stringX[25]=-50.0;
	stringY[25]=-50.0;
	stringX[26]=-250.0;
	stringY[26]=-150.0;
	stringX[27]=-150.0;
	stringY[27]=-50.0;
	stringX[28]=-250.0;
	stringY[28]=-50.0;
	stringX[29]=-350.0;
	stringY[29]=-50.0;
	stringX[30]=-350.0;
	stringY[30]=50.0;
	stringX[31]=-250.0;
	stringY[31]=50.0;
	stringX[32]=-250.0;
	stringY[32]=150.0;
	stringX[33]=-150.0;
	stringY[33]=50.0;
	stringX[34]=-150.0;
	stringY[34]=150.0;
	stringX[35]=-50.0;
	stringY[35]=250.0;
	stringX[36]=-50.0;
	stringY[36]=350.0;
	stringX[37]=-50.0;
	stringY[37]=250.0;
	stringX[38]=-50.0;
	stringY[38]=150.0;
	stringX[39]=-50.0;
	stringY[39]=50.0;
	
	//readin the file names
	
	ifstream runList;
	
	runList.open("NCDreactorIn.txt");
	//runList.open("testNCD.txt");
	//runList.open("NCDfiss.txt");
	//runList.open("NCDalphaNd.txt");
	
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	ofstream cList;
	cList.open("cListNin.txt");
	ofstream tCoins;
	tCoins.open("tCoinsNin.txt");
	TH1F * h1= new TH1F("h1","NCD Position",1300,0,1300);
	TH1F * h2= new TH1F("h2","PMT position",1300,0,1300);
	TH1F * h3= new TH1F("h3","PMT energy",150,0,15);
	h3->SetXTitle("reconstructed energy [MeV]");
	TH1F * h4= new TH1F("h4","NCD energy",140,0,1.4);
	h4->SetXTitle("reconstructed energy [MeV]");
	TH1F * h5= new TH1F("h5","NCD-PMT Time Difference",500,0,.5);
	TH1F * h6= new TH1F("h6","NCD-NCD Time difference",500,0,.5);
	TH1F * h7= new TH1F("h7","Distance Separation neutron-positron",1300,0,1300);
	TH1F * h8= new TH1F("h8","Distance Separation neutron-neutron",1300,0,1300);
	TH1F * h9= new TH1F("h9","NCD Position",1300,0,1300);
	TH1F * h10= new TH1F("h10","Beta14",200,-0.5,1.5);
	TH1F * h11= new TH1F("h11","Theta_ij",200,0,2.0);
	TH1F * h12= new TH1F("h12","ITR",200,0,2.0);
	TH1F * h13= new TH1F("h13","Egen",1000,0,10.0);
	TH1F * h14= new TH1F("h14","Nhits",100,0,100);
	TH1F * h15= new TH1F("h15","Egen Positron",100,0,10.0);
	TH1F * h16= new TH1F("h16","PMT energy",150,0,15);
	h16->SetXTitle("reconstructed energy [MeV]");
	TH1F * h17= new TH1F("h17","NCD energy",140,0,1.4);
	TH1F * h18= new TH1F("h18","PMT nhits",100,0,100);
	
	
	h5->SetXTitle("Delta [s]");
	h6->SetXTitle("time [s]");
	h7->SetXTitle("distance [cm]");
	h8->SetXTitle("distance [cm]");
	h5->SetLineColor(2);
	h6->SetLineColor(3);
	
	h1->SetLineColor(7);
	h2->SetLineColor(6);
	h3->SetLineColor(8);
	h4->SetLineColor(6);	
	
	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k){//loop over ntuples/root files
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("t370");
		nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		//get the necessary information from ntuple
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Enerau",&eneu);
		Td->SetBranchAddress("Enerou",&eneuo);
		Td->SetBranchAddress("Egen",&egen);
		Td->SetBranchAddress("Itru",&itr);
		Td->SetBranchAddress("Ev_tid",&evid);	
		Td->SetBranchAddress("Time_jdy",&jdy); //julian day
		Td->SetBranchAddress("Time_s",&ut1); //time in seconds
		Td->SetBranchAddress("Time_us",&ut2); //time in nanoseconds
		Td->SetBranchAddress("Time_ns",&ut3); //time in nanoseconds
		Td->SetBranchAddress("Xfu",&xfit);
		Td->SetBranchAddress("Yfu",&yfit);
		Td->SetBranchAddress("Zfu",&zfit);
		Td->SetBranchAddress("Rfu",&rf);
		Td->SetBranchAddress("Xe",&xe);
		Td->SetBranchAddress("Ye",&ye);
		Td->SetBranchAddress("Ze",&ze);
		Td->SetBranchAddress("Thetaiju",&thetaij);
		Td->SetBranchAddress("Beta1u",&beta1);
		Td->SetBranchAddress("Beta4u",&beta4);
		Td->SetBranchAddress("Rdamn1",&rdamn1);
		Td->SetBranchAddress("Rdamn2",&rdamn2);
		Td->SetBranchAddress("Rdamn3",&rdamn3);
		Td->SetBranchAddress("Rdamn4",&rdamn4);
		Td->SetBranchAddress("Ncd_stat",&ncdstat);
		Td->SetBranchAddress("Nes_en",&nes);
		Td->SetBranchAddress("Nes_str",&nes_str);
		Td->SetBranchAddress("Nmux_str",&mux_str);
		
		//lastClean=0;
		Td->GetEntry(0);
		ut1_i=ut1;
		ut2_i=ut2;
		ut3_i=ut3;
		jdy_i=jdy;
		startDay=jdy;
		lastN=0;
		lastpe=0;
		startTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			//totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			beta14=beta1+4*beta4;
			//h14->Fill(nhits);
			h13->Fill(egen);
			ncdFlag=0;
			ncdFlag2=0;
			if (ncdstat==1 && nes >nesmin && nes < nesmax) {
				if((Int_t(rdamn1)& 0xCC44)==0 && (Int_t(rdamn2)& 0x0001)==0){
					damnFlag0 =1;  
				}
				else {
					damnFlag0=1;
				}
				if((Int_t(rdamn3)& 0x00FF)==0 && (Int_t(rdamn4)& 0xDDFE)==0){
					damnFlag1 =1;  
				}
				else {
					damnFlag1=1;
				} 
				
				
				if(damnFlag0==1 && damnFlag1==1){
					if((nes_str!=-999999 && nes_str!=3 && nes_str!=10 && nes_str!=20 && nes_str!=30 && nes_str!=0 && nes_str!=1 && nes_str!=8 && nes_str!=18 && nes_str!=26 && nes_str!=31)|| (mux_str!=-999999 && mux_str!=3 &&mux_str!=10 && mux_str!=20&& mux_str!=30&& mux_str!=0 && mux_str!=1 &&mux_str!=8 && mux_str!=18 && mux_str!=26 && mux_str!=31 )){
						ncdFlag=1;
						//cerr<<nes_str<<" "<<nes<<" "<<ut1<<endl;
						//cerr<<nes<<" "<<evid<<" "<<ut1<<endl;
					}
				}	
			}
			
			//cerr<<ncdstat<<" "<<evid<<endl;
			if (ncdFlag==1 ) {
				//h13->Fill(egen);
				h1->Fill(rf);
				numNCD++;
				//cerr<<nes_str<<" "<<xfit<<" "<<yfit<<" "<<zfit<<endl;
				tID=evid;
				h4->Fill(nes);
				ut1_i=ut1;
				ut2_i=ut2;
				ut3_i=ut3;
				jdy_i=jdy;
				ri=rf;
				xeI=xe;
				yeI=ye;
				zeI=ze;
				tEn=nes;
				
				if(nes_str!=-999999){
					xi=stringX[int(nes_str)];
					yi=stringY[int(nes_str)];
					tempSt=nes_str;
				}
				else if(mux_str!=-999999){
					xi=stringX[int(mux_str)];
					yi=stringY[int(mux_str)];
					tempSt=mux_str;
				}
				else {
					xi=0;
					yi=0;
					//cerr<<"no string pos"<<endl;
				}
				
				piTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
				numPred[numPi]=0;
				numFN[numPi]=0;
				l=lastpe;
				while (l<j) {
					Td->GetEntry(l);
					//get the time difference
					day_diff=jdy_i-jdy;
					beta14=beta1+4*beta4;
					if(Int_t(day_diff)==-1){
						dt1=(86400-ut1_i)+ut1;
						if(Int_t(dt1)==0){
							dt2=ut2/1.E6-ut2_i/1.E6;
							dt3=ut3/1.E9-ut3_i/1.E9;
							timediff=dt1+dt2 +dt3;
							
						}
						else if(Int_t(dt1)==1){
							dt2=(1-ut2/1.E6-ut3/1.E9)+ut2_i/1.E6+ut3_i/1.E9;
							dt3=(0);
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
							dt2=(1-ut2/1.E6-ut3/1.E9)+ut2_i/1.E6+ut3_i/1.E9;
							dt3=(0);
							timediff=dt2+dt3;
						}
						else if(Int_t(dt1)==-1){
							dt2=(1-ut2_i/1.E6-ut3_i/1.E9)+ut2/1.E6+ut3/1.E9;
							dt3=(0);
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
							dt2=(1-ut2/1.E6-ut3/1.E9)+ut2_i/1.E6+ut3_i/1.E9;
							dt3=(0);
							timediff=dt2+dt3;
						}
						else timediff=10;											
					}
					//if the time difference is within the window, check to see if it is clean
					if (ncdstat==0 && thetaij>thetamin && thetaij<thetamax && itr>itrmin && itr<itrmax && beta14>betamin && beta14<betamax) {
						if (eneu>1000|| eneu<0 ||eneu!=eneu) {
							ene=eneuo;
						}
						else{
							ene=eneu;
						}
						h15->Fill(egen);
						if (TMath::Abs(timediff)<timeWindow && ene>emin && ene<emax && rf<fiducialVolume) {
							h10->Fill(beta14);
							h11->Fill(thetaij);
							h12->Fill(itr);
							h5->Fill(timediff);
							
							if((Int_t(rdamn1)& 0x0CF7)==0 && (Int_t(rdamn2)& 0x6FE1)==0){
								damnFlag2 =1;  
								dr=TMath::Sqrt((xi-xfit)*(xi-xfit)+(yi-yfit)*(yi-yfit));
								
							}
							else {
								damnFlag2=1;
								dr=TMath::Sqrt((xi-xfit)*(xi-xfit)+(yi-yfit)*(yi-yfit));
								
							}
							if (damnFlag2==1 && dr<distanceWindow) {
								if (numPred[numPi]==0) {
									predSt[0]=tempSt;
									predID[numPred[numPi]]=tID;
									predEn[numPred[numPi]]=tEn;
									predDt[numPred[numPi]]=0;
									predDr[numPred[numPi]]=0;
									predX[numPred[numPi]]=xi;
									predY[numPred[numPi]]=yi;
									predZ[numPred[numPi]]=0;
								}
								h7->Fill(dr);
								predID[numPred[numPi]+1]=evid;
								predEn[numPred[numPi]+1]=ene;
								predDt[numPred[numPi]+1]=timediff;
								predDr[numPred[numPi]+1]=dr;
								predX[numPred[numPi]+1]=xfit;
								predY[numPred[numPi]+1]=yfit;
								predZ[numPred[numPi]+1]=zfit;
								//cerr<<"distance separation "<<dr<<" from ("<<xi<<","<<yi<<") to ("<<xfit<<","<<yfit<<")"<<endl;
								numPred[numPi]++;
								numPMT++;
								h3->Fill(ene);
								h2->Fill(rf);
								peID=evid;
								peTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
								lastpe=l;
							}
							
						}
					}
					l++;
				}//end while lookback loop
				//cerr<<numPred[numPi]<<endl;
				
				int lbEnd=0;
				if (j+5<nd) {
					lbEnd=j+5;
				}
				else {
					lbEnd=nd-1;
				}
				int m=j;
				
				while (m<=lbEnd) {
					//cerr<<"looking forward"<<endl;
					Td->GetEntry(m);
					if (ncdstat==1) {
						//cerr<<"forward neutron"<<endl;
					
						day_diff=jdy_i-jdy;
						beta14=beta1+4*beta4;
						if(Int_t(day_diff)==-1){
							dt1=(86400-ut1_i)+ut1;
							if(Int_t(dt1)==0){
								dt2=ut2/1.E6-ut2_i/1.E6;
								dt3=ut3/1.E9-ut3_i/1.E9;
								timediff=dt1+dt2 +dt3;
								
							}
							else if(Int_t(dt1)==1){
								dt2=(1-ut2/1.E6-ut3/1.E9)+ut2_i/1.E6+ut3_i/1.E9;
								dt3=(0);
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
								dt2=(1-ut2/1.E6-ut3/1.E9)+ut2_i/1.E6+ut3_i/1.E9;
								dt3=(0);
								timediff=dt2+dt3;
							}
							else if(Int_t(dt1)==-1){
								dt2=(1-ut2_i/1.E6-ut3_i/1.E9)+ut2/1.E6+ut3/1.E9;
								dt3=(0);
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
								dt2=(1-ut2/1.E6-ut3/1.E9)+ut2_i/1.E6+ut3_i/1.E9;
								dt3=(0);
								timediff=dt2+dt3;
							}
							else timediff=10;											
						}
						ncdFlag2=0;
						if (ncdstat==1 && nes >nesmin && nes < nesmax) {
							if((Int_t(rdamn1)& 0xCC44)==0 && (Int_t(rdamn2)& 0x0001)==0){
								damnFlag3 =1;  
							}
							else {
								damnFlag3=1;
							}
							if((Int_t(rdamn3)& 0x00FF)==0 && (Int_t(rdamn4)& 0xDDFE)==0){
								damnFlag4 =1;  
							}
							else {
								damnFlag4=1;
							} 
							if(damnFlag3==1 && damnFlag4==1){
								if((nes_str!=-999999 && nes_str!=3 && nes_str!=10 && nes_str!=20 && nes_str!=30 && nes_str!=0 && nes_str!=1 && nes_str!=8 && nes_str!=18 && nes_str!=26 && nes_str!=31)|| (mux_str!=-999999 && mux_str!=3 &&mux_str!=10 && mux_str!=20&& mux_str!=30&& mux_str!=0 && mux_str!=1 &&mux_str!=8 && mux_str!=18 && mux_str!=26 && mux_str!=31 )){
									ncdFlag2=1;
									//cerr<<nes_str<<" "<<nes<<" "<<ut1<<endl;
									//cerr<<nes<<" "<<evid<<" "<<timediff<<endl;
								}
							}
						}
						
						if(TMath::Abs(timediff)<timeWindow && ncdFlag2==1){
							//cerr<<nes_str<<endl;
							int index;
							if (nes_str!=-999999) {
								index=nes_str;
							}
							else if (mux_str!=-999999){
								index=mux_str;
							}
							else {
								index=999;
							}
							
							
							if(index!=999){
								x2=stringX[index];
								y2=stringY[index];
								dr2=TMath::Sqrt((x2-xi)*(x2-xi)+(y2-yi)*(y2-yi));
								h8->Fill(dr2);
							}
							
							else {
								dr2=99999;
							}

							if (dr2<dwNCD) {
								//cerr<<nes<<endl;
								if (numPred[numPi]==0) {
									predSt[0]=tempSt;
									predID[numPred[numPi]]=tID;
									predEn[numPred[numPi]]=tEn;
									predDt[numPred[numPi]]=0;
									predDr[numPred[numPi]]=0;
									predX[numPred[numPi]]=xi;
									predY[numPred[numPi]]=yi;
									predZ[numPred[numPi]]=0;
								}
								predSt[numFN[numPi]+1]=nes_str;
								predID[numPred[numPi]+numFN[numPi]+1]=evid;
								predEn[numPred[numPi]+numFN[numPi]+1]=nes;
								predDt[numPred[numPi]+numFN[numPi]+1]=timediff;
								predDr[numPred[numPi]+numFN[numPi]+1]=dr2;
								predX[numPred[numPi]+numFN[numPi]+1]=x2;
								predY[numPred[numPi]+numFN[numPi]+1]=y2;
								predZ[numPred[numPi]+numFN[numPi]+1]=0;
								numFN[numPi]++;
								fnIn=m;
								fnID=evid;
								fnTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
								h6->Fill(timediff);
							}
						}
					}
					m++;
				}
				
				if (numPred[numPi] + numFN[numPi]>2) {
					cerr<<numPred[numPi]<<"(pmt)+"<<numFN[numPi]+1<<"(ncd)-fold event"<<endl;
					cList<<tID<<" :"<<numPred[numPi]<<"(pmt)+"<<numFN[numPi]+1<<"(ncd)-fold event"<<endl;
					if (numFN[numPi]==1) {
						j=fnIn;
					}
					lastN=j;
					
				}
				
				if (numPred[numPi]==1 || numFN[numPi]==1) {
					if (numPred[numPi]==1 && numFN[numPi]==1) {
						//cerr<<"3-fold: "<<peID<<" "<<tID<<" "<<fnID<<endl;
						//cerr<<"At times: "<<peTime<<" "<<piTime<<" "<<fnTime<<endl;
						cList<<"3-fold"<<endl;
						cList<<predID[0]<<" "<<predEn[0]<<" string "<<predSt[0] <<" ("<<predX[0]<<","<<predY[0]<<")"<<endl;
						cList<<predID[1]<<" "<<predEn[1]<<" ("<<predX[1] <<","<<predY[1]<<","<<predZ[1]<<") "<<predDr[1]<<" "<<predDt[1]<<endl;
						cList<<predID[2]<<" "<<predEn[2]<<" string "<<predSt[1] <<" ("<<predX[2]<<","<<predY[2]<<") "<<predDr[2]<<" "<<predDt[2]<<endl;
						num3fold++; 
					}
					else if (numPred[numPi]==1 && numFN[numPi]==0) {
						numNe++;
						cList<<"2-fold"<<endl;
						cList<<predID[0]<<" "<<predEn[0]<<" string "<<predSt[0] <<" ("<<predX[0]<<","<<predY[0]<<")"<<endl;
						cList<<predID[1]<<" "<<predEn[1]<<" ("<<predX[1] <<","<<predY[1]<<","<<predZ[1]<<") "<<predDr[1]<<" "<<predDt[1]<<endl;
						
					}
					else if (numPred[numPi]==0 && numFN[numPi]==1) {
						numNN++;
						cList<<"2-fold"<<endl;
						cList<<predID[0]<<" "<<predEn[0]<<" string "<<predSt[0] <<" ("<<predX[0]<<","<<predY[0]<<")"<<endl;
						cList<<predID[1]<<" "<<predEn[1]<<" string "<<predSt[1] <<" ("<<predX[1]<<","<<predY[1]<<") "<<predDr[1]<<" "<<predDt[1]<<endl;
						
					}
					
					if (numFN[numPi]==1) {
						j=fnIn;
					}
					lastN=j;
				}
				if (numPred[numPi] + numFN[numPi]>0) {
					coins++;
					numPi++;
				}
				//cerr<<nes<<" "<<evid<<" "<<ut1<<endl;
				
			}
			endDay=jdy;
			endTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			
		}//end loop over events 
		delete Td;
		delete f1;
		cerr<<"next ntuple"<<endl;
		//end loop over ntuples
		}
		cerr<<numNCD<<" ncd events and "<<numPMT<<" pmt events"<<endl;
		cerr<<coins<<" coincidences"<<endl;
		cerr<<num3fold<<" 3-folds"<<endl;
		cerr<<numNN<<" 2-fold (n,n)"<<endl;
		cerr<<numNe<<" 2-fold (n,e+)"<<endl;
		tCoins<<numNCD<<" ncd events and "<<numPMT<<" pmt events"<<endl;
		tCoins<<coins<<" coincidences"<<endl;
		tCoins<<num3fold<<" 3-folds"<<endl;
		tCoins<<numNN<<" 2-fold (n,n)"<<endl;
		tCoins<<numNe<<" 2-fold (n,e+)"<<endl;
		tCoins<<"Run started on day "<<startDay<<" at time "<<startTime<<" and ended on day "<<endDay<<" at time "<<endTime<<endl;
		h4->DrawNormalized();
		h3->DrawNormalized("same");
		cList.close();
		tCoins.close();
		h5->SetXTitle("Delta t [s]");
		h5->SetYTitle("N(Delta t) [events/.001s]");	
		h5->Draw("");
		h5->Fit("expo","","",0,.25);
		
		h6->SetXTitle("Delta t [s]");
		h6->SetYTitle("N(Delta t) [events/.001s]");	
		h6->Draw();
		h6->Fit("expo","","",0,.18);
		//h5->Fit("expo","","",0,.5);
		//h6->Fit("expo","","",0,.5);
		
	}//end code