/*
 *  ncd_test.C
 *  This code will test ncd triggering and measure event rates in ncds and pmts.  
 *
 *  
 */

{ 
	
	gROOT->SetStyle("Plain");
	//snoman variables and temp snoman variables
	Float_t nhits;
	Float_t ene;
	Float_t eneu;
	Float_t eneo;
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
	Int_t numPMTout=0;
	Int_t steps=25;
	Float_t stepSize=0.1;
	Int_t numIn[steps]={0.0};
	Int_t numOut[steps]={0.0};


	//Flags
	Int_t ncdFlag=0;
	Int_t pmtFlag=0;
	Int_t damnFlag=0;
	Int_t damnFlag0=0;
	Int_t damnFlag1=0;
	//files
	Int_t numfiles;
	TString files[10000];
	
	//parameters
	Float_t nhitCut=20;
	Float_t nhitMax=160;
	Float_t fiducialVolume=600;
	Float_t timeWindow=0.100;
	Float_t distanceWindow=350;
	Float_t startWindow=100;
	Float_t betamin=-0.12;
	Float_t betamax=0.95;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=2.7+.511;
	Float_t ethresh=3.5 +.511;
	Float_t emax=50.0;
	Float_t nesmin=0.65;
	Float_t nesmax=0.85;
	Float_t outerVolume=650;
	Float_t innerVolume=550;
	
	//readin the file names
	
	ifstream runList;
	
	runList.open("runlistNCD.runlist");
	//cerr<<"files open"<<endl;
	numfiles =0;
	while(!runList.eof()){
	  //cerr<<numfiles<<endl;
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	//cerr<<"files closed"<<endl;

	TH1F * h1= new TH1F("h1","NCD Position",1300,0,1300);
	TH1F * h2= new TH1F("h2","PMT position",1300,0,1300);
	TH1F * h3= new TH1F("h3","PMT energy",150,0,15);
	TH1F * h4= new TH1F("h4","NCD energy",150,0,15);
	TH1F * h5= new TH1F("h5","PMT energy (outside)",150,0,15);
	TH1F * h6= new TH1F("h6","PMT position (outside",1300,0,1300);
	h1->SetLineColor(7);
	h2->SetLineColor(6);
	h3->SetLineColor(4);
	h4->SetLineColor(5);	
	double livetime=385.17*24*3600;
	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k){//loop over ntuples/root files
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("h370");
		Int_t nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		//get the necessary information from ntuple
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Enerau",&eneu);
		Td->SetBranchAddress("Enerou",&eneo);
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
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			//totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			beta14=beta1+4*beta4;

			if (ncdstat==1 && nes >nesmin && nes < nesmax) {
			  if((Int_t(rdamn1)& 0xCC44)==0 && (Int_t(rdamn2)& 0x0001)==0){
			    damnFlag0 =1;  
			  }
			  else {
			    damnFlag0=0;
			  }
			  if((Int_t(rdamn3)& 0x00FF)==0 && (Int_t(rdamn4)& 0xDDFE)==0){
			    damnFlag1 =1;  
			  }
			  else {
			    damnFlag1=0;
			  } 
			  
			  
			  if(damnFlag0==1 && damnFlag1==1){
			    if((nes_str!=-999999 && nes_str!=3 && nes_str!=10 && nes_str!=20 && nes_str!=30 && nes_str!=0 && nes_str!=1 && nes_str!=8 && nes_str!=18 && nes_str!=26 && nes_str!=31)|| (mux_str!=-999999 && mux_str!=3 &&mux_str!=10 && mux_str!=20&& mux_str!=30&& mux_str!=0 && mux_str!=1 &&mux_str!=8 && mux_str!=18 && mux_str!=26 && mux_str!=31 )){
			      numNCD++;
			      h4->Fill(nes);
			      //cerr<<nes_str<<" "<<nes<<" "<<ut1<<endl;
			      //cerr<<nes<<" "<<evid<<" "<<ut1<<endl;
			    }
			  }	
			}

			else if (ncdstat==0){
			  if((Int_t(rdamn1)& 0x0CF7)==0 && (Int_t(rdamn2)& 0x6FE1)==0){
			    damnFlag =1;  
			  }
			  else {
			    damnFlag=0;
			  }
			  if(eneu>1000 || ene<0 || eneu!=eneu){
			    ene=eneo;
			  }
			  else{ ene=eneu;}
			  
			  if ( rf<600 &&damnFlag==1 && ene>emin && ene<emax){
				h2->Fill(rf);
				h3->Fill(ene);
				//cerr<<"pmt "<<ene<<" "<<ut1<<endl;
				numPMT++;
				for (int in=0; in<steps; in++) {
					if(ene>emin+in*stepSize){
						numIn[in]++;
					}
				}
			  }
			
			  if (rf<outerVolume && rf>innerVolume && ene>emin &&damnFlag==1 && ene<emax) {
				h6->Fill(rf);
				h5->Fill(ene);
				//cerr<<"pmt "<<ene<<" "<<ut1<<endl;
				numPMTout++;
				for (int ot=0; ot<steps; ot++) {
				  if(ene>(emin+ot*stepSize)){
				    numOut[ot]++;
				  }
				}
			  }
			}

		}//end loop over events 
		delete Td;
		delete f1;
		cerr<<"next ntuple"<<endl;
		//end loop over ntuples
		}
		//h4->DrawNormalized();
		//h3->DrawNormalized("same");
		cerr<<numNCD<<" ncd events "<<endl;


		for(int n=0;n<steps;n++){
		  cerr<<"At "<<emin+n*stepSize-.511<<" MeV there were " <<numIn[n]<<" positron inside and "<<numOut[n]<<" positrons outside"<<endl;
		}
		for(int n=0;n<steps;n++){
			cerr<<"At "<<emin+n*stepSize-.511<<" MeV there were " <<numIn[n]/livetime<<" Hz of positron inside and "<<numOut[n]/livetime<<" Hz positrons outside"<<endl;
		}
		
	}//end code
