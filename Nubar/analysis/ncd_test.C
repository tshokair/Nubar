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
	Float_t evid;
	Float_t ncdstat;
	Float_t nes;
	Float_t nes_str;
	Float_t mux_str;
	//counts;
	Int_t numNCD=0;
	Int_t numPMT=0;
	Int_t numPMTout=0;

	//Flags
	Int_t ncdFlag=0;
	Int_t pmtFlag=0;
	//files
	Int_t numfiles;
	TString files[1000];
	
	//parameters
	Float_t nhitCut=20;
	Float_t nhitMax=160;
	Float_t fiducialVolume=600;
	Float_t timeWindow=0.350;
	Float_t distanceWindow=500;
	Float_t startWindow=100;
	Float_t betamin=-0.12;
	Float_t betamax=0.95;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=3.5;
	Float_t ethresh=3.5;
	Float_t outerVolume=650;
	Float_t innerVolume=550;
	
	//readin the file names
	
	ifstream runList;
	
	runList.open("testNCD.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	TH1F * h1= new TH1F("h1","NCD Position",1300,0,1300);
	TH1F * h2= new TH1F("h2","PMT position",1300,0,1300);
	TH1F * h3= new TH1F("h3","PMT energy",150,0,15);
	TH1F * h4= new TH1F("h4","NCD energy",150,0,15);

	h1->SetLineColor(7);
	h2->SetLineColor(6);
	h3->SetLineColor(4);
	h4->SetLineColor(5);	
	
	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k){//loop over ntuples/root files
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("t370");
		nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		//get the necessary information from ntuple
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Enerau",&ene);
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
		/*
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			//totalTime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
			beta14=beta1+4*beta4;
			//cerr<<ncdstat<<" "<<evid<<endl;
			if (ncdstat==1) {
				h1->Fill(rf);
				numNCD++;
				h4->Fill(nes);
				//cerr<<nes_str<<" "<<nes<<" "<<ut1<<endl;
				//cerr<<nes<<" "<<evid<<" "<<ut1<<endl;
				
			}
			else if (ncdstat==0 && rf<600 && ene>2.7) {
				h2->Fill(rf);
				h3->Fill(ene);
				//cerr<<"pmt "<<ene<<" "<<ut1<<endl;
				numPMT++;
			}
			if (ncdstat==0 && rf<650 && rf>550 && ene>2.7) {
				h2->Fill(rf);
				h3->Fill(ene);
				//cerr<<"pmt "<<ene<<" "<<ut1<<endl;
				numPMTout++;
			}

		}*/
		//end loop over events 
		 
		delete Td;
		delete f1;
		cerr<<"next ntuple"<<endl;
		//end loop over ntuples
		}
		cerr<<numNCD<<" ncd events and "<<numPMT<<" pmt events in and  "<<numPMTout<<" out"<<endl;
		//h4->DrawNormalized();
		//h3->DrawNormalized("same");
		
	}//end code