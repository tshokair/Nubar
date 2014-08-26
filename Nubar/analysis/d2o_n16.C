/*
 *  d2o_in.C
 *  
 *
 *  Created by Tim Shokair on 11/8/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *  Code will be used for a mock analysis, rewrite of the look back search.
 *  Looks at the light water and acrylic events 
 *  12/9 fixed problem with duplicates two fold events (problem is with the sorted root file,
 *  but the analysis code was fixed to make sure no coincidences occured with different event numbers and the same gtid)
 *  1/31 Fixed problem with smaller window yeilding more events, it was a problem with duplicates. 
 *  2/1 Version verified against forward search. Ready for analysis.  
 *  4/19 updated energy to account for total vs. kinetic, also placing muon follower cut. 
 *	5/15 fixed bug that overcounts twofold
 */

{ 

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	//snoman variables and temp snoman variables
	Float_t nhits;
	Float_t ene;
	Float_t eneu;
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
	Float_t neutron;
	//counts;
	Int_t numCoin=0;
	Int_t numPred[500000]={0.0};
	Int_t coin_events[500000]={0.0};
	Int_t numTwoFold=0;
	Int_t numThreeFold=0;
	Int_t ct=0;
	Int_t doubleCounts=0;
	Int_t nd;
	Int_t numAcc;
	Int_t numCut;
	Int_t numPi=0;
	Int_t numOnePred=0;
	Int_t numTwoPred=0;
	
	//flags
	Int_t damnFlag1=0;
	Int_t damnFlag2=0;
	Int_t cleanFlag1=0;
	Int_t cleanFlag2=0;
	Int_t checkR=0;
	Int_t check2[500000]={0.0};
	Int_t twoFlag=0;
	Int_t threeFlag=0;
	Int_t duplicateFlag=0;
	Int_t duplicateFlag2=0;

	//times
	double totalTime;
	Int_t day_diff;
	double timediff;
	Double_t timeOfLast[500000][4];
	Double_t timeOfFirst[500000][4];
	//temps
	Int_t lastClean;
	Int_t tempID[500000];
	Int_t tempNhits;
	double tempTime[500000];
	Float_t tempEn[500000];
	double tempDt[500000];
	Int_t tID;
	Float_t tEn;
	double tTime;
	Int_t l=0;
	Float_t predX[100];
	Float_t predY[100];
	Float_t predZ[100];
	Int_t predID[100];
	Float_t predEn[100];
	double predDt[100];
	double evTime[100];
	Float_t deltaR[3];
	Float_t tempR;
	Float_t piEn;
	Float_t pe1En;
	Float_t pe2En;
	double piTime;
	double pe1Time;
	double pe2Time;
	Int_t piID;
	Int_t pe1ID;
	Int_t pe2ID;
	Int_t piNhit;
	Int_t pe1Nhit;
	Int_t pe2Nhit;
	Int_t lastCoinID1=0;
	Int_t lastCoinID2=0;

	//files
	Int_t numfiles;
	TString files[1000];
	
	//parameters
	Float_t nhitCut=20;
	Float_t nhitMax=160;
	Float_t fiducialVolume=600;
	Float_t timeWindow=0.150;
	Float_t distanceWindow=500;
	Float_t startWindow=100;
	Float_t betamin=-0.12;
	Float_t betamax=0.95;
	Float_t thetamin=0.75;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=4.5+.511;
	Float_t ethresh=3.9+.511;
	Float_t outerVolume=650;
	Float_t innerVolume=550;

	//readin the file names
	ifstream runList;
	//runList.open("fissD2O.txt");
	//runList.open("D2Ob8.txt");
	runList.open("D2Otest.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	TH1F * h1= new TH1F("h1","Time Differences D2O",500,0,.5);
	TH1F * h2= new TH1F("h2","Primary Energy",100,0,20);
	TH1F * h3= new TH1F("h3","Predicessor Energy",100,0,20);
	TH1F * h4= new TH1F("h4","Coin IDs",30000,0,30000);
	TH1F * h5= new TH1F("h5","Coin Ev nums",30000,0,30000);
	TH1F * h6= new TH1F("h6","Delta R",1300,0,1300);
	TH1F * h7= new TH1F("h7","Reconstructed radius",1300,0,1300);
	TH1F * h8= new TH1F("h8","Generated radius",1300,0,1300);
	TH1F * h9= new TH1F("h9","Nhits",250,0,250);
	TH1F * h10= new TH1F("h10","recon Energy",1000,0,10);
	TH1F * h11= new TH1F("h11","gen Energy",1000,0,10);
	//TProfile *px = h2->ProfileX("px", 0, 100);
	//TProfile *prof = a_2d_histogram->ProfileX();
	TH2F * h= new TH2F("h"," Primary Energy vs Predecessor Enegry",100,0,10,100,0,10);
	TH2F * p2= new TH2F("p2"," Nhit vs Reconstrcted Enegry",120,0,120,200,0,20);
	TH2F * p= new TH2F("p"," Predessor Time Separation vs Predecessor Space Separation",120,0,120,500,0,50);
	h->SetXTitle("Primary Energy [MeV]");
	h->SetYTitle("Predecessor Energy [MeV]");
	h9->SetXTitle("Nhits");
	h9->SetYTitle("Energy [MeV]");
	p->SetXTitle("Delta R [cm]");
	p->SetYTitle("Delta T [s]");
	p2->SetXTitle("Nhits");
	p2->SetYTitle("Reconstructed Energy [MeV]");
	h2->SetLineColor(2);
	h3->SetLineColor(4);
	h3->SetXTitle("Energy [MeV]");
	//graph points
	Double_t dt[100000];
	Double_t dr[100000];
	Int_t nGr=0;

	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k)//loop over ntuples/root files
    {	
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("h360");
		nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		//get the necessary information from ntuple
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Enekp",&ene);
		Td->SetBranchAddress("Enerau",&eneu);
		Td->SetBranchAddress("Egen",&egen);
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
		Td->SetBranchAddress("Xe",&xe);
		Td->SetBranchAddress("Ye",&ye);
		Td->SetBranchAddress("Ze",&ze);
		Td->SetBranchAddress("Thetaijp",&thetaij);
		Td->SetBranchAddress("Beta1p",&beta1);
		Td->SetBranchAddress("Beta4p",&beta4);
		Td->SetBranchAddress("Rdamn1",&rdamn1);
		Td->SetBranchAddress("Rdamn2",&rdamn2);
		lastClean=0;
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			if (rf>600.5 && rf<605.5 && ene>0) {
				p2->Fill(nhits,ene);
				h9->Fill(nhits,ene);
				//prof->Fill(nhits,ene);
				
				//cerr<<nhits<<" "<<ene<<endl;
			}
			
		}//end loop over events 
		delete Td;
		delete f1;

	}//end loop over ntuples
		TProfile *prof = p2->ProfileX();
		TF1 *fit = new TF1("fit","pol1",0,100);
		prof->Fit("fit");
		prof->Draw();
		//p2.Draw("colz");

}//end code