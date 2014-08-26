/*
 8/18/2011 rdist.C
 Get the distance distribution of random events 
 
 
 
 */
{
	//declare variables
	Float_t nhits;
	Float_t ene;
	Float_t itr;
	Float_t evid;
	Float_t xfit;
	Float_t yfit;
	Float_t zfit;
	Float_t rf;
	Float_t thetaij;
	Float_t rdamn1;
	Float_t rdamn2;
	Float_t xi;
	Float_t yi;
	Float_t zi;
	TRandom2 rn;
	Int_t rIndex;
	Double_t deltaR;
	Int_t numfiles;
	TString files[10000];
	Float_t norm;
	//flags
	Int_t damnflag=0;
	Int_t cutFlag=0;
	Int_t damnflag2=0;
	Int_t cutFlag2=0;

	
	//parameters used in cuts
	Float_t nhitCut=5;
	Float_t nhitMax=160;
	Float_t fiducialVolume=600;
	Float_t distanceWindow=200;
	Float_t thetamin=0.63;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	Float_t emin=3.5;
	Float_t ethresh=2.7;
	Float_t outerVolume=650;
	Float_t fV=0;
	Int_t rstep=20;
	Int_t numstep=7;
	
	//counter
	Int_t totatEV=0;
	Float_t totalClean=0;
	Float_t numInWindow=0;
	
	//readin the file names
	ifstream runList;
	runList.open("test.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	
	//initialize the counted variables
	
	//make histograms
	TH1F * h1= new TH1F("h1","Delta R Distribution",650,0,(outerVolume*2+50));
	

	//loop over ntuples
	for (Int_t k=0; k<(numfiles-1);++k)//loop over ntuples/root files
    {	
		//read in events in ntuple
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("t360");
		nd = Td->GetEntries();
		norm=nd;
		cerr<<nd<<" events in file "<<k<<endl;
		
		//get the necessary information from ntuple
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Enekp",&ene);
		Td->SetBranchAddress("Itrp",&itr);
		Td->SetBranchAddress("Ev_tid",&evid);	

		Td->SetBranchAddress("Xfp",&xfit);
		Td->SetBranchAddress("Yfp",&yfit);
		Td->SetBranchAddress("Zfp",&zfit);
		Td->SetBranchAddress("Rfp",&rf);
		Td->SetBranchAddress("Thetaijp",&thetaij);
		Td->SetBranchAddress("Rdamn1",&rdamn1);
		Td->SetBranchAddress("Rdamn2",&rdamn2);
		rIndex = Int_t(rn.Rndm()*norm);
		Td->GetEntry(rIndex);
		xi=xfit;
		yi=yfit;
		zi=zfit;
		cerr<<"starting with index "<<rIndex<<" and random number "<<rn.Rndm()<<endl;

		
		
		//begin loop over event
		for(Int_t j=0; j<nd; ++j){
			totatEV++;
			rIndex = Int_t(rn.Rndm()*norm);
			cerr<<"getting info for index "<<rIndex<<endl;
			Td->GetEntry(rIndex);
			xi=xfit;
			yi=yfit;
			zi=zfit;
			if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
				damnflag =1;  
			}
			else {
				damnflag=0;
			}
			rIndex = Int_t(rn.Rndm()*norm);
			//cerr<<"getting info for index "<<rIndex<<endl;
			Td->GetEntry(rIndex);
			//cerr<<"event "<<evid<<" started at "<<ut1<<"+"<<ut2<<endl;
			//set the damn flag for this event (UNIDOC)
			

			if (nhits>nhitCut && nhits<nhitMax && rf<outerVolume&& damnflag==1 &&ene>ethresh) {
				cutFlag=1;
			}
			
			else {
				cutFlag=0;
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
			
			
			//***Look for a second random particle
			rIndex = Int_t(rn.Rndm()*norm);
			//cerr<<"getting info for index "<<rIndex<<endl;
			Td->GetEntry(rIndex);
			
			if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
				damnflag2 =1;  
			}
			else {
				damnflag2=0;
			}
			
			//establish a flag for the cuts
			if (nhits>nhitCut && nhits<nhitMax && rf<outerVolume&& damnflag2==1 &&ene>ethresh) {
				cutFlag2=1;
			}
			
			else {
				cutFlag2=0;
			}
			if (cutFlag2==1 &&cutFlag==1) {
				deltaR=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi));
				h1->Fill(deltaR);
				totalClean++;
				if (deltaR<distanceWindow) {
					numInWindow++;
				}
			}
			
		}//end loop over events
		delete Td;
		delete f1;

		}//end loop over ntuples
		h1->Draw();
		cerr<<totatEV<<" total events: "<<numInWindow/totalClean<<"% were separated by less than "<<distanceWindow<<" cm" <<endl;
		
}//end code
		
