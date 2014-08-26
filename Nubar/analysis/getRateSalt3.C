//Code to calculate the accidental coincidence rates from the Salt Phase of SNO

//begin
{
  //declare the SNOMAN variables
  Float_t nhits;
  Float_t ene;
  Float_t evid;
  Float_t jdy;
  Float_t ut1;
  Float_t ut2;
  Float_t ut3;
  Float_t xfit;
  Float_t yfit;
  Float_t zfit;
  Float_t rf;
  Float_t itr;
  Float_t thetaij;
  Float_t beta1;
  Float_t beta4;
  Float_t beta14;
  Float_t rdamn1;
  Float_t rdamn2;
  Int_t nd;

  //declare the counted variables
  Int_t numClean[1000];
  Int_t numpi[1000][50];
  Float_t estep;
  Int_t numcl;
  Int_t numtag[50][50]={0.0};
  Double_t livetime;

  //declare the initializing variables
  Float_t xi;
  Float_t yi;
  Float_t zi;
  Float_t ri;
  Float_t jdyi;
  Float_t ut1i;
  Float_t ut2i;
  Float_t ut3i;

  //declare the arrays associated with rates
  Float_t cleanRate;
  Float_t rate[1000][50];
  Float_t runtime;
  Int_t dayi;
  Int_t day;
  Double_t totaltimei;
  Double_t totaltime;
  Int_t totaldays;
  Float_t tagrate[50];

  //declare the flags used for cuts
  Int_t cleanFlag;
  Int_t damnFlag;
  Int_t piFlag;

  //array for file names
  TString files[1000];
  Int_t numfiles;
  Float_t scale;
  //const char *hist[10];
  
  //declare and initialize cut parameters;
  Float_t nhitmin=20;
  Float_t nhitmax=160;
  Float_t piFV=555;
  Float_t cleanFV=620;
  Float_t timeWindow=.050;
  Float_t distanceWindow =300;
  Float_t betamin=-.12;
  Float_t betamax=.95;
  Float_t thetamin=.63;
  Float_t thetamax=1.45;
  Float_t itrmin=.55;
  Float_t itrmax=.95;
  Float_t emin=3.5;
  Int_t plotPoints=12;
  Float_t stepsize=.1;
  Float_t fvCut;
  Float_t fvStep=15;
  Int_t numFVsteps=3;

  //read in the file names
  ifstream runList;
  runList.open("runlist_salt.txt");
  numfiles=0;
  while(!runList.eof()){
    runList>>files[numfiles];
    numfiles++;
  }
  runList.close();
  
  //make the histogram
  TH1F * h1= new TH1F("h1","Accidental rate as function of energy FV cuts from 555-600 in Salt",30,3,6);
  TH1F * h2= new TH1F("h2","Accidental rate as function of energy for Nhit >20 and fv of 560",30,3,6);
  TH1F * h3= new TH1F("h3","Accidental rate as function of energy for Nhit >20 and fv of 570",30,3,6);
  TH1F * h4= new TH1F("h4","Accidental rate as function of energy for Nhit >20 and fv of 580",30,3,6);
  TH1F * h5= new TH1F("h5","Accidental rate as function of energy for Nhit >20 and fv of 590",30,3,6);
  TH1F * h6= new TH1F("h6","Accidental rate as function of energy for Nhit >20 and various FV Cuts",30,3,6);
  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(3);
  h4->SetLineColor(4);
  h5->SetLineColor(5);
  h6->SetLineColor(6);
  h1->SetYTitle("Number of accidentals");
  h1->SetXTitle("Energy [MeV]");

  //loop over different fidcucial volume cuts
  for(int h=0;h<=numFVsteps;h++){ 
  fvCut= piFV + fvStep*h;
  numcl=0;

  //loop over all ntuples
  for(Int_t i=0; i<(numfiles-1);i++){
    //read in the ntuples
    TFile *f1 = new TFile(files[i]);
    TTree *Td = (TTree*)f1->Get("h330");
    nd = Td->GetEntries();
    if(h==0){
      cerr<<nd<<" events in file "<<i<<endl;
    }
    //reset the number of clean events to 0
    numClean[i]=0;
    estep=emin;
    if(h==0&& i==0){
      livetime=0;
    }

    //Pull info from ntuple
    Td->SetBranchAddress("Nhits",&nhits);
    Td->SetBranchAddress("Enep",&ene);
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
    Td->SetBranchAddress("Beta1p",&beta1);
    Td->SetBranchAddress("Beta4p",&beta4);
    Td->SetBranchAddress("Rdamn1",&rdamn1);
    Td->SetBranchAddress("Rdamn2",&rdamn2);
    
    //loop over mininum accepted energies
    for(int j=0;j<plotPoints;j++){
      //get first entry
      Td->GetEntry(0);
      //initialize the time variables
      if(j==0){
	dayi=jdy;
	totaltimei=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
      }
      //initialize the counted variable
      numpi[i][j]=0;
      numtag[h][j]=0;
      //loop over events in ntuple
      for(int k=0;k<nd;k++){
	Td->GetEntry(k);
	totaltime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
	day=jdy;
	beta14=beta1+4*beta4;
	//set the damn flag use UNIDOC
	if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
	  damnFlag =1;  
	}
	else {
	  damnFlag=0;
	}
	
	//set the clean flag
	if(nhits>nhitmin && damnFlag==1 && rf<cleanFV){
	  cleanFlag=1;
	  if(j==0 && h==0){
	    numcl++;
	  }
	}
	else cleanFlag=0;
	
	//determine if event is a primary
	if( damnFlag==1 && cleanFlag==1 && ene>estep && rf<fvCut /* addition cuts */){
	  piFlag=1;
	  numpi[i][j]++;
	  numtag[h][j]++;
	}

      }//end loop over events in ntuple
      //find the total runtime on only the first iteration of the stepping loop
      if(j==0 && h==0){
	totaldays=day-dayi;
	if(totaldays==0){
	  runtime=totaltime-totaltimei;
	}
	else if(totaldays>0){
	  runtime= (86400.0-totaltimei)+(totaldays-1)*86400.0+ totaltime;
	}
	else runtime=0;
	
	livetime+=runtime;
      }
      estep=estep+stepsize;
    }//end loop over minimum accepted energies
 
  }//end loop over ntuples

 
  //print out acceptence rates and fill histogram

  if(h==0){
    if(livetime !=0.0){
      cleanRate=numcl/livetime;
    }
    else cleanRate=999;
    cerr<<"The rate of clean events was "<<cleanRate<<" in a livetime of "<<livetime<<" s"<<endl;
  }
  cerr<<"There were "<<numtag[h][5]<<" tagged events at energy "<<emin+stepsize*5<<endl;
    for(int m=0;m<plotPoints;m++){
      cerr<<numtag[h][m]<<" tagged events at fv cut of "<<fvCut<<" and energy "<<emin+stepsize*m<<endl;
      if(h==0){
	h1->Fill(emin+stepsize*m,numtag[h][m]*cleanRate*timeWindow/livetime);
      }
      else if(h==1){
	h2->Fill(emin+stepsize*m,numtag[h][m]*cleanRate*timeWindow/livetime);
      }
      else if(h==2){
	h3->Fill(emin+stepsize*m,numtag[h][m]*cleanRate*timeWindow/livetime);
      }
      else if(h==3){
	h4->Fill(emin+stepsize*m,numtag[h][m]*cleanRate*timeWindow/livetime);
      }
      else if(h==4){
	h5->Fill(emin+stepsize*m,numtag[h][m]*cleanRate*timeWindow/livetime);
      }
      else if(h==5){
	h6->Fill(emin+stepsize*m,numtag[h][m]*cleanRate*timeWindow/livetime);
      }
    }

 
  if(h==0){
    h1->Draw();
  }
  else if( h==1){
    h2->Draw("same");
  }
  else if( h==2){
  
    h3->Draw("same");
  }
  else if( h==3){

    h4->Draw("same");
  }
  else if( h==4){

    h5->Draw("same");
  }
  else if( h==5){
    h6->Draw("same");
  }
  }//end loop over different fv cuts

}//end code
