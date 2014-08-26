//This code looks at atmospheric neutrino monte carlo and counts 2 and three fold coincidences

{//start code
  gROOT->SetStyle("Plain");
  ifstream runList;
  int numfiles;
  TString files[1000];
  Float_t nhits;
  Float_t ene;
  Float_t jdy;
  Float_t dt;
  Float_t dt0;
  Float_t ut3;
  Float_t rf;
  Float_t xfit;
  Float_t yfit;
  Float_t zfit;
  Float_t xi;
  Float_t yi;
  Float_t zi;
  Float_t xp[1000];
  Float_t yp[1000];
  Float_t zp[1000];
	Float_t dtp[1000];
  Float_t theta;
  Float_t phi;
  Float_t thetaij;
  Float_t evid;
  Float_t rdamn1;
  Float_t rdamn2;
  Float_t itr;
  Float_t neutron;
  Float_t xgen;
  Float_t ygen;
  Float_t zgen;
  Float_t rgen;
  Float_t egen;
  Float_t enerau;
  Float_t beta1;
  Float_t beta14;
  Float_t mult;
  Float_t frati;
  Float_t nburst;
  Float_t retrig;
  Int_t mb=0;
	Int_t bList[100];
	Int_t piNum;
  
  

  Float_t farray[12];

  Int_t gt[1000000]={0.0};
  Int_t index[1000000]={0.0};
  Int_t min;
  Int_t numEV =0;
  //out.open("mc_sort.root");
  //cut
//   Float_t nhitCut=20;
//   Float_t nhitMax=150;
//   Float_t fiducialVolume=600;
//   Float_t timeWindow=0.150;
//   Float_t distanceWindow=500;
//   Float_t startWindow=100;
//   Float_t betamin=-0.12;
//   Float_t betamax=0.95;
//   Float_t thetamin=0.75;
//   Float_t thetamax=1.45;
//   Float_t itrmin= 0.55;
//   Float_t itrmax=0.95;
//   Float_t emin=4.5+.511;
//   Float_t ethresh=3.9+.511;
//   Float_t outerVolume=650;
//   Float_t innerVolume=550;


  Int_t nhitmax=150;
  Int_t nhitmin=25;
  Float_t emin=4.5+.511;
  Float_t etrig=3.9+.511;
  Float_t OV=650.0;
  Float_t IV=550.0;
  Float_t FV=600.0;
  Float_t betamin=-.12;
  Float_t betamax=.95;
  Float_t thetamin=.75;
  Float_t thetamax=1.45;
  Float_t itrmin=.55;
  Float_t itrmax=.95;
  Float_t distanceWindow=500.0;
  Float_t timeWindow=.150;
  //flags
  Int_t cleanFlag=0;
  Int_t damnFlag=0;
  Int_t damnFlag2=0;
  Int_t hlcFlag=0;
  Int_t fratiFlag=0;
  Int_t piFlag=0;
  Int_t posIsPi=0;
  Int_t posFlag=0;
  //primary information
  Int_t PiID[1000000]={0.0};
  Float_t PiEn[1000000]={0.0};
  Int_t fratiBit[1000000]={0.0};
  Int_t tempPiID=0;
  Float_t tempPiEn=0;
  Int_t posID=0;
  Float_t posTime=0;
  Float_t fTime=0;
  Float_t tTime=0;
  Int_t lastID=0;
  Int_t numMMF1=0;
  Int_t numMMF2=0;
  //coincidence info
  Int_t coinNum=0;
  Int_t num2fold=0;
  Int_t num3fold=0;
  Int_t numClean=0;
  Int_t bursts=0;
  Int_t numPred[1000000]={0.0}; 
  Int_t preds=0;
  Int_t PeID[1000000][50]={0.0};
  Float_t PeEn[1000000][50]={0.0};
  Float_t PeDt[1000000][50]={0.0};
  Float_t PeDr[1000000][50]={0.0};
  Float_t tempR;
  Float_t dR;
  Float_t dR1;
  Float_t dR2;
  Float_t dR3[3];

  Float_t timediff;

  //TFile *f1 = new TFile("~shokair/snodata/atm/forMelissa/atm_d2o_leta.root");
  //TFile *f1 = new TFile("/home/shokair/ntuples/sort/atm_d2o_leta_sort.root");
  	TFile *f1 = new TFile("/Volumes/shokair@nu.hep.upenn.edu/ntuples/sort/atm_d2o_leta_sort.root");

	TTree *Td = (TTree*)f1->Get("h102");
  //TTree *Td1 = (TTree*)f1->Get("h201");

  Int_t nd = Td->GetEntries();
  cerr<<nd<<" events in file "<<0<<endl;
  h1= new TH1F("h1","Delta R Atmospheric BG D2O",1300,0,1300);
  h2= new TH1F("h2","Beta14",400,-2,2);
  h3= new TH1F("h3","ITR",400,-2.0,2);
  h4= new TH1F("h4","ThetaIJ",400,-2,2);
  h5= new TH1F("h5","Multiplicity D2O",20,0,20);
  h6= new TH1F("h6","Delta t",100,0,.1);
  h7= new TH1F("h7","Nhits",100,0,100);
  h8= new TH1F("h8","Rgen",100,0,100);
  h9= new TH1F("h9","Delta R 2D(XY)",1300,0,1300);
  h10= new TH1F("h10","Frati Time",100,0,1);
  h11= new TH1F("h11","Rfit",650,0,650); 
  h12= new TH1F("h12","Second particle travel length",1300,0,1300);
  h13= new TH1F("h13","Pos energy",200,0,20);
  h14= new TH1F("h14","neutron energy",200,0,20);
  TH2F * p= new TH2F("p"," Predessor Time Separation vs Predecessor Space Separation",250,0,500,45,0,.15);
  h1->SetXTitle("delta r [cm]");
  h1->SetLineColor(4);
  h9->SetXTitle("delta r [cm]");
  h9->SetLineColor(3);
  h5->SetXTitle("Multiplicty of Burst");

  //get the necessary information from ntuple
  Td->SetBranchAddress("Nhits",&nhits);
  Td->SetBranchAddress("Enep",&ene);
  Td->SetBranchAddress("Itrp",&itr);
  Td->SetBranchAddress("Ev_tid",&evid);	
  Td->SetBranchAddress("Beta_14",&beta14);
  Td->SetBranchAddress("Rfp",&rf);
  Td->SetBranchAddress("T_dif",&dt); //julian day
  Td->SetBranchAddress("T0_dif",&dt0); //time in seconds
  Td->SetBranchAddress("Nu_ener",&egen);
  Td->SetBranchAddress("Thetaijp",&thetaij);
  Td->SetBranchAddress("E_rdamn1",&rdamn1);
  Td->SetBranchAddress("E_rdamn2",&rdamn2);
  Td->SetBranchAddress("E_thetai",&theta);
  Td->SetBranchAddress("E_phiprb",&phi);
  Td->SetBranchAddress("Pid1",&xfit);
  Td->SetBranchAddress("Pid2",&yfit);;
  Td->SetBranchAddress("Pid3",&zfit);
  Td->SetBranchAddress("Multnc",&mult);
  Td->SetBranchAddress("Frati_f",&frati);
  Td->SetBranchAddress("Burst",&nburst);
  Td->SetBranchAddress("Retrig",&retrig);

  Float_t eI;
  Int_t endJ=0;
  Int_t setFlag=0;
  Int_t tMult;
  Int_t tBit=0;
  Int_t tNhits;

  for (Int_t j=0;j<nd;j++){
    Td->GetEntry(j);
    //Td1->GetEntry(j);
    tTime=dt0;
    fTime=0;
    tMult=mult;
    tNhits=nhits;
    Int_t n;
    Int_t lastN;
    tBit=frati;
    eI=egen;
    //cerr<<frati<<" "<<nhits<<" "<<dt0<<endl;
    //cerr<<dt0<<" "<<j<<endl;
    if (nhits>nhitmax){
      //find start of burst
      //fratiBit[j]=1;
      n=1;
      //cerr<<j<<": "<<evid<<" "<<dt0<<endl;
      while(n<=tMult+50){
	Td->GetEntry(j+n);
	//cerr<<j-n<<": "<<evid<<" "<<dt0<<endl;
	if(tTime==99999 && egen==eI){
	  fTime=dt0;
	}
	else if (egen==eI){
	  fTime=dt0-tTime;
	}
	h10->Fill(fTime);
	if(fTime<=.250 && fratiBit[j+n]!=1 && egen==eI){
	  fratiBit[j+n]=1;
	  numMMF1++;
	  //cerr<<j-n<<" "<<evid<<" "<<egen<<" "<<eI<<endl;
	  
	  if(Int_t(frati)!=1){

	    //cerr<<j-n<<" "<<evid<<" "<<egen<<" "<<eI<<" "<<tNhits<<" "<<fTime<<endl;
	  }
	  
	}
	n++;
      }

    }

  }



  cerr<<numMMF1<<endl;
  numMMF1=0;

  Int_t numDiff=0;
  for (Int_t i=0;i<nd;i++){
    Td->GetEntry(i);
    //Td1->GetEntry(i);
    h7->Fill(nhits);
    fratiFlag=0;
    //mb=0;
    //cerr<<fratiBit[i]<<endl;
    if(Int_t(frati)-fratiBit[i]!=0){
      numDiff++;
      //cerr<<Int_t(frati)-fratiBit[i]<<endl;
    }
    if(frati==1){
      numMMF2++;
      //cerr<<i<<endl;
    }
    if(fratiBit[i]==1){
      numMMF1++;
      //cerr<<i<<endl;
    }
    //cerr<<"Delta T: "<<dt0<<" ID: "<<evid<<" energy: "<<ene<<" multiplicity: "<<mult<<endl; 
    //cerr<<"Itr: " <<itr<<" Beta14: "<<beta14<<" Thetaij: "<<thetaij<<endl;
    /*
    if(nhits>nhitmax){
      int p=i+1;
      tTime=dt0;
      fTime=0;
      while(fTime<.250 && fTime>=0 && dt0!=99999){
	Td->GetEntry(p);
	fTime=tTime-dt0;
	//cerr<<fTime<<" "<<i<<" "<<p<<endl;
	p++;
      }
      if(p>i+2){
	i=p;
      }
      else{ Td->GetEntry(i);}
      fratiFlag=1;
    }
    */  
    if(retrig!=1 &&nburst==0 && fratiBit[i]==0){
	damnFlag = 1;  
	h2->Fill(beta14);
	h3->Fill(itr);
	h4->Fill(theta);
	
    }
    if(itr<itrmax && itr>itrmin && theta>thetamin && theta<thetamax ){
      hlcFlag=1;
    }
    else{hlcFlag=0;}
    if( rf>IV && rf<OV && ene>etrig && damnFlag==1 && hlcFlag==1 && nhits<nhitmax){
      cleanFlag=1; 
      numClean++;
      //cerr<<dt0<<endl;
      //cerr<<"Itr: " <<itr<<" Beta14: "<<beta14<<" Thetaij: "<<thetaij/100<<endl;
    }
    else{ 
      cleanFlag=0;
      //cerr<<"Itr: " <<itr<<" Beta14: "<<beta14<<" Energy : "<<ene<<" Nhits: "<<nhits<<" Rf: "<<rf<<endl;

      //cerr<<"NOT CLEAN!"<<endl;
    }
    if(dt0==99999){
      //cerr<<"Delta T: "<<dt0<<" ID: "<<evid<<" energy: "<<ene<<" multiplicity: "<<mult<<endl;
      tempPiEn=ene;
      tempR=rf;
      preds=0;
      piFlag=0;
      posIsPi=0;
      posFlag=0;
      mb=0;
      piNum=999;
      tempPiID= -999;

      if(ene<emin && cleanFlag==1){
		  xp[mb]=xfit;
		  yp[mb]=yfit;
		  zp[mb]=zfit;
		  dtp[mb]=0;
		  mb++;
		  posID= evid;
		  posTime=0;
		  posFlag=1;
	
      }
	  if(ene>emin && cleanFlag==1){
		xp[mb]=xfit;
		yp[mb]=yfit;
		zp[mb]=zfit;
		  dtp[mb]=0;
		  piNum=mb;
		mb++;
		tempPiID= evid;
		xi=xfit;
		yi=yfit;
		zi=zfit;
		piFlag=1;
		posIsPi=1;
	}

	
      //dt0=0;
      do{
	i++;
	Td->GetEntry(i);
	//Td1->GetEntry(i);
	if(frati==1){
	  numMMF2++;
	  //cerr<<i<<endl;
	}
	if(fratiBit[i]==1){
	  numMMF1++;
	  //cerr<<i<<endl;
	}
    //
	if(retrig!=1 &&nburst==0 && fratiBit[i]==0){
	  damnFlag2=1;  
	  h6->Fill(dt0);
	}
	else {
	  damnFlag2=0;
	}
	if(itr>itrmin && itr<itrmax && theta>thetamin && theta<thetamax ){
	  hlcFlag=1;
	}
	else{hlcFlag=0;}
	if( rf<FV && ene>etrig && hlcFlag==1 && damnFlag2==1 &&nhits<nhitmax ){
	  
	  /*
	  xfit=rf*TMath::Sin(theta)*TMath::Cos(phi);
	  yfit=rf*TMath::Sin(theta)*TMath::Sin(phi);
	  zfit=rf*TMath::Cos(theta);
	  */
	  numClean++;
	      
	  //if(dt0<timeWindow ){
	    
	    // h1->Fill(dR);
	    //h9->Fill(dR2);
	    if(piFlag==0 && ene>emin){
	      	xi=xfit;
		yi=yfit;
		zi=zfit;
		piFlag=1;
		piTime=dt0;
		xp[mb]=xfit;
		yp[mb]=yfit;
		zp[mb]=zfit;
			dtp[mb]=dt0;
		piNum=mb;
		mb++;
	    }
	    else{
	      xp[mb]=xfit;
	      yp[mb]=yfit;
	      zp[mb]=zfit;
			dtp[mb]=dt0;
	      mb++;
	    }
	    PeID[coinNum][preds]=evid;
	    PeEn[coinNum][preds]=ene;
	    PeDt[coinNum][preds]=dt0;
	    PeDr[coinNum][preds]=dR;
	    h11->Fill(rf);
	    //preds++;
	    mb++;

	    //cerr<<"r: "<<rf<<" theta: "<<theta<<" phi: "<<phi<<endl;
	    //cerr<<"("<<xfit<<", "<<yfit<<", "<<zfit<<")"<<endl;
	    //cerr<<"Thetaij: "<<thetaij/100<<endl;
	    //cerr<<evid <<" is pred with nhit:"<<nhits<<" energy "<<ene<<" dt "<<dt0<<" for "<<tempPiID<<endl;
	 // }
	}

      }while(dt0!=99999 && i<nd);//end while
      //cerr<<tempPiID<<" had "<<numPred[coinNum]<<" predecessors"<<endl;
      i=i-1;
		preds=0;
  for (int b=0; b<mb; b++) {
	  if (b!=piNum) {
		  if (TMath::Abs(dtp[b]-dtp[piNum])) {
			  bList[preds]=b;
			  preds++;
		  }
	  }
  }
	if(piFlag==1){
	h5->Fill(preds);
      
	if(preds>0){
	  bursts++;
	}
	if(preds==2 ||preds==1){
	  coinNum++;
	  //cerr<<tempPiID<<" had "<<preds<<"predecessors"<<endl;

	}
	if(preds==1){
	  for(int ct=0;ct<preds;ct++){
	    if(ct!= piNum){
	      dR=TMath::Sqrt((xi-xp[bList[ct]])*(xi-xp[bList[ct]])+(yi-yp[bList[ct]])*(yi-yp[bList[ct]])+(zi-zp[bList[ct]])*(zi-zp[bList[ct]]));
	    }
	  }
	  h1->Fill(dR);
	  if( dR<distanceWindow){
	    num2fold++;
	    h12->Fill(PeDr[coinNum-1][0]);
	    p->Fill(PeDr[coinNum-1][0],PeDt[coinNum-1][0]);
	    h13->Fill(tempPiEn);
	  }
	}
	if(preds==2){
	  for(int ct=0;ct<preds;ct++){
	    if(ct!= piNum){
	      dR3[ct]=TMath::Sqrt((xi-xp[bList[ct]])*(xi-xp[bList[ct]])+(yi-yp[bList[ct]])*(yi-yp[bList[ct]])+(zi-zp[bList[ct]])*(zi-zp[bList[ct]]));
	    }
	    else{
	      dR3[ct]=0;
	    }
	    h1->Fill(dR3[ct]);
	  }
	  if(dR3[0]<distanceWindow && dR3[0]<distanceWindow && dR3[0]<distanceWindow){
	    num3fold++;
	    p->Fill(PeDr[coinNum-1][0],PeDt[coinNum-1][0]);
	    p->Fill(PeDr[coinNum-1][1],PeDt[coinNum-1][1]);
	    h13->Fill(tempPiEn);
	    h12->Fill(PeDr[coinNum-1][0]);
	    //cerr<<"3-fold event "<<endl;
	    //cerr<<tempPiID<<" "<<tempPiEn<<endl;
	    //cerr<<PeID[coinNum-1][0]<<" "<<PeEn[coinNum-1][0]<<" "<<PeDt[coinNum-1][0]<<endl;
	    //cerr<<PeID[coinNum-1][1]<<" "<<PeEn[coinNum-1][1]<<" "<<PeDt[coinNum-1][1]<<endl;
	  }
	}
	if(mb>3){
	  //cerr<<mb<<"-fold event"<<endl;
	}
      }
    }//end actions on clean
  }//end loop over events
  cerr<<"Frati cut removed "<<numMMF2<<" my version removed "<<numMMF1<<" -  "<<numDiff<<"discrepencies "<<endl;
  cerr<<bursts<<" bursts"<<endl;
  cerr<<coinNum<<" coincidences"<<endl;
  cerr<<num2fold<<" 2-fold"<<endl;
  cerr<<num3fold<<" 3-fold"<<endl;

  h5->Draw();
}//end code
