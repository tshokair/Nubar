//Code to calculate the accidental coincidence rates from the Salt Phase of SNO
//latest modification on 5/18/2011 
//no longer using CINT but compiling c++
//changed the getRateSalt.cpp code so that it is quicker 6/9/11
//7/7/11 look at more OV cuts
//8/22/11 took out the loop over fv cuts to make the code run faster



//begin
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGaxis.h>
#include <TApplication.h>
#include <TLegend.h>

using namespace std;

int main(int argc, char** argv)
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
  Int_t numClean[2000];
  Int_t numpi0[2000][50];
  Int_t numpi1[2000][50];
  Int_t numpi2[2000][50];
  Int_t numpi3[2000][50];
  Float_t estep;
  Int_t numcl[10];
  Int_t numtag0[50][50]={0.0};
  Int_t numtag1[50][50]={0.0};
  Int_t numtag2[50][50]={0.0};
  Int_t numtag3[50][50]={0.0};
  Double_t livetime;
  Int_t tC;
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
  Float_t cleanRate[10];
  Float_t rate[10000][50];
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
  Int_t piFlag[50];
  Int_t cleanFlag0;
  Int_t cleanFlag1;
  Int_t cleanFlag2;
  Int_t cleanFlag3;
  //array for file names
  TString files[10000]={0.0};
  Int_t numfiles=0;
  Float_t scale;
  //const char *hist[10];
  
  //declare and initialize cut parameters  Float_t nhitmin=20;
  Float_t nhitmax=160;
  Float_t piFV=550;
  Float_t cleanFV=675;
  Float_t timeWindow=.150;
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
  Int_t numFVsteps=4;
  Float_t etrig=2.7;
  Float_t avStart=599.5;
  Float_t avEnd=605;
  Float_t cleanStep=25;
  
  //read in the file names
  ifstream rL;
  rL.open("b8filesd2o.runlist");
  //runList.open("runlist.txt");
  numfiles=0;
  while(!rL.eof()){
    rL>>files[numfiles];
    numfiles++;
  }
  rL.close();
  TApplication* rootapp = new TApplication("ex", &argc, argv);
  //make the histogram(rf<avStart ||(rf>avEnd && rf<cleanFV))
  std::stringstream caph1;
  std::stringstream caph2;
  std::stringstream caph3;
  std::stringstream caph4;
  caph1<<"Accidental rates for an Outer Volume of "<<cleanFV - cleanStep*0<<" cm in D2O";
  caph2<<"Accidental rates for an Outer Volume of "<<cleanFV - cleanStep*1<<" cm in D2O";
  caph3<<"Accidental rates for an Outer Volume of "<<cleanFV - cleanStep*2<<" cm in D2O";
  caph4<<"Accidental rates for an Outer Volume of "<<cleanFV - cleanStep*3<<" cm in D2O";
  TString sth1=caph1.str();
  TString sth2=caph2.str();
  TString sth3=caph3.str();
  TString sth4=caph4.str();
  
  TH1F * h1= new TH1F("h1",sth1,24,3,6);
  TH1F * h2= new TH1F("h2",sth1,24,3,6);
  TH1F * h3= new TH1F("h3",sth1,24,3,6);
  TH1F * h4= new TH1F("h4",sth1,24,3,6);
  TH1F * h5= new TH1F("h5",sth1,24,3,6);
  TH1F * h6= new TH1F("h6",sth1,24,3,6);

  TH1F * h7= new TH1F("h7","Accidental events as function of energy FV cuts from 555-600 in D2o OV=650",24,3,6);
  TH1F * h8= new TH1F("h8","Accidental rate as function of energy for Nhit >20 and fv of 560",24,3,6);
  TH1F * h9= new TH1F("h9","Accidental rate as function of energy for Nhit >20 and fv of 570",24,3,6);
  TH1F * h10= new TH1F("h10",sth2,24,3,6);
  TH1F * h11= new TH1F("h11","Accidental rate as function of energy for Nhit >20 and fv of 590",24,3,6);
  TH1F * h12= new TH1F("h12","Accidental rate as function of energy for Nhit >20 and various FV Cuts",24,3,6);

  TH1F * h13= new TH1F("h13","Accidental events as function of energy FV cuts from 555-600 in D2o OV=625",24,3,6);
  TH1F * h14= new TH1F("h14","Accidental rate as function of energy for Nhit >20 and fv of 560",24,3,6);
  TH1F * h15= new TH1F("h15","Accidental rate as function of energy for Nhit >20 and fv of 570",24,3,6);
  TH1F * h16= new TH1F("h16",sth3,24,3,6);
  TH1F * h17= new TH1F("h17","Accidental rate as function of energy for Nhit >20 and fv of 590",24,3,6);
  TH1F * h18= new TH1F("h18","Accidental rate as function of energy for Nhit >20 and various FV Cuts",24,3,6);

  TH1F * h19= new TH1F("h19","Accidental events as function of energy FV cuts from 555-600 in D2o OV=600",24,3,6);
  TH1F * h20= new TH1F("h20","Accidental rate as function of energy for Nhit >20 and fv of 560",24,3,6);
  TH1F * h21= new TH1F("h21","Accidental rate as function of energy for Nhit >20 and fv of 570",24,3,6);
  TH1F * h22= new TH1F("h22",sth4,24,3,6);
  TH1F * h23= new TH1F("h23","Accidental rate as function of energy for Nhit >20 and fv of 590",24,3,6);
  TH1F * h24= new TH1F("h24","Accidental rate as function of energy for Nhit >20 and various FV Cuts",24,3,6);

  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(3);
  h4->SetLineColor(4);
  h5->SetLineColor(5);
  h6->SetLineColor(6);
  h1->SetYTitle("Number of accidentals [Hz]");
  h1->SetXTitle("Energy [MeV]");

  h7->SetLineColor(1);
  h8->SetLineColor(2);
  h9->SetLineColor(3);
  h10->SetLineColor(4);
  h11->SetLineColor(5);
  h12->SetLineColor(6);
  h7->SetYTitle("Number of accidentals [Hz]");
  h7->SetXTitle("Energy [MeV]");

  h13->SetLineColor(1);
  h14->SetLineColor(2);
  h15->SetLineColor(3);
  h16->SetLineColor(4);
  h17->SetLineColor(5);
  h18->SetLineColor(6);
  h13->SetYTitle("Number of accidentals [Hz]");
  h13->SetXTitle("Energy [MeV]");

  h19->SetLineColor(1);
  h20->SetLineColor(2);
  h21->SetLineColor(3);
  h22->SetLineColor(4);
  h23->SetLineColor(5);
  h24->SetLineColor(6);
  h19->SetYTitle("Number of accidentals [Hz]");
  h19->SetXTitle("Energy [MeV]");


  numcl[0]=0;
  numcl[1]=0;
  numcl[2]=0;
  numcl[3]=0;

  //loop over all ntuples
  for(Int_t i=0; i<(numfiles-1);i++){
    //read in the ntuples
    TFile *f1 = new TFile(files[i]);
    TTree *Td = (TTree*)f1->Get("h360");
    nd = Td->GetEntries();
    tC+=nd;
    
      cerr<<nd<<" events in file "<<i<<" and "<<tC<<" events total"<<endl;
    
    //reset the number of clean events to 0
    numClean[i]=0;
    //estep=emin;
    if(i==0){
      livetime=0;
    }

    //Pull info from ntuple
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
    Td->SetBranchAddress("Beta1p",&beta1);
    Td->SetBranchAddress("Beta4p",&beta4);
    Td->SetBranchAddress("Rdamn1",&rdamn1);
    Td->SetBranchAddress("Rdamn2",&rdamn2);
    

 
      //get first entry
      Td->GetEntry(0);
      //initialize the time variables
	dayi=jdy;
	totaltimei=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
      //initialize the counted variable
	for (int h=0;h<numFVsteps;h++){
	  for(int j=0;j<plotPoints;j++){
	    piFlag[j]=0;
	    numpi0[i][j]=0;
	    numtag0[h][j]=0;
	    numpi1[i][j]=0;
	    numtag1[h][j]=0;
	    numpi2[i][j]=0;
	    numtag2[h][j]=0;
	    numpi3[i][j]=0;
	    numtag3[h][j]=0;
	  }
	}
	//loop over events in ntuple
      for(int k=0;k<nd;k++){
	Td->GetEntry(k);
	totaltime=(ut1)+(ut2)/1.E6+(ut3)/1.E9;
	day=jdy;
	beta14=beta1+4*beta4;
	cleanFlag0=0;cleanFlag1=0;cleanFlag2=0;cleanFlag3=0;
	//set the damn flag use UNIDOC
	if((Int_t(rdamn1)& 0x01b5)==0 && (Int_t(rdamn2)&& 0x6FE1)==0){
	  damnFlag =1;  
	}
	else {
	  damnFlag=0;
	}
	
	//set the clean flag
	
	if(ene>etrig && damnFlag==1){
	  if (rf<cleanFV){
	    cleanFlag0=1;
	    numcl[0]++;

	  }
	  if(rf<cleanFV-cleanStep*1){
	    cleanFlag1=1;
	    numcl[1]++;

	  }
	  if(rf<cleanFV-cleanStep*2){
	    cleanFlag2=1;
	    numcl[2]++;

	  }
	  if(rf<cleanFV-cleanStep*3){
	    cleanFlag3=1;
	    numcl[3]++;
	  }
	}
	else {cleanFlag0=0;cleanFlag1=0;cleanFlag2=0;cleanFlag3=0;}
	
	//determine if event is a primary given a series of energies
	for(int h=0;h<numFVsteps;h++){ 
	  for(int j=0;j<plotPoints;j++){
	    
	    fvCut= piFV + fvStep*h;
	    if(cleanFlag0==1 && ene>emin+stepsize*j && rf<fvCut /* addition cuts */){
	      piFlag[j]=1;
	      numpi0[i][j]++;
	      numtag0[h][j]++;
	    }
	    if(cleanFlag1==1 && ene>emin+stepsize*j && rf<fvCut /* addition cuts */){
	      piFlag[j]=1;
	      numpi1[i][j]++;
	      numtag1[h][j]++;
	    }
	    if(cleanFlag2==1 && ene>emin+stepsize*j && rf<fvCut /* addition cuts */){
	      piFlag[j]=1;
	      numpi2[i][j]++;
	      numtag2[h][j]++;
	    }
	    if(cleanFlag3==1 && ene>emin+stepsize*j && rf<fvCut /* addition cuts */){
	      piFlag[j]=1;
	      numpi3[i][j]++;
	      numtag3[h][j]++;
	    }
	  }
	}
	
      }//end loop over events in ntuple
      //find the total runtime on only the first iteration of the stepping loop
      
	totaldays=day-dayi;
	if(totaldays==0){
	  runtime=totaltime-totaltimei;
	}
	else if(totaldays>0){
	  runtime= (86400.0-totaltimei)+(totaldays-1)*86400.0+ totaltime;
	}
	else runtime=0;
	
	livetime+=runtime;
      

    delete Td;
    delete f1;
  }//end loop over ntuples

 
  //print out acceptence rates and fill histogram

 
    if(livetime !=0.0){
      cleanRate[0]=numcl[0]/livetime;
      cleanRate[1]=numcl[1]/livetime;
      cleanRate[2]=numcl[2]/livetime;
      cleanRate[3]=numcl[3]/livetime;
      
    }
    else {cleanRate[0]=999;cleanRate[1]=999;cleanRate[2]=999;cleanRate[3]=999;}
    cerr<<"The clean rate for OV= "<<cleanFV-cleanStep*0<<"  was "<<cleanRate[0]<<" in a livetime of "<<livetime<<" s"<<endl;
    cerr<<"The clean rate for OV= "<<cleanFV-cleanStep*1 << " was "<<cleanRate[1]<<endl;
    cerr<<"The clean rate for OV= "<<cleanFV-cleanStep*2 <<" was "<<cleanRate[2]<<endl;
    cerr<<"The clean rate for OV= "<<cleanFV-cleanStep*3 <<"  was "<<cleanRate[3]<<endl;
  
  //cerr<<"There were "<<numtag0[h][5]<<" tagged events at energy "<<emin+stepsize*5<<endl;
    for(int m=0;m<plotPoints;m++){
      /*
      cerr<<numtag0[h][m]<<" tagged events at fv cut of "<<fvCut<<" and energy "<<emin+stepsize*m<<endl;
      cerr<<numtag1[h][m]<<" tagged events at fv cut of "<<fvCut<<" and energy "<<emin+stepsize*m<<endl;
      cerr<<numtag2[h][m]<<" tagged events at fv cut of "<<fvCut<<" and energy "<<emin+stepsize*m<<endl;
      cerr<<numtag3[h][m]<<" tagged events at fv cut of "<<fvCut<<" and energy "<<emin+stepsize*m<<endl;
      */
      
	h1->Fill(emin+stepsize*m,numtag0[0][m]*cleanRate[0]*timeWindow);
	h7->Fill(emin+stepsize*m,numtag1[0][m]*cleanRate[1]*timeWindow);
	h13->Fill(emin+stepsize*m,numtag2[0][m]*cleanRate[2]*timeWindow);
	h19->Fill(emin+stepsize*m,numtag3[0][m]*cleanRate[3]*timeWindow);
      
      
	h2->Fill(emin+stepsize*m,numtag0[1][m]*cleanRate[0]*timeWindow);
	h8->Fill(emin+stepsize*m,numtag1[1][m]*cleanRate[1]*timeWindow);
	h14->Fill(emin+stepsize*m,numtag2[1][m]*cleanRate[2]*timeWindow);
	h20->Fill(emin+stepsize*m,numtag3[1][m]*cleanRate[3]*timeWindow);
      
      
	h3->Fill(emin+stepsize*m,numtag0[2][m]*cleanRate[0]*timeWindow);
	h9->Fill(emin+stepsize*m,numtag1[2][m]*cleanRate[1]*timeWindow);
	h15->Fill(emin+stepsize*m,numtag2[2][m]*cleanRate[2]*timeWindow);
	h21->Fill(emin+stepsize*m,numtag3[2][m]*cleanRate[3]*timeWindow);
      
      
	h4->Fill(emin+stepsize*m,numtag0[3][m]*cleanRate[0]*timeWindow);
	h10->Fill(emin+stepsize*m,numtag1[3][m]*cleanRate[1]*timeWindow);
	h16->Fill(emin+stepsize*m,numtag2[3][m]*cleanRate[2]*timeWindow);
	h22->Fill(emin+stepsize*m,numtag3[3][m]*cleanRate[3]*timeWindow);
      
      
	h5->Fill(emin+stepsize*m,numtag0[4][m]*cleanRate[0]*timeWindow);
	h11->Fill(emin+stepsize*m,numtag1[4][m]*cleanRate[1]*timeWindow);
	h17->Fill(emin+stepsize*m,numtag2[4][m]*cleanRate[2]*timeWindow);
	h23->Fill(emin+stepsize*m,numtag3[4][m]*cleanRate[3]*timeWindow);
      
      
	h6->Fill(emin+stepsize*m,numtag0[5][m]*cleanRate[0]*timeWindow);
	h12->Fill(emin+stepsize*m,numtag1[5][m]*cleanRate[1]*timeWindow);
	h18->Fill(emin+stepsize*m,numtag2[5][m]*cleanRate[2]*timeWindow);
	h24->Fill(emin+stepsize*m,numtag3[5][m]*cleanRate[3]*timeWindow);
      
    }
  
    std::stringstream cap1;
    std::stringstream cap2;
    std::stringstream cap3;
    std::stringstream cap4;
    cap1<<piFV + fvStep*0<<" cm";
    cap2<<piFV + fvStep*1<<" cm";
    cap3<<piFV + fvStep*2<<" cm";
    cap4<<piFV + fvStep*3<<" cm";
    TString st1=cap1.str();
    TString st2=cap2.str();
    TString st3=cap3.str();
    TString st4=cap4.str();
    TLegend *leg = new TLegend(0.6,0.7,0.89,0.89);
    leg->AddEntry(h4,st4,"l");
    leg->AddEntry(h3,st3,"l");
    leg->AddEntry(h2,st2,"l");
    leg->AddEntry(h1,st1,"l");
    TCanvas *c1 =new TCanvas("c1","c1",400,100,520,400);
    h4->Draw();
    h3->Draw("same");
    h2->Draw("same");
    h1->Draw("same");
    leg->Draw();
    
    TCanvas *c2 =new TCanvas("c2","c2",400,100,520,400);
    h10->Draw();
    h9->Draw("same");
    h8->Draw("same");
    h7->Draw("same");
    leg->Draw();
    
    TCanvas *c3 =new TCanvas("c3","c3",400,100,520,400);
    h16->Draw();
    h15->Draw("same");
    h14->Draw("same");
    h13->Draw("same");
    leg->Draw();
    TCanvas *c4 =new TCanvas("c4","c4",400,100,520,400);
    h22->Draw();
    h21->Draw("same");
    h20->Draw("same");
    h19->Draw("same");
    leg->Draw();
    rootapp->Run();

  return 0;
}//end code
