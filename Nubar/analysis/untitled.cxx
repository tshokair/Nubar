/*
 *  untitled.h
 *  
 *
 *  Created by Tim Shokair on 5/18/11.
 *  Copyright 2011 University of Pennsylvania. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TRoot.h>
#include <TRint.h>
#include <TMath.h>
#include <TTree.h>
#include <TH1F.h>
#include <TFile.h>

using namespace std;

int main(int argc, char** agrv)
{
	Float_t nhits;
	Float_t enep;
	Float_t eneg;
	Float_t time;
	Float_t jdy;
	Float_t ut1;
	Float_t ut2;
	Float_t jdy_i;
	Float_t ut1_i;
	Float_t ut2_i;
	Int_t numev;
	Int_t nump;
	Int_t mcNumP;
	Int_t numg1;
	Int_t numg2;
	Int_t mcNumG;
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
	Float_t rmag;
	Float_t r0[10000];
	Float_t r1[10000];
	Float_t r2[10000];
	Float_t dt1;
	Float_t dt2;
	Float_t temp_pos;
	Float_t ediff;
	Float_t thetaij;
	Int_t numfiles;
	TString files[1000];
	Int_t nd;
	Int_t days;
	Float_t day_diff;
	Int_t eventID[10000][3];
	Float_t evid;
	Int_t partNum;
	Int_t tempEventID;
	Float_t idPos;
	Float_t idNeutron;
	Int_t isFirst;
	Int_t isSecond;
	Int_t isThird;
	Float_t tempTime;
	Int_t numCoincidence;	
	
	//files[0]="../ntuples/anticc/run_nuebar_b8_cc_rsp15949_combined.root";
	//files[1]="../ntuples/anticc/run_nuebar_b8_cc_rsp15958_combined.root";
	files[0]="run_nuebar_b8_cc_rsp15978_combined.root";
	// files[1]="../ntuples/anticc/run_nuebar_b8_cc_rsp15997_combined.root";
	//files[4]="../ntuples/anticc/run_nuebar_b8_cc_rsp15998_combined.root";
	//files[5]="../ntuples/anticc/run_nuebar_b8_cc_rsp16002_combined.root";
	//files[6]="../ntuples/anticc/run_nuebar_b8_cc_rsp16003_combined.root";
	//files[7]="../ntuples/anticc/run_nuebar_b8_cc_rsp16013_combined.root";
	
	numfiles =1;
	
	TH1F * h1= new TH1F("h1","Nhits",100,0,200);
	TH1F * h2= new TH1F("h2","Positron Energy",100,0,20);
	TH1F * h3= new TH1F("h3","First Gamma Energy",100,0,20);
	TH1F * h5= new TH1F("h5","Second Gamma Energy",100,0,20); 
	TH1F * h4= new TH1F("h4","Event Time",100000,0,90000);
	TH1F * h6= new TH1F("h6","Thetaij positrons",100,0,3);
	TH1F * h7= new TH1F("h7","Thetaij gammas",100,0,3);
	TH1F * h8= new TH1F("h8","Nhit positrons",300,0,150);
	TH1F * h9= new TH1F("h9","Nhit gammas",150,0,150);
	TH1F * h10= new TH1F("h10","Event Energy",150,0,20);
	TH1F * h11= new TH1F("h11","Event Day",1000,0,1);
	TH1F * h12= new TH1F("h12","Time Difference seconds",100000,0,10);
	TH1F * h13= new TH1F("h13","Time Difference nanoseconds",1000,0,1000000000);
	TH1F * h14= new TH1F("h14","Vertex position cm",100000,0,600);
	h2->SetLineColor(1);
	h3->SetLineColor(2);
	h5->SetLineColor(3);
	h6->SetLineColor(3);
	h7->SetLineColor(2);
	h8->SetLineColor(7);
	h9->SetLineColor(2);
	//cerr<<files[0]<<endl;
	
	numev=0;
	nump=0;
	mcNumP=0;
	mcNumG=0;
	numg1=0;
	numg2=0;
	days=0;
	numCoincidence=0;
	
	
	for (Int_t k=0; k<numfiles;++k)
    {
		xi=0;
		yi=0;
		zi=0;
		ri=0;
		is_pos=0;
		is_gamma=0;
		isFirst=0;
		isSecond=0;
		isThird=0;
		tempTime=0;
		
		//cout<<files[k]<<endl;
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("Combined");
		nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		Td->SetBranchAddress("h149_Nhits",&nhits);
		Td->SetBranchAddress("h149_Enep",&enep);
		Td->SetBranchAddress("h149_Eneg",&eneg);
		Td->SetBranchAddress("h149_tid",&evid);	
		Td->SetBranchAddress("h149_Ev_jdy",&jdy); //julian day
		Td->SetBranchAddress("h149_Ev_ut1",&ut1); //time in seconds
		Td->SetBranchAddress("h149_Ev_ut2",&ut2); //time in nanoseconds
		Td->SetBranchAddress("h149_Xf",&xfit);
		Td->SetBranchAddress("h149_Yf",&yfit);
		Td->SetBranchAddress("h149_Zf",&zfit);
		Td->SetBranchAddress("h149_Thetaij",&thetaij);
		Td->SetBranchAddress("h147_Id_in",&idPos); //Positron Truth information
		Td->SetBranchAddress("h148_Id_in",&idNeutron); //Neutron Truth information
		Td->GetEntry(0);
		jdy_i=jdy;
		ut1_i=ut1;
		ut2_i=ut2;
		
		for(Int_t j=0; j<nd; ++j){
			Td->GetEntry(j);
			rfit=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi));
			rmag=TMath::Sqrt(xfit*xfit+yfit*yfit+zfit*zfit);
			//get the day difference
			day_diff = jdy-jdy_i;
			//h11->Fill(day_diff);
			//FIXME: take care of events that span two days
			if(day_diff>=1 && day_diff<2){
				dt1=(86400-ut1_i)+ut1;
				if(dt1<1){
					dt2=ut2-ut2_i;
				}
				else if(dt1==1){
					dt2=(1000000000-ut2_i)+ut2;
				}
				cerr<<day_diff<<" days "<<dt1<<" s  "<<dt2<<" ns between events"<<endl;
			}
			else if(day_diff<1 &&day_diff>=0){
				dt1=ut1-ut1_i;
				if(dt1<1){
					dt2=ut2-ut2_i;
				}
				else if(dt1>1 &&dt1<2){
					dt2=(1000000000-ut2_i)+ut2;
					//cerr<<dt2<<endl;
				}
			}
			
			if (nhits>10 &&nhits<150 &&rmag<550&&day_diff<=1){
				numev++;
				h1->Fill(nhits);
				h10->Fill(enep);
				//dt1=ut1-ut1_i;
				//dt2=ut2-ut2_i;
				h4->Fill(ut1);
				h12->Fill(dt1);
				h13->Fill(dt2);
				//FIXME find out a real way to tag positrons
				if(idPos==21){
					mcNumP++;
					is_pos=1;
				}
				else{
					mcNumG++;
					is_pos=0;
				}		
				/*eventID[numCoincidence][0]=0;
				 eventID[numCoincidence][1]=0;
				 eventID[numCoincidence][2]=0;
				 */
				if((isFirst==0 && isSecond==0)||isThird==1||rfit>350||dt1>1||(dt1<1 &&dt2>150000000)){
					tempEventID=evid;
					temp_pos=enep;
					ut1_i=ut1;
					ut2_i=ut2;
					xi=xfit;
					yi=yfit;
					zi=zfit;
					r0[numev]=rfit;
					isFirst=1;
					isSecond=0;
					isThird =0;
					is_pos=1;
					nump++;
					//cerr<<"found "<<nump<<" positrons and this is the first particle "<<isFirst<<endl;
					is_gamma=0;
					h6->Fill(thetaij);
					h14->Fill(rmag);
					h8->Fill(nhits);
					
					jdy_i=jdy;
					
					
				}
				else if ( /*is_pos==0 &&*/ rfit<350 &&dt1<=1){
					if(isFirst==1 && isSecond==0 && dt2<150000000 && dt2>0 ){
						eventID[numCoincidence][0]=tempEventID;
						//store the event ID of the initial positron (This id is only for events in one file
						eventID[numCoincidence][1]=evid;
						numCoincidence++;
						//store the event ID of the first coincidence particle ie first gamma detected
						h2->Fill(temp_pos); 
						//only store the enegy for a "positron" if it is followed by a neutron
						h3->Fill(enep); 
						//store the energy of the neutron/gamma
						isSecond=1;
						isFirst=0;
						/*isThird=0;
						 ut1_i=ut1;*/
						tempTime=ut2;
						numg1++;
						//cerr<<"found "<<numg1<<" first neutrons  and this is the second particle "<<isSecond<<endl;
						h7->Fill(thetaij);
						//h4->Fill(ut2);
						//cerr<<dt1<<endl;
						h9->Fill(nhits);
						jdy_i=jdy;
						cerr<<"neutron is "<<idNeutron<<" and is neutron 1 detected "<<dt2<<" ns after the positron "<<endl;
					}
					else if (isSecond==1 && dt2<150000000 && dt2>0 &&tempTime<ut2 ){
						eventID[numCoincidence][2]=evid;//store event ID of second coincidence particle
						h5->Fill(enep);
						r2[numev]=rfit;
						isFirst=0;
						isSecond=0;
						isThird=1;
						numg2++;
						//h4->Fill(ut2);
						h9->Fill(nhits);
						//cerr<<"event occured on day "<<jdy<<" "<<day_diff<<" days after previous "<<endl;
						jdy_i=jdy;
						cerr<<"neutron is "<<idNeutron<<" and is neutron 2 detected "<<dt2<<" ns after the positron"<<endl;
					}
					else{
						isFirst=0;
						isSecond=0;
						isThird=0;
						cerr<<"Nothing fits parameters and code is resetting"<<endl;
					}
					
				}
				else{
					isFirst=0;
					isSecond=0;
				}
			}
			else jdy_i=jdy;
		}
    }
	h8->Draw();
	h9->Draw("same");
	cerr<<numCoincidence<<" candidate antinuetrinos"<<endl;
	for (Int_t n=0; n<=numCoincidence; n++) {
		if(eventID[n][2]==0){
			cerr<<"2Fold event with event ids "<<eventID[n][0]<<" "<<eventID[n][1]<<" "<<endl;
		}
		if(eventID[n][2]!=0){
			cerr<<"3Fold event with event ids "<<eventID[n][0]<<" "<<eventID[n][1]<<" "<<eventID[n][2]<<endl;
		}
	}
	cerr<<numev<<" events"<<endl;
	cerr<<nump<<" positrons "<<mcNumP<<" mc positrons" <<endl;
	cerr<<numg1<<" first gammas "<<mcNumG<<" from mc neutrons" <<endl;
	cerr<<numg2<<" second gammas"<<endl;
	
	return 0;
	
}

