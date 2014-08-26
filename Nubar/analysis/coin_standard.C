//#include <string>
//#include <iostream>
//#include <fstream>
{
	Float_t nhits;
	Float_t enep;
	Float_t eneg;
	Float_t time;
	Float_t jdy;
	Float_t ut1;
	Float_t ut2;
	Float_t ut3;
	Float_t jdy_i;
	Float_t ut1_i;
	Float_t ut2_i;
	Float_t ut3_i;	
	Float_t itr;
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
	Float_t rf;
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
	Double_t tempTime;
	Int_t numCoincidence;
	Double_t totalTime;	
	Double_t totalTime_i;
	Double_t timediff;
	Float_t rdamn1;
	Float_t rdamn2;
	Int_t damnflag;
	Int_t numRawEvents;
	
	//parameters used in cuts
	Float_t nhitCut=27.5;
	Float_t nhitMax=150;
	Float_t fiducialVolume=550;
	Float_t timeWindow=0.150;
	Float_t distanceWindow=500;
	Float_t thetamin=0.75;
	Float_t thetamax=1.45;
	Float_t itrmin= 0.55;
	Float_t itrmax=0.95;
	
	ifstream runList;
	runList.open("runlist.txt");
	numfiles =0;
	while(!runList.eof()){
		runList>>files[numfiles];
		numfiles++;
	}
	runList.close();
	
	
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
	TH1F * h12= new TH1F("h12","Time Difference seconds",100000,0,100);
	TH1F * h13= new TH1F("h13","Time Difference nanoseconds",1000,0,1000000000);
	TH1F * h14= new TH1F("h14","Vertex position cm",100000,0,600);
	h2->SetLineColor(1);
	h3->SetLineColor(2);
	h5->SetLineColor(3);
	h6->SetLineColor(3);
	h7->SetLineColor(2);
	h8->SetLineColor(7);
	h9->SetLineColor(2);
	
	numev=0;
	nump=0;
	mcNumP=0;
	mcNumG=0;
	numg1=0;
	numg2=0;
	days=0;
	numCoincidence=0;
	numRawEvents=0;
	
	
	for (Int_t k=0; k<(numfiles-1);++k)//loop over ntuples/root files
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
		
		TFile *f1 = new TFile(files[k]);
		TTree *Td = (TTree*)f1->Get("h330");
		nd = Td->GetEntries();
		cerr<<nd<<" events in file "<<k<<endl;
		
		Td->SetBranchAddress("Nhits",&nhits);
		Td->SetBranchAddress("Enep",&enep);
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
		Td->GetEntry(0);
		//initialize times
		jdy_i=jdy;
		ut1_i=Int_t(ut1);
		ut2_i=Int_t(ut2);
		ut3_i=Int_t(ut3);
		for(Int_t j=0; j<nd; ++j){ //start looping over events in files
			Td->GetEntry(j);
			rfit=TMath::Sqrt((xfit-xi)*(xfit-xi)+(yfit-yi)*(yfit-yi)+(zfit-zi)*(zfit-zi)); //get distance between verticies of current and previous event
			rmag=TMath::Sqrt(xfit*xfit+yfit*yfit+zfit*zfit); //get distance from center of detector to check if event is in fiducial volume
			
			day_diff = jdy-jdy_i;//get the day difference
			totalTime_i=Int_t(ut1_i)+Int_t(ut2_i)/1.E6-Int_t(ut3_i)/1.E9;
			if((Int_t(rdamn1)& 0x0b32)==0 && (Int_t(rdamn2)&& 0x1ec1)==0){
				damnflag =1;  
			}
			else {
				damnflag=0;
			}
			

			if(Int_t(day_diff)==1){
				dt1=(86400-ut1_i)+ut1;
				if(Int_t(dt1)==0){
					dt2=ut2-ut2_i;
					dt3=ut3-ut3_i;
					timediff=dt1+dt2/1.E6 +dt3/1.E9;
					
				}
				else if(Int_t(dt1)==1){
					dt2=(1-ut2_i/1.E6)+ut2/1.E6;
					dt3=(1-ut3_i/1.E9)+ut3/1.E9;
					timediff=dt2+dt3;
				}
				else timediff=10;
				
				//cerr<<day_diff<<" days "<<dt1<<" s  "<<dt2<<" ns between events"<<endl;
			}
			else if(Int_t(day_diff)==0){
				dt1=ut1-ut1_i;
				if(Int_t(dt1)==0){
					dt2=ut2-ut2_i;
					dt3=ut3-ut3_i;
					timediff=dt1+dt2/1.E6+dt3/1.E9;
				}
				else if(Int_t(dt1)==1){
					dt2=(1-ut2_i/1.E6)+ut2/1.E6;
					dt3=(1-ut3_i/1.E9)+ut3/1.E9;
					timediff=dt2+dt3;
				}
				else timediff =10;
			}
			numRawEvents++;
			if (nhits>nhitCut &&nhits<nhitMax &&rmag<fiducialVolume && thetaij>thetamin &&thetaij<thetamax && itr<itrmax&&itr>itrmin &&damnflag==1){
				numev++;
				h1->Fill(nhits);
				h10->Fill(enep);
				totalTime=Int_t(ut1)+Int_t(ut2)/1.E6-Int_t(ut3)/1.E9;


				h4->Fill(ut1);
				h12->Fill(dt1);
				h13->Fill(dt2);

				/*
				 eventID[numCoincidence][0]=0;
				 eventID[numCoincidence][1]=0;
				 eventID[numCoincidence][2]=0;
				 */
				if((isFirst==0 && isSecond==0)||isThird==1||rfit>distanceWindow||timediff>timeWindow){
					tempEventID=evid;
					temp_pos=enep;
					ut1_i=ut1;
					ut2_i=ut2;
					ut3_i=ut3;
					xi=xfit;
					yi=yfit;
					zi=zfit;
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
					tempTime=totalTime;
					jdy_i=jdy;
					
					
				}
				else if (rfit<distanceWindow &&timediff<timeWindow){
					if(isFirst==1 && isSecond==0 &&timediff>0 ){
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
						//restart clock
						/*
						ut1_i=ut1;
						ut2_i=ut2;
						ut3_i=ut3;
						*/
						//cerr<<tempTime<<endl;
						//cerr<<totalTime<<endl;
						numg1++;
						//cerr<<"found "<<numg1<<" first neutrons  and this is the second particle "<<isSecond<<endl;
						h7->Fill(thetaij);
						//h4->Fill(ut2);
						//cerr<<dt1<<endl;
						h9->Fill(nhits);
						cerr<<timediff<<" ms time between initial and first at "<<rfit<<endl;
						cerr<<itr<<" itr "<<thetaij<<" thetaij "<< nhits<<" nhit"<<endl;
						
						jdy_i=jdy; 
						//cerr<<"neutron is "<<idNeutron<<" and is neutron 1 detected "<<dt2<<" ns after the positron "<<endl;
					}
					else if (isSecond==1 && timediff>0 ){
						eventID[numCoincidence][2]=evid;//store event ID of second coincidence particle
						h5->Fill(enep);
						isFirst=0;
						isSecond=0;
						isThird=1;
						xi=xfit;
						yi=yfit;
						zi=zfit;
						numg2++;
						//cerr<<totalTime<<endl;
						//h4->Fill(ut2);
						h9->Fill(nhits);
						cerr<<timediff<<" ms time between initial and second "<<endl;
						//cerr<<"event occured on day "<<jdy<<" "<<day_diff<<" days after previous "<<endl;
						jdy_i=jdy;
						//cerr<<"neutron is "<<idNeutron<<" and is neutron 2 detected "<<dt2<<" ns after the positron"<<endl;
					}
					/*else{
						isFirst=0;
						isSecond=0;
						isThird=0;
						xi=xfit;
						yi=yfit;
						zi=zfit;
						cerr<<"Nothing fits parameters and code is resetting"<<endl;
					}
					 */
					
				}
				/*else{
					//isFirst=0;
					//isSecond=0;
					//xi=xfit;
					//yi=yfit;
					//zi=zfit;
					jdy_i=jdy;
				}*/
			}
			else { 
				if (timediff>=.150) {
					isFirst=0;
				}
				jdy_i=jdy;

				
			}



		}
		}//end loop over ntuples/root files
		h8->Draw();
		h9->Draw("same");
		cerr<<numCoincidence<<" candidate antinuetrinos"<<endl;
		for (Int_t n=0; n<numCoincidence; n++) {
		 if(eventID[n][2]==0){
		 cerr<<"2Fold event with event ids "<<eventID[n][0]<<" "<<eventID[n][1]<<" "<<endl;
		 }
		 if(eventID[n][2]!=0){
		 cerr<<"3Fold event with event ids "<<eventID[n][0]<<" "<<eventID[n][1]<<" "<<eventID[n][2]<<endl;
		 }
		 }
		 
		//cerr<<numRawEvents<<" raw events"<<endl;
		cerr<<numev<<" events"<<endl;
		cerr<<nump<<" positrons "<<mcNumP<<" mc positrons" <<endl;
		cerr<<numg1<<" first gammas "<<mcNumG<<" from mc neutrons" <<endl;
		cerr<<numg2<<" second gammas"<<endl;
		
		}
