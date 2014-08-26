#include <TRandom.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TVector3.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TH1.h>
#include <complex>
#include <TNtupleD.h>
#include <TNtuple.h>
#include <TSpline.h>

using namespace std;
void open_the_box_salt_50ms(){

 //Deuglify plots
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptTitle(1);
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetTextFont(132);
  gStyle->SetOptStat(1);
  gStyle->SetMarkerStyle(20);
  gStyle->SetTitleFont(132,"XYZ");
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetOptFit(1111);  // Show fit parameters in plot

  gStyle->SetTitleBorderSize(0);
  gStyle->SetFillColor(0);
  gStyle->SetPadBorderSize(1);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderSize(1);

  // Basic outline of code:  There are two main loops here.  The first loop 
  //scans through the datafile and fills an array showing the time difference 
  //between that event and the one before it.  That loop also fills an array 
  //indicating whether or not each event passes the high level cuts.  The 
  //second main loop applies the high NHIT cut and the NEF cut, which need the 
  //info from the timing array in order to cut events that follow a "bad" 
  //event.  In the second loop, the search window is being applied to look for 
  //bursts.  If there is a burst of multiplicity 2 or more, then the deltaR 
  //cut is applied before the burst is counted or rejected.  If the window 
  //contains a "burst" of multiplicity 1 or more, and the first event of the 
  //burst is not the first event of the window, then the window slides over to 
  //make that event the first event of the window.  Otherwise, the window 
  //moves to the next unexamined event after it has been scanned.

  ifstream goodruns;
  goodruns.open("salt_data.runlist");
  char filename[100];
  
  //flags to define which extra cuts are being turned on; 0=off, 1=on
  Int_t flag_fom=1;
  Int_t flag_beta14=1;
  Int_t flag_nhit=1;
  Int_t flag_rcut=1;

  Int_t multiplicity=0;//Main variable used to count multiplicity in window
  Double_t window_limit=0.05;//window length in seconds

  Float_t Nhits, Itrp, Thetaijp, Angksap, Enekp, Phiksap, Angkspp, Phikspp, Timksp, Qpt_nwp, Qpt_prbp, Qpt_nhit, Enekump, Enekupp, Enerap, Qtijp, Time, Time_us, Time_s, Time_ns, Time_jdy, Date, Beta1p, Beta4p, Rfp, Rdamn1, Rdamn2, Xfp, Yfp, Zfp,Ev_tid;

  TH1D *deltaR2 = new TH1D("deltaR2","Antibox DeltaR mult=2",1200,0.0,1200.);
  deltaR2->SetXTitle("DeltaR (cm)");

  TH1D *deltaR3 = new TH1D("deltaR3","Antibox DeltaR mult=3",1200,0.0,1200.);
  deltaR3->SetXTitle("DeltaR (cm)");

  TH2D *total_drdt=new TH2D("total_drdt","deltaT vs deltaR, mult=2",1200,0.0,1200.0,10005,-5.0,10000.0);
  total_drdt->SetXTitle("deltaR (cm)");
  total_drdt->SetYTitle("deltaT (ms)");
  total_drdt->SetMarkerSize(0.1);

  TH2D *total_drdt3=new TH2D("total_drdt3","deltaT vs deltaR, mult=3",1200,0.0,1200.0,10005,-5.0,10000.0);
  total_drdt3->SetXTitle("deltaR (cm)");
  total_drdt3->SetYTitle("deltaT (ms)");
  total_drdt3->SetMarkerSize(0.1);

  Int_t counts_my_mult2=0;
  Int_t counts_my_mult3=0;
  Int_t counts_my_mult4=0;
  Int_t counts_my_mult5=0;
  Int_t counts_my_mult6=0;

  //parameters for high_nhit cut
  Double_t param0=7.424958;
  Double_t param1=0.027350;
  Double_t param2=-0.001024;
  Double_t param3=0.000007;
  Double_t upper_e_thresh=80.0;
  Double_t nhit_limit=param0+upper_e_thresh*param1+pow(upper_e_thresh,2)*param2+pow(upper_e_thresh,3)*param3;

  cout.precision(8);

  ofstream mylog;
  if(abs(window_limit-0.05)<0.01){mylog.open("box_current_results/box_salt_6.5MeV_0.05s.dat");}
  if(abs(window_limit-1.0)<0.01){mylog.open("box_current_results/box_salt_8.5MeV_1s.dat");}
  if(abs(window_limit-10.)<0.01){mylog.open("box_current_results/box_salt_4.5MeV_10s.dat");}
  //mylog.open("antibox_salt_7.5MeV_0.4s_NEFoff_NHIToff.dat");
  //mylog.open("temp.dat");
  mylog.precision(8);

  mylog<<window_limit<<" "<<flag_fom<<" "<<flag_beta14<<" "<<flag_nhit<<endl;

  //////// Begin scanning through data files
  //for(Int_t m=0;m<1;m++){
  while(!goodruns.eof()){

    goodruns >> filename;
    TFile f(filename,"UPDATE");
    //TFile f("salt_processed_data/sno_024782_p6.root","UPDATE");    

    TTree *t=(TTree*) f.Get("h360");
    //TTree *t2=(TTree*) f.Get("h363");//so far I don't use anything from here
    
    t->SetBranchAddress("Nhits",&Nhits);
    t->SetBranchAddress("Itrp",&Itrp);
    t->SetBranchAddress("Thetaijp",&Thetaijp);
    t->SetBranchAddress("Angksap",&Angksap);
    t->SetBranchAddress("Enekp",&Enekp);
    t->SetBranchAddress("Phiksap",&Phiksap);
    t->SetBranchAddress("Angkspp",&Angkspp);
    t->SetBranchAddress("Phikspp",&Phikspp);
    t->SetBranchAddress("Timksp",&Timksp);
    t->SetBranchAddress("Qpt_nwp",&Qpt_nwp);
    t->SetBranchAddress("Qpt_prbp",&Qpt_prbp);
    t->SetBranchAddress("Qpt_nhit",&Qpt_nhit);
    t->SetBranchAddress("Enekump",&Enekump);
    t->SetBranchAddress("Enekupp",&Enekupp);
    t->SetBranchAddress("Enerap",&Enerap);
    t->SetBranchAddress("Qtijp",&Qtijp);
    t->SetBranchAddress("Time",&Time);    
    t->SetBranchAddress("Time_s",&Time_s);
    t->SetBranchAddress("Time_us",&Time_us);
    t->SetBranchAddress("Time_ns",&Time_ns);
    t->SetBranchAddress("Time_jdy",&Time_jdy);
    t->SetBranchAddress("Date",&Date);
    t->SetBranchAddress("Beta1p",&Beta1p);
    t->SetBranchAddress("Beta4p",&Beta4p);
    t->SetBranchAddress("Rfp",&Rfp);
    t->SetBranchAddress("Rdamn1",&Rdamn1);
    t->SetBranchAddress("Rdamn2",&Rdamn2);
    t->SetBranchAddress("Xfp",&Xfp);
    t->SetBranchAddress("Yfp",&Yfp);
    t->SetBranchAddress("Zfp",&Zfp);
    t->SetBranchAddress("Ev_tid",&Ev_tid);

    Int_t num_entries=t->GetEntries();
    Int_t hlc_pass_array[num_entries];
    Float_t time_array[num_entries+1];
    Double_t skipping_clock[num_entries];
    Float_t skipping_index[num_entries];
    //some cuts look at time of "next event;" last event needs a "next time"
    time_array[num_entries]=999999.0;//bogus entry so I don't fall off stack

    t->GetEvent(0);
    //Time is not an accurate variable to use; use Time_s, Time_us, ect.
    Float_t day_prev=Time_jdy;
    Float_t hour_prev=Int_t(Time_s)/60/60;
    Float_t minute_prev=Int_t(Time_s - hour_prev*60*60)/60;
    Float_t second_prev=Int_t(Time_s - hour_prev*60*60 - minute_prev*60);
    Float_t usecond_prev=Time_us;
    Float_t nsecond_prev=Time_ns;
    Float_t time_in_sec_prev=Time_s;

    cout<<"Run "<<filename<<endl;
    cout<<"Number of Events "<<num_entries<<endl;
    mylog<<"Run "<<filename<<endl;
    mylog<<"Number of Events "<<num_entries<<endl;
    cout<<"Filling time and HLC arrays"<<endl;
    if(0xB7 != 183 || 0x6FE1 != 28641){cout<<"ERROR: wrong damn masks!"<<endl;}

    ///// Begin scanning through events:
    ///// Count number of bursts; fill the time array; fill the HLC array
    for(Int_t n=0;n<num_entries;n++){
      t->GetEvent(n);
      //t2->GetEvent(n);
          
      //Fill time array with the time difference between events in seconds
      //Breaking the time down this way helps when the day changes
      Float_t day=Time_jdy;
      Float_t hour=Int_t(Time_s)/60/60;
      Float_t minute=Int_t(Time_s - hour*60*60)/60;
      Float_t second=Int_t(Time_s - hour*60*60 - minute*60);
      Float_t usecond=Time_us;
      Float_t nsecond=Time_ns;

      Float_t time_diff=0.0;//time diff in seconds
      Float_t day_diff=(day-day_prev);
      Float_t hour_diff=(hour-hour_prev);
      Float_t minute_diff=(minute-minute_prev);
      Float_t second_diff=second+usecond/1.E6+nsecond/1.E9-(second_prev+usecond_prev/1.E6+nsecond_prev/1.E9);
      time_diff = day_diff*24*60*60+hour_diff*60*60+minute_diff*60+second_diff;
      //time_dif_s is simpler and works for all events except day-change events
      Float_t time_dif_s=Time_s+usecond/1.E6+nsecond/1.E9-(time_in_sec_prev+usecond_prev/1.E6-nsecond_prev/1.E9);
      if(day==day_prev){time_array[n]=time_dif_s;}
      else{time_array[n]=time_diff;}

      if(fabs(time_dif_s-time_diff)>0.0001 && day_diff==0){mylog<<"Time ERROR 1: "<<n<<" "<<time_dif_s<<" "<<time_diff<<" "<<day_diff<<endl;}
      if(time_diff<0.0||day_diff<0.0){mylog<<"Time ERROR 2: "<<n<<" "<<time_diff<<" "<<day_diff<<endl;}

      day_prev=day;
      hour_prev=hour;
      minute_prev=minute;
      second_prev=second;
      usecond_prev=usecond;   
      nsecond_prev=nsecond;
      time_in_sec_prev=Time_s; 

      hlc_pass_array[n]=0;
      Int_t itrp_flag=0;
      Int_t thetaijp_flag=0;
      Int_t angksap_flag=0;
      Int_t angkspp_flag=0;
      Int_t timksp_flag=0;
      Int_t qpt_flag=0;
      Int_t eneku_flag=0;
      Int_t other_flag=0;
      Int_t qtijp_flag=0;
      Int_t damn_flag=0;
      Int_t rfp_flag=0;
      Int_t energy_flag=0;
      if( (0.74-Itrp)<2.7*0.43/sqrt(Qpt_nhit) ){itrp_flag=1;}
      if( Thetaijp>0.89 && Thetaijp<1.6 ){thetaijp_flag=1;}
      if( Angksap>=0.902/pow(Enekp,4) || Phiksap>=0.0422 ){angksap_flag=1;}
      if( Angkspp>=16.71/pow(Enekp,4) || Phikspp>=0.0066 ){angkspp_flag=1;}
      if( (Timksp>0.0001 && Timksp<=1.0)||((Timksp-2.0)>0.0001 && Timksp<=3.0) ){timksp_flag=1;}    
      if( (TMath::Poisson(Qpt_nwp,Nhits*0.0062)>0.002 && (1-pow(1-Qpt_prbp,Qpt_nwp))>0.01) || Qpt_nwp==0 ){qpt_flag=1;}
      if( Enekp>13 || (Enekump/(-0.13288+0.3904*sqrt(Enekp-0.511)+0.0251*(Enekp-0.511))>(-0.867-0.0239*Enekp+0.0009*pow(Enekp,2)) && Enekupp/(-0.13288+0.3904*sqrt(Enekp-0.511)+0.0251*(Enekp-0.511))<(1.13+0.0164*Enekp-0.0009*pow(Enekp,2))) ){eneku_flag=1;}
      if( (Enerap/Enekp)/Itrp < 1.5 || Qtijp>0.65){other_flag=1;}
      if(Qtijp>0.4){qtijp_flag=1;}
      if( (Int_t(Rdamn1) & 0x04B7)==0 && (Int_t(Rdamn2) & 0x6FE1)==0 ){damn_flag=1;}
      if(Rfp<550.0){rfp_flag=1;}
      if(Enekp>8.5 && abs(window_limit-1.0)<0.01){energy_flag=1;}//energy threshold corresponds to window
      else if(Enekp>6.5 && abs(window_limit-0.05)<0.01){energy_flag=1;}//energy threshold corresponds to window
      else if(Enekp>4.5 && abs(window_limit-10.)<0.01){energy_flag=1;}//energy threshold corresponds to window
      //Apply all HLCs for events below 15MeV
      if(itrp_flag==1 && thetaijp_flag==1 && rfp_flag==1 && Enekp>15.0 && damn_flag==1 && energy_flag==1){hlc_pass_array[n]=1;}
      else if(itrp_flag==1 && thetaijp_flag==1 && angksap_flag==1 && angkspp_flag==1 && timksp_flag==1 && qpt_flag==1 && eneku_flag==1 && qtijp_flag==1 && other_flag==1 && rfp_flag==1 && damn_flag==1 && energy_flag==1){hlc_pass_array[n]=1;}//assign each event a pass or fail for HLCs

    }//ends for loop that counts bursts, fills time array, & fills HLC array

    cout<<"Applying high NHIT, NEF, and DeltaR cuts"<<endl;
    //Apply High NHIT Cut and NEF Cut, along with deltaR
    Int_t max_burst=100;
    Int_t burst_exceeding_max=0;
    Float_t radius_array[max_burst];
    Float_t x_array[max_burst];
    Float_t y_array[max_burst];
    Float_t z_array[max_burst];
    Float_t event_energy[max_burst];
    Float_t event_time_diff[max_burst];
    Float_t event_num[max_burst];
    Float_t event_evtid[max_burst];
    Int_t radius_array_index=0;
    Int_t pass_test_beta14=0;
    Int_t pass_test_fom=0;
    Int_t pass_test_nhit=0;
    Int_t new_burst_flag=1;//tells window to start with this event
    Int_t skipping_some_events=0;
    Double_t time_clock=0.0;
    Double_t rcut=0.0;
    Double_t tcutoff=0.0;

    Double_t alpha=1.0;//parameter for drdt cut, sqrt model; mult2
    Double_t alpha2=45.0;//parameter for drdt cut, box model; mult2
    Double_t beta=40.;//mult999
    Double_t beta2=5.0;//mult999
    Double_t beta3=200.0;//mult999

    if(abs(window_limit-1.0)<0.01){rcut=410.0;tcutoff=1000.;}
    else if(abs(window_limit-0.05)<0.01){rcut=385.0;tcutoff=50.0;}
    else if(abs(window_limit-10.)<0.01){rcut=400.0;tcutoff=10000.;}
    else{cout<<"ERROR with window_limit and rcut"<<endl;}
    Double_t rtest=0.0;
    Double_t window=0.0;
    Int_t first_burst_event=0;//stores event number of first event in burst
    Int_t first_window_event=0;//event number of where the window starts
    for(Int_t n=0;n<num_entries;n++){
      t->GetEvent(n);
      window=0.0;
      if(new_burst_flag==1){
        radius_array_index=0;
        //Do not reset these 2 in case an atm burst overlaps 2 windows
        //skipping_some_events=0;//flag for followers within 200ms
        //time_clock=0.0;
        new_burst_flag=0;
        multiplicity=0;
        first_window_event=n;
        while(window<window_limit){
	  pass_test_beta14=0;
          pass_test_fom=0;
          pass_test_nhit=0;
          //Which tests are passed?
          Double_t beta_14=Beta1p + 4.0*Beta4p;
          if(Nhits<25){pass_test_fom=1;}//nhit too low to trigger a fail
          else if((Angksap>0.902/pow(Enekp,4) || Phiksap>0.0422) && (Angkspp>16.71/pow(Enekp,4) || Phikspp>0.0066) ){pass_test_fom=1;}
          else{if(flag_fom==0){pass_test_fom=1;}}//ignore flag if it's off
          if(Nhits<25){pass_test_beta14=1;}//nhit too low to trigger a fail
          else if(beta_14>-0.12 && beta_14<0.95){pass_test_beta14=1;}
          else{if(flag_beta14==0){pass_test_beta14=1;}}//ignore flag if off
          if(Nhits<nhit_limit*upper_e_thresh){pass_test_nhit=1;}
          else{if(flag_nhit==0){pass_test_nhit=1;}}//ignore flag if it's off
          if(pass_test_nhit==0 || skipping_some_events==1.0 || pass_test_fom==0 || pass_test_beta14==0){
            //if(skipping_some_events==1.0){time_clock += time_array[n+1];}
            //else{time_clock=time_array[n+1]; skipping_some_events=1.0;}
            skipping_some_events=1;
            if(pass_test_nhit==0 || pass_test_fom==0 || pass_test_beta14==0){time_clock=0.0;}
            else{time_clock += time_array[n];}
            if(time_clock>0.2){//200msec
              skipping_some_events=0;
              time_clock=0.0;
            }
          }
	  skipping_clock[n]=time_clock;//record in case window moves
          skipping_index[n]=skipping_some_events;

          if(hlc_pass_array[n]==1 && pass_test_beta14==1 && pass_test_fom==1 && pass_test_nhit==1 && skipping_some_events==0){
            multiplicity++;
            if(radius_array_index<max_burst){
              radius_array[radius_array_index]=Rfp;
              x_array[radius_array_index]=Xfp;
              y_array[radius_array_index]=Yfp;
              z_array[radius_array_index]=Zfp;
              event_energy[radius_array_index]=Enekp;
              event_time_diff[radius_array_index]=window;
              //event_time_diff[radius_array_index]=time_array[n];
              event_num[radius_array_index]=n;
              event_evtid[radius_array_index]=Ev_tid;
            }
            radius_array_index++;
            if(radius_array_index>max_burst){mylog<<"WARNING: mult>"<<max_burst<<"; Mult & total overspill bursts "<<radius_array_index<<" "<<burst_exceeding_max<<endl;}
            if(multiplicity==1){first_burst_event=n;}
          }
          window += time_array[n+1];        
          if(window<window_limit){n++; t->GetEvent(n);} //t2->GetEvent(n);}
          else{new_burst_flag=1;}
        }//end while window<0.4 loop
      }//end "if new burst" statement

      //Apply rcut
      Double_t rcut_sum=0.0;
      Double_t mC2=0;//number of non-redundant combinations
      for(Int_t j=0;j<multiplicity-1;j++){
        for(Int_t k=j+1;k<multiplicity;k++){
          rcut_sum += sqrt( pow(x_array[j]-x_array[k],2) + pow(y_array[j]-y_array[k],2) + pow(z_array[j]-z_array[k],2) );
          if(j<k){mC2++;}
        }
      }
      if(mC2>0.0){rtest=rcut_sum/mC2;} else{rtest=0.0;}
      if(multiplicity>max_burst){rtest=0.0;}//fails rtest otherwise
      //if(multiplicity==2 && rtest>rcut){cout<<radius_array[0]<<" "<<radius_array[1]<<" "<<mC2<<endl;}
      //cout<<n<<" "<<rtest<<" "<<rcut<<endl;

      //If window had one event, start window at that event and rescan
      if(multiplicity>=1 && first_burst_event != first_window_event && new_burst_flag==1){n=first_burst_event-1;time_clock=skipping_clock[first_burst_event-1];skipping_some_events=skipping_index[first_burst_event-1];}
      else if(flag_rcut==1){
        Double_t deltat=event_time_diff[1]*1000.;
        Int_t flag_drdt=0;
        //if(deltat<tcutoff && (rtest>rcut || deltat>alpha2)){flag_drdt=1;}//box
        //if(deltat<tcutoff && (rtest<rcut && deltat<alpha2)){flag_drdt=1;}//antibox
        if(deltat>beta || rtest>rcut || (rtest>beta3 && deltat>sqrt(beta*beta+beta2*beta3-rtest*beta2))){flag_drdt=1;}//box
	//if(deltat<beta && rtest<rcut && (rtest<beta3 && deltat<sqrt(beta*beta+beta2*beta3-rtest*beta2))){flag_drdt=1;}//antibox
        if(multiplicity>=2 && flag_drdt==1){
          counts_my_mult2++;
          if(multiplicity==2){deltaR2->Fill(rtest);total_drdt->Fill(rtest,deltat);}
          mylog<<"Burst event num "<<event_num[0]<<" "<<event_num[1]<<" "<<event_num[2]<<" "<<event_num[3]<<" "<<event_num[5]<<endl;
          mylog<<"Burst ev_tid "<<event_evtid[0]<<" "<<event_evtid[1]<<" "<<event_evtid[2]<<" "<<event_evtid[3]<<" "<<event_evtid[5]<<endl;
          mylog<<"Burst event energy "<<event_energy[0]<<" "<<event_energy[1]<<" "<<event_energy[2]<<" "<<event_energy[3]<<" "<<event_energy[5]<<endl;
          mylog<<"Burst event time diff "<<event_time_diff[0]<<" "<<event_time_diff[1]<<" "<<event_time_diff[2]<<" "<<event_time_diff[3]<<" "<<event_time_diff[5]<<endl;
          mylog<<"Burst event deltaR "<<rtest<<endl;
        }
        if(multiplicity==3 && flag_drdt==1){
          Double_t deltat=(event_time_diff[1]+event_time_diff[2]+event_time_diff[1]+event_time_diff[2])/3.;
          total_drdt->Fill(rtest,deltat);
        }
        if(multiplicity>=3 && flag_drdt==1){
          counts_my_mult3++;
          if(multiplicity==3){deltaR3->Fill(rtest);}
        }
        if(multiplicity>=4 && flag_drdt==1){counts_my_mult4++;}
        if(multiplicity>=5 && flag_drdt==1){counts_my_mult5++;}
        if(multiplicity>=6 && flag_drdt==1){counts_my_mult6++;;if(multiplicity>max_burst){burst_exceeding_max++;}}
      }
      else{
        if(multiplicity>=2){counts_my_mult2++;}
        if(multiplicity>=3){counts_my_mult3++;}
        if(multiplicity>=4){counts_my_mult4++;}
        if(multiplicity>=5){counts_my_mult5++;}
        if(multiplicity>=6){counts_my_mult6++;}
      }

      for(Int_t j=0;j<max_burst;j++){radius_array[j]=0.0;x_array[j]=0.0;y_array[j]=0.0;z_array[j]=0.0;event_num[j]=0.0;event_energy[j]=0.0;event_time_diff[j]=0.0;event_evtid[j]=0.0;}//reset position arrays

    }// end loop through entries

    //TFile f1(rootfilename,"RECREATE");
    TFile f1("box_current_results/box_salt_6.5MeV_0.05s.root","RECREATE");
    //TFile f1("antibox_current_results/antibox_salt_8.5MeV_1s.root","RECREATE");
    //TFile f1("antibox_current_results/antibox_salt_4.5MeV_10s.root","RECREATE");
    //TFile f1("temp.root","RECREATE");
    deltaR2->Write();
    deltaR3->Write();
    total_drdt->Write();
    total_drdt3->Write();

    mylog<<counts_my_mult2<<" "<<counts_my_mult3<<" "<<counts_my_mult4<<" "<<counts_my_mult5<<" "<<counts_my_mult6<<endl;

    t->Delete();
    //t2->Delete();


  }//closes big loop that reads datafiles

  cout<<counts_my_mult2<<" "<<counts_my_mult3<<" "<<counts_my_mult4<<" "<<counts_my_mult5<<" "<<counts_my_mult6<<endl;

}//ends program
