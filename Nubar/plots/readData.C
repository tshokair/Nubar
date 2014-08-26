{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    //ifstream myFile ("testFreqs.txt");
    //ifstream myFile ("Mode_Map_room_temp_Y14M01D16.txt");
    ifstream myFile ("drdtDataFile.txt");
    //ofstream outFile("Mode_Map_Y14M06D13_fine_1.txt");
    TString x;
    int numEntries=0;
    TString line;
    int numdr=0;
    int numdt=0;
    int headerFlag=0;
    Double_t dr[10000000]={0};
    Double_t dt[100000000]={0};
    int numLines=0;
    TString deltaR;
    TString deltaT;
    int i=0;
    char ch;
    TH2F *p = new TH2F("p","Dr Dt ",500,0,500,100,0,0.1);
    p->SetXTitle("Delta R [cm]");
    p->SetYTitle("Delta T [s]");
    
    while(!myFile.eof()){
		myFile>>x>>deltaR>>deltaT;
        //cerr<<x<<endl;
		numEntries++;
        dt[i]=deltaT.Atof();
        dr[i]=deltaR.Atof();
        if(i%1000==0){
            cerr<<"Event "<<x.Atof()<<endl;
            cerr<<dr[i]<<" "<<dt[i]<<endl;
        }
        p->Fill(dr[i],dt[i]);

        i++;

	}
    cerr<<numEntries<<" lines of data"<<endl;
    p->Draw("colz");
    TLine *lineX = new TLine(200,0,200,0.05);
    lineX->SetLineColor(kRed);
    lineX->SetLineWidth(4);
    lineX->Draw();
    TLine *lineY = new TLine(0,0.05,200,0.05);
    lineY->SetLineColor(kRed);
    lineY->SetLineWidth(4);
    lineY->Draw();
    //cerr<<numEntries<<" entries"<<endl;
	myFile.close();



}
