{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    //ifstream myFile ("testFreqs.txt");
    //ifstream myFile ("Mode_Map_room_temp_Y14M01D16.txt");
    ifstream myFile ("accDist.txt");
    //ofstream outFile("Mode_Map_Y14M06D13_fine_1.txt");
    TString x;
    int numEntries=0;
    TString line;
    int numdr=0;
    int numdt=0;
    int headerFlag=0;
    Double_t dr[10000000]={0};
    Double_t dt[100000000]={0};
    Double_t wt[100000000]={0};
    int numLines=0;
    TString deltaR;
    TString deltaT;
    TString weight;
    int i=0;
    char ch;
    TH2F *p = new TH2F("p","Dr Dt ",650,0,1300,150,0,0.150);
    p->SetXTitle("Delta R [cm]");
    p->SetYTitle("Delta T [s]");
    TRandom3 r; // generate a number in interval ]0,1] (0 is excluded)
    r.Rndm();
    
    while(!myFile.eof()){
        
		myFile>>deltaR>>weight;
        //cerr<<x<<endl;
		numEntries++;
        wt[i]=weight.Atof();
        dr[i]=deltaR.Atof();
        cerr<<dr[i]<<" "<<wt[i]<<endl;
        i++;

	}
    for(int j=0;j<i;j++){
        for(int k=0;k<wt[j];k++){
            p->Fill(dr[j],r.Rndm()*.150);
        }
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
