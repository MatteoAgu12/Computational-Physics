\\ Resistivity as function of temperature for metals and semiconductorss
\\ 2-wires and 4-wires measurements techniques


int RT()
{
    TFile *time = new TFile("RT_time.txt", "read");
    
    TFile *pot1 = new TFile("RT_potential1.txt", "read");
    TFile *pot2 = new TFile("RT_potential2.txt", "read");
    TFile *pot3 = new TFile("RT_potential3.txt", "read");
    TFile *pot4 = new TFile("RT_potential4.txt", "read");

    TFile *R1 = new TFile("RT_relative_Ge.txt", "read");
    TFile *R2 = new TFile("RT_relative_Ni.txt", "read");
    TFile *R3 = new TFile("RT_relative_Cu.txt", "read");

    TFile *temperature = new TFile("RT_temperature.txt", "read");

    TFile *files[9] = {time, pot1, pot2, pot3, pot4, R1, R2, R3, temperature};

    TCanvas *c = new TCanvas("c", "Voltage Regimes", 1200, 600);
    c->Divide(3, 1);

    for (int i{0}; i < 3; ++i)
    {
        TGraph *gS = new TGraph(files[i]->GetName(), "%lg %lg %*lg %*lg %*lg");
        gS->SetTitle("Source; Time (s); Voltage (V)");
        TGraph *gL = new TGraph(files[i]->GetName(), "%lg %*lg %lg %*lg %*lg");
        gL->SetTitle("Inductor; Time (s); Voltage (V)");
        TGraph *gC = new TGraph(files[i]->GetName(), "%lg %*lg %*lg %lg %*lg");
        gC->SetTitle("Capacitor; Time (s); Voltage (V)");
        TGraph *gR = new TGraph(files[i]->GetName(), "%lg %*lg %*lg %*lg %lg");
        gR->SetTitle("Resistance; Time (s); Voltage (V)");

        gS->SetMarkerStyle(7);

        gL->SetMarkerStyle(7);
        gL->SetMarkerColor(kRed);
        //gL->SetLineColor(kRed);

        gC->SetMarkerStyle(7);
        gC->SetMarkerColor(kGreen);
        //gC->SetLineColor(kGreen);

        gR->SetMarkerStyle(7);
        gR->SetMarkerColor(kBlue);
       // gR->SetLineColor(kBlue);


        TMultiGraph *mg = new TMultiGraph();
        mg->Add(gL);
        mg->Add(gR);
        mg->Add(gC);
        mg->Add(gS);

        mg->SetTitle(title[i]);
        mg->GetXaxis()->SetTitle("Time (s)");
        mg->GetYaxis()->SetTitle("Voltage (V)");
        c->cd(i + 1);
        mg->Draw("APL");
        c->cd(i + 1)->BuildLegend();
    }

    return 0;
}
