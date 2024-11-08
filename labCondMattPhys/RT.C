int RT()
{
    // Read data files
    TFile *data_V = new TFile("RT_data_V.txt", "read");
    TFile *data_R = new TFile("RT_data_R.txt", "read");

    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "Temperature vs Time", 1200, 600);
    TCanvas *c2 = new TCanvas("c2", "Potential vs Time", 1200, 600);
    TCanvas *c3 = new TCanvas("c3", "Relative resistances vs Temperature", 1200, 600);

    // Temperature graph
    TGraph *temp = new TGraph("RT_data_V.txt", "%lg %*lg %*lg %*lg %lg");
    temp->SetTitle("T(t); Time (s); Temperature (K)");
    temp->SetLineColor(kGreen); 

    // Potentials graphs
    TGraph *pot1 = new TGraph("RT_data_V.txt", "%lg %lg %*lg %*lg %*lg");
    pot1->SetTitle("V1(t); Time (s); Voltage (V)");
    pot1->SetLineColor(kBlue);
    
    TGraph *pot2 = new TGraph("RT_data_V.txt", "%lg %*lg %lg %*lg %*lg");
    pot2->SetTitle("V2(t); Time (s); Voltage (V)");
    pot2->SetLineColor(kRed); 
   
    TGraph *pot3 = new TGraph("RT_data_V.txt", "%lg %*lg %*lg %lg %*lg");
    pot3->SetTitle("V3(t); Time (s); Voltage (V)");
    pot3->SetLineColor(kViolet); 

    TGraph *pot4 = new TGraph("RT_data_V.txt", "%lg %*lg %*lg %*lg %lg");
    pot4->SetTitle("V4(t); Time (s); Voltage (V)");
    pot4->SetLineColor(kGreen); 

    // Create multigraph for potentials
    TMultiGraph *mg_pot = new TMultiGraph();
    mg_pot->Add(pot1);
    mg_pot->Add(pot2);
    mg_pot->Add(pot3);
    mg_pot->Add(pot4);

    // Cosmetics for multigraph
    mg_pot->SetTitle("Potentials versus time");
    mg_pot->GetXaxis()->SetTitle("Time (s)");
    mg_pot->GetYaxis()->SetTitle("Voltage (V)");
    
    // Relative resistances graphs
    TGraph *R1 = new TGraph("RT_data_R.txt", "%lg %lg %*lg %*lg");
    R1->SetTitle("R1(T)/R1_0; Temperature (K); Relative resistance (a.u.)");
    R1->SetLineColor(kBlue); // Set color for first resistance graph

    TGraph *R2 = new TGraph("RT_data_R.txt", "%lg %*lg %lg %*lg");
    R2->SetTitle("R2(T)/R2_0; Temperature (K); Relative resistance (a.u.)");
    R2->SetLineColor(kRed); // Set color for second resistance graph

    TGraph *R3 = new TGraph("RT_data_R.txt", "%lg %*lg %*lg %lg");
    R3->SetTitle("R3(T)/R3_0; Temperature (K); Relative resistance (a.u.)");
    R3->SetLineColor(kViolet); // Set color for third resistance graph

    // Create multigraph for relative resistances
    TMultiGraph *mg_res = new TMultiGraph();
    mg_res->Add(R1);
    mg_res->Add(R2);
    mg_res->Add(R3);
    
    // Cosmetics for multigraph
    mg_res->SetTitle("Relative resistances versus temperature");
    mg_res->GetXaxis()->SetTitle("Temperature (K)");
    mg_res->GetYaxis()->SetTitle("Relative resistances (a.u.)");
    
    // Draw graphs on canvas 
    c1->cd();
    temp->Draw("APL");
    
    c2->cd();
    mg_pot->Draw("APL");
    c2->BuildLegend();

    c3->cd();
    mg_res->Draw("APL");
    c3->BuildLegend();

    return 0;
}