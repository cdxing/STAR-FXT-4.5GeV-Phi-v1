

void syserrbar()
{
    TCanvas *c1 = new TCanvas("c1","c1",200,10,1024,768);
    c1->DrawFrame(0., -0.5, 6., 2);
    double x[5]    = {1, 2, 3, 4, 5};
    double zero[5] = {0, 0, 0, 0, 0};
    // data set (1) with stat and sys errors
    double py1[5]      = {1.2, 1.15, 1.19, 0.9, 1.4};
    double ey_stat1[5] = {0.2, 0.18, 0.17, 0.2, 0.4};
    double ey_sys1[5]  = {0.5, 0.71, 0.76, 0.5, 0.45};
    // Now draw data set (1)
    // We first have to draw it only with the stat errors
    gStyle->SetOptDate(0);
    TGraphErrors *graph1 = new TGraphErrors(5, x, py1, zero, ey_stat1);
    graph1->SetMarkerStyle(20);
    graph1->Draw("PZ");
    // Now we have to somehow depict the sys errors
    TGraphErrors *graph1_sys = new TGraphErrors(5, x, py1, zero, ey_sys1);
    graph1_sys->Draw("[]");
}
