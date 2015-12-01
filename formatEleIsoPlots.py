from ROOT import *
ElectronJetC = TCanvas("ElectronJetC", "Overlapping electron ET / jet ET", 800, 800)
#gStyle.SetOptTitle("E_{T}^{/mu} / E_{T}^{jet}")
gStyle.SetOptTitle(0)
ElectronJetC.SetLogy()
ElectronJetC.cd()
isoPlotsFile = TFile("no2015CflatTuple.root")
isoPlotsDir = isoPlotsFile.Get("Isolation")
ElectronJet1 = isoPlotsDir.Get("ElectronJetIso1")
ElectronJet2 = isoPlotsDir.Get("ElectronJetIso2")
ElectronJet3 = isoPlotsDir.Get("ElectronJetIso3")
ElectronJet4 = isoPlotsDir.Get("ElectronJetIso4")
ElectronJet1.Rebin(2)
ElectronJet2.Rebin(2)
ElectronJet3.Rebin(2)
ElectronJet4.Rebin(2)
ElectronJet1.GetXaxis().SetRangeUser(0,1.05)
ElectronJet1.GetXaxis().SetTitle("E_{T}^{e}/E_{T}^{jet}")
ElectronJet1.GetYaxis().SetTitle("Events")
ElectronJet1.SetLineColor(kBlue)
ElectronJet2.SetLineColor(kGreen-2)
ElectronJet3.SetLineColor(kRed)
ElectronJet4.SetLineColor(kBlack)
ElectronJet1.Draw()
gPad.Update()
st1 = ElectronJet1.FindObject("stats")
st1.SetX1NDC(0.660804-.2)
st1.SetX2NDC(0.771357-.2)
st1.SetY1NDC(0.901682)
st1.SetY2NDC(0.994825)
ElectronJetC.Update()
ElectronJet2.Draw("SAMES")
gPad.Update()
st2 = ElectronJet2.FindObject("stats")
st2.SetX1NDC(0.771357-.2)
st2.SetX2NDC(0.901682-.2)
st2.SetY1NDC(0.901682)
st2.SetY2NDC(0.994825)
ElectronJetC.Update()
ElectronJet3.Draw("SAMES")
gPad.Update()
st3 = ElectronJet3.FindObject("stats")
st3.SetX1NDC(0.660804-.2)
st3.SetX2NDC(0.771357-.2)
st3.SetY1NDC(0.901682)
st3.SetY2NDC(0.8085391)
ElectronJetC.Update()
ElectronJet4.Draw("SAMES")
gPad.Update()
st4 = ElectronJet4.FindObject("stats")
st4.SetX1NDC(0.771357-.2)
st4.SetX2NDC(0.901682-.2)
st4.SetY1NDC(0.901682)
st4.SetY2NDC(0.8085391)
gPad.Update()
ElectronJetC.Update()
legend = TLegend(0.5, 0.85, 0.88, 0.65, "Overlapping E_{T}^{e} / E_{T}^{jet} for E_{T}^{jet} ranges:", "brNDC")
legend.SetTextSize(0.02);
legend.SetLineColor(1);
legend.SetLineStyle(1);
legend.SetLineWidth(1);
legend.SetFillStyle(1001);
legend.SetFillColor(10);
legend.AddEntry(ElectronJet1,"50 < Jet E_{T} < 150","lp");
legend.AddEntry(ElectronJet2,"150 < Jet E_{T} < 250","lp");
legend.AddEntry(ElectronJet3,"250 < Jet E_{T} < 400","lp");
legend.AddEntry(ElectronJet4,"400 < Jet E_{T}","lp");
legend.Draw()
gPad.Update()
ElectronJetC.Update()
legend.SetX1NDC(0.278894-.12)
legend.SetX2NDC(0.658291-.2)
legend.SetY1NDC(0.829237)
legend.SetY2NDC(0.997413)
gPad.Update()
ElectronJetC.Update()
gStyle.SetOptTitle(0)
gPad.Update()
ElectronJetC.Update()

outfile = TFile("output/EleIsoPlotsOutFile.root", "RECREATE")
ElectronJetC.Write()
ElectronJetC.SaveAs("output/ElectronIsoPlot.pdf")

#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]
