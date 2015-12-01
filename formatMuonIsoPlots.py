from ROOT import *
MuonJetC = TCanvas("MuonJetC", "Overlapping muon ET / jet ET", 800, 800)
#gStyle.SetOptTitle("E_{T}^{/mu} / E_{T}^{jet}")
gStyle.SetOptTitle(0)
MuonJetC.SetLogy()
MuonJetC.cd()
isoPlotsFile = TFile("no2015CflatTuple.root")
isoPlotsDir = isoPlotsFile.Get("Isolation")
MuonJet1 = isoPlotsDir.Get("MuonJetIso1")
MuonJet2 = isoPlotsDir.Get("MuonJetIso2")
MuonJet3 = isoPlotsDir.Get("MuonJetIso3")
MuonJet4 = isoPlotsDir.Get("MuonJetIso4")
MuonJet1.Rebin(2)
MuonJet2.Rebin(2)
MuonJet3.Rebin(2)
MuonJet4.Rebin(2)
MuonJet1.GetXaxis().SetRangeUser(0,1.05)
MuonJet1.GetXaxis().SetTitle("E_{T}^{#mu}/E_{T}^{jet}")
MuonJet1.GetYaxis().SetTitle("Events")
MuonJet1.SetLineColor(kBlue)
MuonJet2.SetLineColor(kGreen-2)
MuonJet3.SetLineColor(kRed)
MuonJet4.SetLineColor(kBlack)
MuonJet1.Draw()
gPad.Update()
st1 = MuonJet1.FindObject("stats")
st1.SetX1NDC(0.660804-.2)
st1.SetX2NDC(0.771357-.2)
st1.SetY1NDC(0.901682)
st1.SetY2NDC(0.994825)
MuonJetC.Update()
MuonJet2.Draw("SAMES")
gPad.Update()
st2 = MuonJet2.FindObject("stats")
st2.SetX1NDC(0.771357-.2)
st2.SetX2NDC(0.901682-.2)
st2.SetY1NDC(0.901682)
st2.SetY2NDC(0.994825)
MuonJetC.Update()
MuonJet3.Draw("SAMES")
gPad.Update()
st3 = MuonJet3.FindObject("stats")
st3.SetX1NDC(0.660804-.2)
st3.SetX2NDC(0.771357-.2)
st3.SetY1NDC(0.901682)
st3.SetY2NDC(0.8085391)
MuonJetC.Update()
MuonJet4.Draw("SAMES")
gPad.Update()
st4 = MuonJet4.FindObject("stats")
st4.SetX1NDC(0.771357-.2)
st4.SetX2NDC(0.901682-.2)
st4.SetY1NDC(0.901682)
st4.SetY2NDC(0.8085391)
gPad.Update()
MuonJetC.Update()
legend = TLegend(0.5, 0.85, 0.88, 0.65, "Overlapping E_{T}^{#mu} / E_{T}^{jet} for E_{T}^{jet} ranges:", "brNDC")
legend.SetTextSize(0.02);
legend.SetLineColor(1);
legend.SetLineStyle(1);
legend.SetLineWidth(1);
legend.SetFillStyle(1001);
legend.SetFillColor(10);
legend.AddEntry(MuonJet1,"50 < Jet E_{T} < 150","lp");
legend.AddEntry(MuonJet2,"150 < Jet E_{T} < 250","lp");
legend.AddEntry(MuonJet3,"250 < Jet E_{T} < 400","lp");
legend.AddEntry(MuonJet4,"400 < Jet E_{T}","lp");
legend.Draw()
gPad.Update()
MuonJetC.Update()
legend.SetX1NDC(0.278894-.12)
legend.SetX2NDC(0.658291-.2)
legend.SetY1NDC(0.829237)
legend.SetY2NDC(0.997413)
gPad.Update()
MuonJetC.Update()
gStyle.SetOptTitle(0)
gPad.Update()
MuonJetC.Update()

outfile = TFile("output/IsoPlotsOutFile.root", "RECREATE")
MuonJetC.Write()
MuonJetC.SaveAs("output/MuonIsoPlot.pdf")

#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]
