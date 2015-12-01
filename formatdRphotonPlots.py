from ROOT import *
PhotonJetC = TCanvas("PhotonJetC", "Overlapping muon ET / jet ET", 800, 800)
#gStyle.SetOptTitle("E_{T}^{/mu} / E_{T}^{jet}")
gStyle.SetOptTitle(0)
PhotonJetC.SetLogy()
PhotonJetC.cd()
isoPlotsFile = TFile("no2015CflatTuple.root")
isoPlotsDir = isoPlotsFile.Get("Isolation")
PhotonJet1 = isoPlotsDir.Get("PhotonJetoverlapdR1")
PhotonJet2 = isoPlotsDir.Get("PhotonJetoverlapdR2")
PhotonJet3 = isoPlotsDir.Get("PhotonJetoverlapdR3")
PhotonJet4 = isoPlotsDir.Get("PhotonJetoverlapdR4")
PhotonJet1.Rebin(2)
PhotonJet2.Rebin(2)
PhotonJet3.Rebin(2)
PhotonJet4.Rebin(2)
PhotonJet1.GetXaxis().SetRangeUser(0,1.05)
PhotonJet1.GetXaxis().SetTitle("dR(#gamma,jet)")
PhotonJet1.GetYaxis().SetTitle("Events")
PhotonJet1.SetLineColor(kBlue)
PhotonJet2.SetLineColor(kGreen-2)
PhotonJet3.SetLineColor(kRed)
PhotonJet4.SetLineColor(kBlack)
PhotonJet1.Draw()
gPad.Update()
st1 = PhotonJet1.FindObject("stats")
st1.SetX1NDC(0.660804)
st1.SetX2NDC(0.771357)
st1.SetY1NDC(0.901682)
st1.SetY2NDC(0.994825)
PhotonJetC.Update()
PhotonJet2.Draw("SAMES")
gPad.Update()
st2 = PhotonJet2.FindObject("stats")
st2.SetX1NDC(0.771357)
st2.SetX2NDC(0.901682)
st2.SetY1NDC(0.901682)
st2.SetY2NDC(0.994825)
PhotonJetC.Update()
PhotonJet3.Draw("SAMES")
gPad.Update()
st3 = PhotonJet3.FindObject("stats")
st3.SetX1NDC(0.660804)
st3.SetX2NDC(0.771357)
st3.SetY1NDC(0.901682)
st3.SetY2NDC(0.8085391)
PhotonJetC.Update()
PhotonJet4.Draw("SAMES")
gPad.Update()
st4 = PhotonJet4.FindObject("stats")
st4.SetX1NDC(0.771357)
st4.SetX2NDC(0.901682)
st4.SetY1NDC(0.901682)
st4.SetY2NDC(0.8085391)
gPad.Update()
PhotonJetC.Update()
legend = TLegend(0.3, 0.85, 0.88, 0.65, "dR for photons and jets for E_{T}^{jet} ranges:", "brNDC")
legend.SetTextSize(0.02);
legend.SetLineColor(1);
legend.SetLineStyle(1);
legend.SetLineWidth(1);
legend.SetFillStyle(1001);
legend.SetFillColor(10);
legend.AddEntry(PhotonJet1,"50 < Jet E_{T} < 150 (GeV)","lp");
legend.AddEntry(PhotonJet2,"150 < Jet E_{T} < 250 (GeV)","lp");
legend.AddEntry(PhotonJet3,"250 < Jet E_{T} < 400 (GeV)","lp");
legend.AddEntry(PhotonJet4,"400 < Jet E_{T} (GeV)","lp");
legend.Draw()
gPad.Update()
PhotonJetC.Update()
legend.SetX1NDC(0.278894)
legend.SetX2NDC(0.658291)
legend.SetY1NDC(0.829237)
legend.SetY2NDC(0.997413)
gPad.Update()
PhotonJetC.Update()
gStyle.SetOptTitle(0)
gPad.Update()
PhotonJetC.Update()

outfile = TFile("output/PhotondROutFile.root", "RECREATE")
PhotonJetC.Write()
PhotonJetC.SaveAs("output/PhotondRPlot.pdf")

#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]
