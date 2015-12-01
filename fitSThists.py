# Script to fit ST distributions from the BHflatTuplizer and make them pretty. John Hakala 12/1/2015
# this guy takes two arguments: the lower bound of the fit range and the upper bound. Example:
# python fitSThists.py 1900 2900
from ROOT import *
from sys import argv
STexcComparisons = []
STincComparisons = []
upperExcPads = []
upperIncPads = []
lowerExcPads = []
lowerIncPads = []
PlotsFile = TFile("no2015CflatTuple.root")
PlotsDir = PlotsFile.Get("ST")
OutFile = TFile("output/FormattedFitPlots_Dec1.root", "RECREATE")
for i in range(2,12):
    stExcMEThist=PlotsDir.Get("stExc%02iHist"%i)
    stExc2METhist=PlotsDir.Get("stExc02Hist")
    stIncMEThist=PlotsDir.Get("stInc%02iHist"%i)

    stExcMHThist=PlotsDir.Get("stExc%02iHistMHT"%i)
    stExc2MHThist=PlotsDir.Get("stExc02HistMHT")
    stIncMHThist=PlotsDir.Get("stInc%02iHistMHT"%i)

    if (i == 2):
        STexcComparisons.append(TCanvas("stExc%02iCanvas"%i, "ST, N=%i"%i, 800, 800))
        STexcComparisons[i-2].cd()
        upperExcPads.append(TPad("Exc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
        upperExcPads[i-2].SetBottomMargin(0)
        upperExcPads[i-2].Draw()
        upperExcPads[i-2].cd()
        upperExcPads[i-2].SetLogy()
#        stExcMEThist.GetXaxis().SetRangeUser(1000, 7000)
#        stExcMEThist.Draw()

        stExcMEThist.GetXaxis().SetRangeUser(1000, 7000)
        stExcMEThist.GetXaxis().SetTitle("S_{T} (GeV)")
        stExcMEThist.GetYaxis().SetTitle("Events")
        stExcMEThist.SetLineColor(kBlue)
        stExcMEThist.Draw()

        stExcMHThist.GetXaxis().SetRangeUser(1000, 7000)
        stExcMHThist.GetXaxis().SetTitle("S_{T} (GeV)")
        stExcMHThist.GetYaxis().SetTitle("Events")
        stExcMHThist.SetLineColor(kRed)
        stExcMHThist.Draw("SAME")

        f1 = TF1("f1", "[0]/([1]*x)**[2]", 1000, 7000)
        f1.SetParameters(2.08e9, 4.28e-3, 6.7)
        for j in range(0, 20):
            stExcMEThist.Fit("f1", "", "", float(argv[1]), float(argv[2]) )
        f2 = TF1("f2", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
        f2.SetParameters(5e3, 2, 0.3, 0.2)
        for j in range(0, 20):
            stExcMEThist.Fit("f2", "", "", float(argv[1]), float(argv[2]) )
        f3 = TF1("f3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
        f3.SetParameters(5e9, 3, 5e3, 0.6)
        for j in range(0, 20):
            stExcMEThist.Fit("f3", "", "", float(argv[1]), float(argv[2]) )
        f1.SetLineColor(kCyan-5)
        f1.SetLineStyle(6)
        f1.Draw("SAME")
        f2.SetLineColor(kCyan-5)
        f2.SetLineStyle(6)
        f2.Draw("SAME")
        f3.SetLineColor(kCyan-5)
        f3.SetLineStyle(6)
        f3.Draw("SAME")

        f1MHT = TF1("f1MHT", "[0]/([1]*x)**[2]", 1000, 7000)
        f1MHT.SetParameters(2.08e9, 4.28e-3, 6.7)
        for j in range(0, 20):
            stExcMHThist.Fit("f1MHT", "", "", float(argv[1]), float(argv[2]) )
        f2MHT = TF1("f2MHT", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
        f2MHT.SetParameters(5e3, 2, 0.3, 0.2)
        for j in range(0, 20):
            stExcMHThist.Fit("f2MHT", "", "", float(argv[1]), float(argv[2]) )
        f3MHT = TF1("f3MHT", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
        f3MHT.SetParameters(5e9, 3, 5e3, 0.6)
        for j in range(0, 20):
            stExcMHThist.Fit("f3MHT", "", "", float(argv[1]), float(argv[2]) )
        f1MHT.SetLineColor(kOrange+5)
        f1MHT.SetLineStyle(6)
        f1MHT.Draw("SAME")
        f2MHT.SetLineColor(kOrange+5)
        f2MHT.SetLineStyle(6)
        f2MHT.Draw("SAME")
        f3MHT.SetLineColor(kOrange+5)
        f3MHT.SetLineStyle(6)
        f3MHT.Draw("SAME")


        STexcComparisons[i-2].cd()
        lowerExcPads.append(TPad("Exc%02iratiopad"%i, "ratiopad1", 0, 0.04, 1, 0.3))
        lowerExcPads[i-2].SetTopMargin(2)
        lowerExcPads[i-2].SetBottomMargin(0.2)
        lowerExcPads[i-2].Draw()
        lowerExcPads[i-2].cd()
        lowerExcPads[i-2].Draw()
        STexcComparisons[i-2].Write()

        STincComparisons.append(TCanvas("stInc%02iCanvas"%i, "ST, N>=%i"%i, 800, 800))
        STincComparisons[i-2].cd()
        upperIncPads.append(TPad("Inc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
        upperIncPads[i-2].Draw()
        upperIncPads[i-2].cd()
        upperIncPads[i-2].Draw()
        STincComparisons[i-2].cd()
        lowerIncPads.append(TPad("Inc%02iratiopad"%i, "incratiopad1", 0, 0.04, 1, 0.3))
        lowerIncPads[i-2].SetTopMargin(0)
        lowerIncPads[i-2].SetBottomMargin(0.2)
        lowerIncPads[i-2].Draw()
        STincComparisons[i-2].Write()

    else:
        STexcComparisons.append(TCanvas("stExc%02iCanvas"%i, "ST, N=%i"%i, 800, 800))
        STexcComparisons[i-2].cd()
        upperExcPads.append(TPad("Exc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
        upperExcPads[i-2].SetBottomMargin(2)
        upperExcPads[i-2].Draw()
        upperExcPads[i-2].cd()
        upperExcPads[i-2].SetLogy()
        stExcMEThist.GetXaxis().SetRangeUser(1000, 7000)
        stExcMEThist.GetXaxis().SetTitle("S_{T} (GeV)")
        stExcMEThist.GetYaxis().SetTitle("Events")
        stExcMEThist.Draw()
        stExcMHThist.GetXaxis().SetRangeUser(1000, 7000)
        stExcMHThist.GetXaxis().SetTitle("S_{T} (GeV)")
        stExcMHThist.GetYaxis().SetTitle("Events")
        stExcMHThist.SetLineColor(kRed)
        stExcMHThist.Draw("SAME")
        legend = TLegend(0.35, 0.68, 0.58, 0.85, "Exclusive ST distributions for n=%i"%i, "brNDC")
        legend.SetTextSize(0.02);
        legend.SetLineColor(1);
        legend.SetLineWidth(1);
        legend.SetFillStyle(1001);
        legend.SetFillColor(10);
        legend.AddEntry(stExcMEThist,"ST using MET","lp");
        legend.AddEntry(stExcMHThist,"ST using MHT","lp");
        legend.Draw()

        STexcComparisons[i-2].cd()
        lowerExcPads.append(TPad("Exc%02iratiopad"%i, "ratiopad1", 0, 0.04, 1, 0.3))
        lowerExcPads[i-2].SetTopMargin(5)
        lowerExcPads[i-2].SetBottomMargin(0.2)
        lowerExcPads[i-2].Draw()
        lowerExcPads[i-2].cd()

        stExcMETRatio = stExcMEThist.Clone("stExcMETRatio")
        stExcMETRatio.GetXaxis().SetLabelSize(.08)
        stExcMETRatio.GetXaxis().SetTitle("")
        stExcMETRatio.GetYaxis().SetTitle("Ratio to n=2")
        stExcMETRatio.GetYaxis().SetLabelSize(.075)
        stExcMETRatio.GetYaxis().SetTitleSize(.1)
        stExcMETRatio.GetYaxis().SetTitleOffset(.3)
        stExcMETRatio.Divide(stExc2METhist)
        stExcMETRatio.SetTitle("")
        stExcMETRatio.SetLineColor(kBlue)
        stExcMETRatio.SetStats(0)
        stExcMETRatio.Draw()
        stExcMHTRatio = stExcMHThist.Clone("stExcMHTRatio")
        stExcMHTRatio.GetXaxis().SetLabelSize(.08)
        stExcMHTRatio.GetXaxis().SetTitle("")
        stExcMHTRatio.GetYaxis().SetTitle("Ratio to n=2")
        stExcMHTRatio.GetYaxis().SetLabelSize(.075)
        stExcMHTRatio.GetYaxis().SetTitleSize(.1)
        stExcMHTRatio.GetYaxis().SetTitleOffset(.3)
        stExcMHTRatio.Divide(stExc2MHThist)
        stExcMHTRatio.SetTitle("")
        stExcMHTRatio.SetLineColor(kRed)
        stExcMHTRatio.SetStats(0)
        stExcMHTRatio.Draw("SAME")
        STexcComparisons[i-2].Write()

        STincComparisons.append(TCanvas("stInc%02iCanvas"%i, "ST, N>=%i"%i, 800, 800))
        STincComparisons[i-2].cd()
        upperIncPads.append(TPad("Inc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
        upperIncPads[i-2].SetBottomMargin(2)
        upperIncPads[i-2].Draw()
        upperIncPads[i-2].cd()
        upperIncPads[i-2].SetLogy()
        stIncMEThist.GetXaxis().SetRangeUser(1000, 7000)
        stIncMEThist.GetXaxis().SetTitle("S_{T} (GeV)")
        stIncMEThist.GetYaxis().SetTitle("Events")
        stIncMEThist.SetLineColor(kBlue)
        stIncMEThist.Draw()
        stIncMHThist.GetXaxis().SetRangeUser(1000, 7000)
        stIncMHThist.GetXaxis().SetTitle("S_{T} (GeV)")
        stIncMHThist.GetYaxis().SetTitle("Events")
        stIncMHThist.SetLineColor(kRed)
        stIncMHThist.Draw("SAME")
        legend = TLegend(0.35, 0.68, 0.58, 0.85, "Inclusive ST distributions for n>=%i"%i, "brNDC")
        legend.SetTextSize(0.02);
        legend.SetLineColor(1);
        legend.SetLineStyle(1);
        legend.SetLineWidth(1);
        legend.SetFillStyle(1001);
        legend.SetFillColor(10);
        legend.AddEntry(stExcMEThist,"ST using MET","lp");
        legend.AddEntry(stExcMHThist,"ST using MHT","lp");
        legend.Draw()

        STincComparisons[i-2].cd()
        lowerIncPads.append(TPad("Inc%02iratiopad"%i, "incratiopad1", 0, 0.04, 1, 0.3))
        lowerIncPads[i-2].SetTopMargin(5)
        lowerIncPads[i-2].SetBottomMargin(0.2)
        lowerIncPads[i-2].Draw()
        lowerIncPads[i-2].cd()

        stIncMETRatio = stIncMEThist.Clone("stIncMETRatio")
        stIncMETRatio.GetXaxis().SetLabelSize(.08)
        stIncMETRatio.GetYaxis().SetTitle("Ratio to n=2")
        stIncMETRatio.GetYaxis().SetLabelSize(.075)
        stIncMETRatio.GetYaxis().SetTitleSize(.1)
        stIncMETRatio.GetYaxis().SetTitleOffset(.3)
        stIncMETRatio.GetXaxis().SetTitle("")
        stIncMETRatio.Divide(stExc2METhist)
        stIncMETRatio.SetStats(0)
        stIncMETRatio.SetTitle("")
        stIncMETRatio.SetLineColor(kBlue)
        stIncMETRatio.Draw()
        stIncMHTRatio = stIncMHThist.Clone("stIncMHTRatio")
        stIncMHTRatio.GetXaxis().SetLabelSize(.08)
        stIncMHTRatio.GetYaxis().SetTitle("Ratio to n=2")
        stIncMHTRatio.GetYaxis().SetLabelSize(.075)
        stIncMHTRatio.GetYaxis().SetTitleSize(.1)
        stIncMHTRatio.GetYaxis().SetTitleOffset(.3)
        stIncMHTRatio.GetXaxis().SetTitle("")
        stIncMHTRatio.Divide(stExc2MHThist)
        stIncMHTRatio.SetStats(0)
        stIncMHTRatio.SetTitle("")
        stIncMHTRatio.SetLineColor(kRed)
        stIncMHTRatio.Draw("SAME")
        STincComparisons[i-2].Write()

#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]
