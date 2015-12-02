# Script to fit ST distributions from the BHflatTuplizer and make them pretty. John Hakala 12/1/2015
# this guy takes four arguments: the input BHflatTuple filename, the output filename, the name of the file defining the fit ranges, and either useMET or useMHT
# python fitSThists.py myBHflatTuple.root myOutputFile.root myFitNormRangesFile.txt useMHT
from ROOT import *
from fitAndNormRanges import *
from sys import argv
STexcComparisons = []
STincComparisons = []
upperExcPads = []
upperIncPads = []
lowerExcPads = []
lowerIncPads = []
PlotsFile = TFile(argv[1])
PlotsDir = PlotsFile.Get("ST")
OutFile = TFile(argv[2], "RECREATE")
fitNormRanges = FitAndNormRange(argv[3])
fitNormRanges.showFitRanges()
fitNormRanges.showNormRanges()
if argv[4] == "useMET" or argv[4] == "useBoth":
    f1 = TF1("f1", "[0]/([1]*x)**[2]", 1000, 7000)
    f2 = TF1("f2", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
    f3 = TF1("f3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
if argv[4] == "useMHT" or argv[4] == "useBoth":
    f1MHT = TF1("f1MHT", "[0]/([1]*x)**[2]", 1000, 7000)
    f2MHT = TF1("f2MHT", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
    f3MHT = TF1("f3MHT", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
for i in range(2,12):
    if argv[4] == "useMET" or argv[4] == "useBoth":
        stExcMEThist=PlotsDir.Get("stExc%02iHist"%i)
        stExc2METhist=PlotsDir.Get("stExc02Hist")
        stIncMEThist=PlotsDir.Get("stInc%02iHist"%i)

    if argv[4] == "useMHT" or argv[4] == "useBoth":
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

        if argv[4] == "useMET" or argv[4] == "useBoth":
            stExcMEThist.GetXaxis().SetRangeUser(1000, 7000)
            stExcMEThist.GetXaxis().SetTitle("S_{T} (GeV)")
            stExcMEThist.GetYaxis().SetTitle("Events")
            stExcMEThist.SetLineColor(kBlue)
            stExcMEThist.Draw()

        if argv[4] == "useMHT" or argv[4] == "useBoth":
            stExcMHThist.GetXaxis().SetRangeUser(1000, 7000)
            stExcMHThist.GetXaxis().SetTitle("S_{T} (GeV)")
            stExcMHThist.GetYaxis().SetTitle("Events")
            stExcMHThist.SetLineColor(kRed)
            if argv[4] == "useBoth":
                stExcMHThist.Draw("SAME")
            else:
                stExcMHThist.Draw()

        if argv[4] == "useMET" or argv[4] == "useBoth":
            f1.SetParameters(2.08e9, 4.28e-3, 6.7)
            for j in range(0, 20):
                stExcMEThist.Fit("f1", "0", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
            print "done fitting f1"
            f2.SetParameters(1.7e12, 2.1, .16, 0.62)
            for j in range(0, 20):
                stExcMEThist.Fit("f2", "0", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
            print "done fitting f2"
            f3.SetParameters(2e27, 4e5, -4e2, 3.6)
            for j in range(0, 20):
                stExcMEThist.Fit("f3", "0", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
            print "done fitting f3"
            f1.SetLineColor(kCyan-5)
            f1.SetLineStyle(6)
            f1.Draw("SAME")
            f2.SetLineColor(kCyan-5)
            f2.SetLineStyle(6)
            f2.Draw("SAME")
            f3.SetLineColor(kCyan-5)
            f3.SetLineStyle(6)
            f3.Draw("SAME")

        if argv[4] == "useMHT" or argv[4] == "useBoth":
            f1MHT.SetParameters(2.08e9, 4.28e-3, 6.7)
            for j in range(0, 20):
                stExcMHThist.Fit("f1MHT", "0", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2")  )
            f2MHT.SetParameters(1.7e12, 2.1, .16, 0.62)
            for j in range(0, 20):
                stExcMHThist.Fit("f2MHT", "0", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
            f3MHT.SetParameters(2e27, 4e5, -4e2, 3.6)
            for j in range(0, 20):
                stExcMHThist.Fit("f3MHT", "0", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2")  )
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
        if argv[4] == "useMET" or argv[4] == "useBoth":
            stExcMEThist.GetXaxis().SetRangeUser(1000, 7000)
            stExcMEThist.GetXaxis().SetTitle("S_{T} (GeV)")
            stExcMEThist.GetYaxis().SetTitle("Events")
            stExcMEThist.Draw()
            lowerNormBin = stExcMEThist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("exc%i"%i)))
            upperNormBin = stExcMEThist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("exc%i"%i)))
            lowerNormEdge = stExcMEThist.GetXaxis().GetBinLowEdge(lowerNormBin)
            upperNormEdge = stExcMEThist.GetXaxis().GetBinLowEdge(upperNormBin)
            normBinTotal = 0;
            for normbin in range(lowerNormBin, upperNormBin):
                normBinTotal+=stExcMEThist.GetBinContent(normbin)
            normfactor =  (normBinTotal/f1.Integral(lowerNormEdge, upperNormEdge))*stExcMEThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
            #print "normalization factor is: %f" % normfactor
            f1Normalized = f1.Clone()
            f1Normalized.SetParameter(0, f1.GetParameter(0)*normfactor)
            f1Normalized.Draw("SAME")
            f2Normalized = f2.Clone()
            f2Normalized.SetParameter(0, f2.GetParameter(0)*normfactor)
            f2Normalized.Draw("SAME")
            f3Normalized = f3.Clone()
            f3Normalized.SetParameter(0, f3.GetParameter(0)*normfactor)

            f3Normalized.Draw("SAME")
        if argv[4] == "useMHT" or argv[4] == "useBoth":
            stExcMHThist.GetXaxis().SetRangeUser(1000, 7000)
            stExcMHThist.GetXaxis().SetTitle("S_{T} (GeV)")
            stExcMHThist.GetYaxis().SetTitle("Events")
            stExcMHThist.SetLineColor(kRed)
            if argv[4] == "useBoth":
                stExcMHThist.Draw("SAME")
            else:
                stExcMHThist.Draw()
            lowerNormBin = stExcMHThist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("inc%i"%i)))
            upperNormBin = stExcMHThist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("inc%i"%i)))
            lowerNormEdge = stExcMHThist.GetXaxis().GetBinLowEdge(lowerNormBin)
            upperNormEdge = stExcMHThist.GetXaxis().GetBinLowEdge(upperNormBin)
            normBinTotal = 0;
            for normbin in range(lowerNormBin, upperNormBin):
                normBinTotal+=stExcMHThist.GetBinContent(normbin)
            normfactorMHT =  (normBinTotal/f1MHT.Integral(lowerNormEdge, upperNormEdge))*stExcMHThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
            #print "normalization factor is: %f" % normfactorf1MHT
            f1MHTNormalized = f1MHT.Clone()
            f1MHTNormalized.SetParameter(0, f1MHT.GetParameter(0)*normfactorMHT)
            f1MHTNormalized.Draw("SAME")
            f2MHTNormalized = f2MHT.Clone()
            f2MHTNormalized.SetParameter(0, f2MHT.GetParameter(0)*normfactorMHT)
            f2MHTNormalized.Draw("SAME")
            f3MHTNormalized = f3MHT.Clone()
            f3MHTNormalized.SetParameter(0, f3MHT.GetParameter(0)*normfactorMHT)
            f3MHTNormalized.Draw("SAME")


        legend = TLegend(0.35, 0.68, 0.58, 0.85, "Exclusive ST distributions for n=%i"%i, "brNDC")
        legend.SetTextSize(0.02);
        legend.SetLineColor(1);
        legend.SetLineWidth(1);
        legend.SetFillStyle(1001);
        legend.SetFillColor(10);
        if argv[4] == "useMET" or argv[4] == "useBoth":
            legend.AddEntry(stExcMEThist,"ST using MET","lp");
        if argv[4] == "useMHT" or argv[4] == "useBoth":
            legend.AddEntry(stExcMHThist,"ST using MHT","lp");
        legend.Draw()

        STexcComparisons[i-2].cd()
        lowerExcPads.append(TPad("Exc%02iratiopad"%i, "ratiopad1", 0, 0.04, 1, 0.3))
        lowerExcPads[i-2].SetTopMargin(5)
        lowerExcPads[i-2].SetBottomMargin(0.2)
        lowerExcPads[i-2].Draw()
        lowerExcPads[i-2].cd()

        if argv[4] == "useMET" or argv[4] == "useBoth":
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
        if argv[4] == "useMHT" or argv[4] == "useBoth":
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
            if argv[4]=="useBoth":
                stExcMHTRatio.Draw("SAME")
            else:
                stExcMHTRatio.Draw()
        STexcComparisons[i-2].Write()

        STincComparisons.append(TCanvas("stInc%02iCanvas"%i, "ST, N>=%i"%i, 800, 800))
        STincComparisons[i-2].cd()
        upperIncPads.append(TPad("Inc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
        upperIncPads[i-2].SetBottomMargin(2)
        upperIncPads[i-2].Draw()
        upperIncPads[i-2].cd()
        upperIncPads[i-2].SetLogy()
        if argv[4] == "useMET" or argv[4] == "useBoth":
            stIncMEThist.GetXaxis().SetRangeUser(1000, 7000)
            stIncMEThist.GetXaxis().SetTitle("S_{T} (GeV)")
            stIncMEThist.GetYaxis().SetTitle("Events")
            stIncMEThist.SetLineColor(kBlue)
            stIncMEThist.Draw()
            lowerNormBin = stIncMEThist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("exc%i"%i)))
            upperNormBin = stIncMEThist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("exc%i"%i)))
            lowerNormEdge = stIncMEThist.GetXaxis().GetBinLowEdge(lowerNormBin)
            upperNormEdge = stIncMEThist.GetXaxis().GetBinLowEdge(upperNormBin)
            normBinTotal = 0;
            for normbin in range(lowerNormBin, upperNormBin):
                normBinTotal+=stIncMEThist.GetBinContent(normbin)
            normfactorMET =  (normBinTotal/f1.Integral(lowerNormEdge, upperNormEdge))*stIncMEThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
            #print "normalization factor is: %f" % normfactor
            f1Normalized = f1.Clone()
            f1Normalized.SetParameter(0, f1.GetParameter(0)*normfactorMET)
            f1Normalized.Draw("SAME")
            f2Normalized = f2.Clone()
            f2Normalized.SetParameter(0, f2.GetParameter(0)*normfactorMET)
            f2Normalized.Draw("SAME")
            f3Normalized = f3.Clone()
            f3Normalized.SetParameter(0, f3.GetParameter(0)*normfactorMET)
            f3Normalized.Draw("SAME")

        if argv[4] == "useMHT" or argv[4] == "useBoth":
            stIncMHThist.GetXaxis().SetRangeUser(1000, 7000)
            stIncMHThist.GetXaxis().SetTitle("S_{T} (GeV)")
            stIncMHThist.GetYaxis().SetTitle("Events")
            stIncMHThist.SetLineColor(kRed)
            stIncMHThist.Draw("SAME")
            lowerNormBin = stIncMHThist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("inc%i"%i)))
            upperNormBin = stIncMHThist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("inc%i"%i)))
            lowerNormEdge = stIncMHThist.GetXaxis().GetBinLowEdge(lowerNormBin)
            upperNormEdge = stIncMHThist.GetXaxis().GetBinLowEdge(upperNormBin)
            normBinTotal = 0;
            for normbin in range(lowerNormBin, upperNormBin):
                normBinTotal+=stIncMHThist.GetBinContent(normbin)
            normfactorMHT =  (normBinTotal/f1MHT.Integral(lowerNormEdge, upperNormEdge))*stIncMHThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
            #print "normalization factor is: %f" % normfactor
            f1MHTNormalized = f1MHT.Clone()
            f1MHTNormalized.SetParameter(0, f1MHT.GetParameter(0)*normfactorMHT)
            f1MHTNormalized.Draw("SAME")
            f2MHTNormalized = f2MHT.Clone()
            f2MHTNormalized.SetParameter(0, f2MHT.GetParameter(0)*normfactorMHT)
            f2MHTNormalized.Draw("SAME")
            f3MHTNormalized = f3MHT.Clone()
            f3MHTNormalized.SetParameter(0, f3MHT.GetParameter(0)*normfactorMHT)
            f3MHTNormalized.Draw("SAME")

        legend = TLegend(0.35, 0.68, 0.58, 0.85, "Inclusive ST distributions for n>=%i"%i, "brNDC")
        legend.SetTextSize(0.02);
        legend.SetLineColor(1);
        legend.SetLineStyle(1);
        legend.SetLineWidth(1);
        legend.SetFillStyle(1001);
        legend.SetFillColor(10);
        if argv[4] == "useMET" or argv[4] == "useBoth":
            legend.AddEntry(stExcMEThist,"ST using MET","lp");
        if argv[4] == "useMHT" or argv[4] == "useBoth":
            legend.AddEntry(stExcMHThist,"ST using MHT","lp");
        legend.Draw()

        STincComparisons[i-2].cd()
        lowerIncPads.append(TPad("Inc%02iratiopad"%i, "incratiopad1", 0, 0.04, 1, 0.3))
        lowerIncPads[i-2].SetTopMargin(5)
        lowerIncPads[i-2].SetBottomMargin(0.2)
        lowerIncPads[i-2].Draw()
        lowerIncPads[i-2].cd()

        if argv[4] == "useMET" or argv[4] == "useBoth":
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
        if argv[4] == "useMHT" or argv[4] == "useBoth":
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
            if argv[4] == "useBoth":
                stIncMHTRatio = stIncMHThist.Clone("stIncMHTRatio")
            else:
                stIncMHTRatio.Draw("SAME")
        STincComparisons[i-2].Write()

#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]
