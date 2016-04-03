# Script to fit ST distributions from the BHflatTuplizer, make them pretty, and spit out text files for the combine tool. John Hakala 12/1/2015
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
Chi2List     = []
Funcs	     = {}
PlotsFile = TFile(argv[1])
PlotsDir = PlotsFile.Get("ST")
OutFile = TFile("output/%s"%argv[2], "RECREATE")
fitNormRanges = FitAndNormRange(argv[3])
fitNormRanges.showFitRanges()
fitNormRanges.showNormRanges()
if argv[4] == "useMET" :
    f1 = TF1("f1", "[0]/([1]*x)**[2]", 1000, 7000)
    f2 = TF1("f2", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
    f3 = TF1("f3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
    f1_exc3 = TF1("f1_exc3", "[0]/([1]*x)**[2]", 1000, 7000)
    f2_exc3 = TF1("f2_exc3", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
    f3_exc3 = TF1("f3_exc3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
    fLow=TF1("fLow", "[0]*(2*[1]/([2]*x)**[3]-[4]/([5] + [6]*x + x**2)**[7])", 1000, 7000)
    Funcs ={0:f1,1:f2,2:f3,3:f1_exc3,4:f2_exc3,5:f3_exc3,6:fLow}
if argv[4] == "useMHT" :
    f1MHT = TF1("f1MHT", "[0]/([1]*x)**[2]", 1000, 7000)
    f2MHT = TF1("f2MHT", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
    f3MHT = TF1("f3MHT", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
    f1MHT_exc3 = TF1("f1MHT_exc3", "[0]/([1]*x)**[2]", 1000, 7000)
    f2MHT_exc3 = TF1("f2MHT_exc3", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
    f3MHT_exc3 = TF1("f3MHT_exc3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
    Funcs ={0:f1MHT,1:f2MHT,2:f3MHT,3:f1MHT_exc3,4:f2MHT_exc3,5:f3MHT_exc3}
for i in range(2,4):
    if argv[4] == "useMET" :
        stExcMEThist=PlotsDir.Get("stExc%02iHist"%i)
        stExc2METhist=PlotsDir.Get("stExc02Hist")
        stIncMEThist=PlotsDir.Get("stInc%02iHist"%i)

    if argv[4] == "useMHT" :
        stExcMHThist=PlotsDir.Get("stExc%02iHistMHT"%i)
        stExc2MHThist=PlotsDir.Get("stExc02HistMHT")
        stIncMHThist=PlotsDir.Get("stInc%02iHistMHT"%i)

    if (i == 2 or i==3 ):
        if argv[4] == "useMET" :
            if i==2:
                f1.SetParameters(2.08e9, 4.28e-3, 6.7)
                for j in range(0, 20):
                    stExcMEThist.Fit("f1", "Q0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                stExcMEThist.Fit("f1", "E0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                print "done fitting f1"
		Chi2List.append(f1.GetChisquare())

                f2.SetParameters(1.7e12, 2.1, .16, 0.62)
                for j in range(0, 20):
                    stExcMEThist.Fit("f2", "Q0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                stExcMEThist.Fit("f2", "E0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                print "done fitting f2"
		Chi2List.append(f2.GetChisquare())

                f3.SetParameters(2e27, 4e5, -4e2, 3.6)
                for j in range(0, 20):
                    stExcMEThist.Fit("f3", "Q0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                stExcMEThist.Fit("f3", "E0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                print "done fitting f3"
                print "f3 has parameters %f, %f, %f, %f)"%(f3.GetParameter(0), f3.GetParameter(1), f3.GetParameter(2), f3.GetParameter(3))
		Chi2List.append(f3.GetChisquare())
                f1.SetLineColor(kRed)
                f1.SetLineStyle(6)
                f2.SetLineColor(kBlack)
                f2.SetLineStyle(6)
                f3.SetLineColor(kViolet)
                f3.SetLineStyle(6)
                #fLow.SetParameters(1., f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2), f3.GetParameter(0), f3.GetParameter(1), f3.GetParameter(2), f3.GetParameter(3))
        	#stExcMEThist.Rebin(10)

            if i==3:
                print "got into the exc3 fitting loop!"
                f1_exc3.SetParameters(2.08e9, 4.28e-3, 6.7)
                for j in range(0, 20):
                    stExcMEThist.Fit("f1_exc3", "Q0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                stExcMEThist.Fit("f1_exc3", "E0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                print "done fitting f1_exc3"
		Chi2List.append(f1_exc3.GetChisquare())

                f2_exc3.SetParameters(1.7e12, 2.1, .16, 0.62)
                for j in range(0, 20):
                    stExcMEThist.Fit("f2_exc3", "Q0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                stExcMEThist.Fit("f2_exc3", "E0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                print "done fitting f2_exc3"
		Chi2List.append(f2_exc3.GetChisquare())

                f3_exc3.SetParameters(2e27, 4e5, -4e2, 3.6)
                for j in range(0, 20):
                    stExcMEThist.Fit("f3_exc3", "Q0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
		stExcMEThist.Fit("f3_exc3", "E0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                print "done fitting f3_exc3"
        	
		#stExcMEThist.Rebin(10)
		Chi2List.append(f3_exc3.GetChisquare())
                f1_exc3.SetLineColor(kGreen)
                f1_exc3.SetLineStyle(4)
                f2_exc3.SetLineColor(kCyan)
                f2_exc3.SetLineStyle(4)
                f3_exc3.SetLineColor(kMagenta)
                f3_exc3.SetLineStyle(4)
                #print "f1 has value " + str(f1.Eval(6500)) + " at ST=6500"
                #print "f2 has value " + str(f2.Eval(6500)) + " at ST=6500"
                #print "f3 has value " + str(f3.Eval(6500)) + " at ST=6500"
                #print "f1_exc3 has value " + str(f1_exc3.Eval(6500)) + " at ST=6500"
                #print "f2_exc3 has value " + str(f2_exc3.Eval(6500)) + " at ST=6500"
                #print "f3_exc3 has value " + str(f3_exc3.Eval(6500)) + " at ST=6500"
        if argv[4] == "useMHT" :
            if i==2:
                f1MHT.SetParameters(2.08e9, 4.28e-3, 6.7)
                for j in range(0, 20):
                    stExcMHThist.Fit("f1MHT", "L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2")  )
                f2MHT.SetParameters(1.7e12, 2.1, .16, 0.62)
                for j in range(0, 20):
                    stExcMHThist.Fit("f2MHT", "L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                f3MHT.SetParameters(2e27, 4e5, -4e2, 3.6)
                for j in range(0, 20):
                    stExcMHThist.Fit("f3MHT", "L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2")  )
                print "done fitting f1MHT_exc3"
                f1MHT.SetLineColor(kOrange+5)
                f1MHT.SetLineStyle(6)
                f2MHT.SetLineColor(kOrange+5)
                f2MHT.SetLineStyle(6)
                f3MHT.SetLineColor(kOrange+5)
                f3MHT.SetLineStyle(6)
            if i==3:
                f1MHT_exc3.SetParameters(2.08e9, 4.28e-3, 6.7)
                for j in range(0, 20):
                    stExcMHThist.Fit("f1MHT_exc3", "L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                f2MHT_exc3.SetParameters(1.7e12, 2.1, .16, 0.62)
                for j in range(0, 20):
                    stExcMHThist.Fit("f2MHT_exc3", "L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                print "done fitting f2MHT_exc3"
                f3MHT_exc3.SetParameters(2e27, 4e5, -4e2, 3.6)
                for j in range(0, 20):
                    stExcMHThist.Fit("f3MHT_exc3", "L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                print "done fitting f3MHT_exc3"
                f1MHT_exc3.SetLineColor(kOrange+5)
                f1MHT_exc3.SetLineStyle(6)
                f2MHT_exc3.SetLineColor(kOrange+5)
                f2MHT_exc3.SetLineStyle(6)
                f3MHT_exc3.SetLineColor(kOrange+5)
                f3MHT_exc3.SetLineStyle(6)
print "Printing list of chi2"
for j in range(0,len(Chi2List)):
	print Chi2List[j]
print "The minimum chi2 is %s %.3f" % (Chi2List.index(min(Chi2List)),min(Chi2List))
print "This should be the same chi2 %.3f" % (Funcs[Chi2List.index(min(Chi2List))].GetChisquare() )

for j in range(2,11):
    if argv[4] == "useMET" :
        stExcMEThist=PlotsDir.Get("stExc%02iHist"%j)
        stIncMEThist=PlotsDir.Get("stInc%02iHist"%j)

    if argv[4] == "useMHT" :
        stExcMHThist=PlotsDir.Get("stExc%02iHistMHT"%j)
        stIncMHThist=PlotsDir.Get("stInc%02iHistMHT"%j)

    STexcComparisons.append(TCanvas("stExc%02iCanvas"%j, "ST, N=%i"%j, 800, 800))
    STexcComparisons[j-2].cd()
    upperExcPads.append(TPad("Exc%02ipad"%j, "pad1", 0, 0.3, 1, 1.0))
    upperExcPads[j-2].SetBottomMargin(2)
    upperExcPads[j-2].Draw()
    upperExcPads[j-2].cd()
    upperExcPads[j-2].SetLogy()
    if argv[4] == "useMET" :
        stExcMEThist.GetXaxis().SetRangeUser(1000, 7000)
        stExcMEThist.GetXaxis().SetTitle("S_{T} (GeV)")
        stExcMEThist.GetYaxis().SetTitle("Events")
        stExcMEThist.SetMarkerColor(kBlack)
        stExcMEThist.SetMarkerStyle(8)
        stExcMEThist.SetMarkerSize(0.7)
        stExcMEThist.Draw("EP")
        lowerNormBin = stExcMEThist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("exc%i"%j)))
        upperNormBin = stExcMEThist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("exc%i"%j)))
        lowerNormEdge = stExcMEThist.GetXaxis().GetBinLowEdge(lowerNormBin)
        upperNormEdge = stExcMEThist.GetXaxis().GetBinLowEdge(upperNormBin)
        normBinTotal = 0;
        for normbin in range(lowerNormBin, upperNormBin):
            normBinTotal+=stExcMEThist.GetBinContent(normbin)
        #normfactor =  (normBinTotal/f3.Integral(lowerNormEdge, upperNormEdge))*stExcMEThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
        normfactor =  (normBinTotal/f2.Integral(lowerNormEdge, upperNormEdge))*stExcMEThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
        normfactor_exc3 =  (normBinTotal/f2_exc3.Integral(lowerNormEdge, upperNormEdge))*stExcMEThist.GetXaxis().GetBinWidth(upperNormBin)
        print "normalization factor is: %f" % normfactor
        f1Normalized = f1.Clone()
        f1Normalized.SetParameter(0, f1.GetParameter(0)*normfactor)
        print "f1Normalized has parameters" + " " + str(f1Normalized.GetParameter(0)) + " " + str(f1Normalized.GetParameter(1)) + " " + str(f1Normalized.GetParameter(2))
        f1Normalized.Draw("SAME")
        f3Normalized = f3.Clone()
        f3Normalized.SetParameter(0, f3.GetParameter(0)*normfactor)
        f3Normalized.Draw("SAME")
        f2Normalized = f2.Clone()
        f2Normalized.SetParameter(0, f2.GetParameter(0)*normfactor)
        f2Normalized.Draw("SAME")
        f1_exc3Normalized = f1_exc3.Clone()
        f1_exc3Normalized.SetParameter(0, f1_exc3.GetParameter(0)*normfactor_exc3)
        print "f1_exc3Normalized has parameters" + " " + str(f1_exc3Normalized.GetParameter(0)) + " " + str(f1_exc3Normalized.GetParameter(1)) + " " + str(f1_exc3Normalized.GetParameter(2))
        f1_exc3Normalized.Draw("SAME")
        f2_exc3Normalized = f2_exc3.Clone()
        f2_exc3Normalized.SetParameter(0, f2_exc3.GetParameter(0)*normfactor_exc3)
        f2_exc3Normalized.Draw("SAME")
        f3_exc3Normalized = f3_exc3.Clone()
        f3_exc3Normalized.SetParameter(0, f3_exc3.GetParameter(0)*normfactor_exc3)
        f3_exc3Normalized.Draw("SAME")
        fLow=TF1("fLow", "(2*([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x))))-([4]/([5]*x)**[6])", 1000, 7000)

        fLow.SetParameters(f2.GetParameter(0)*normfactor, f2.GetParameter(1), f2.GetParameter(2), f2.GetParameter(3), f1.GetParameter(0)*normfactor, f1.GetParameter(1), f1.GetParameter(2) )
        if j==2:
            print "fLow has parameters" + " " + str(fLow.GetParameter(0)) + " " + str(fLow.GetParameter(1)) + " " + str(fLow.GetParameter(2)) + " " + str(fLow.GetParameter(3)) + " " + str(fLow.GetParameter(4)) + " " + str(fLow.GetParameter(5))
            print "The chi2 for %s is %s" %( f1Normalized.GetName()	,stExcMEThist.Chisquare(f1Normalized,"R"))
	    print "The chi2 for %s is %s" %( f2Normalized.GetName()	,stExcMEThist.Chisquare(f2Normalized,"R"))
	    print "The chi2 for %s is %s" %( f3Normalized.GetName()	,stExcMEThist.Chisquare(f3Normalized,"R"))
	    print "The chi2 for %s is %s" %( f1_exc3Normalized.GetName()	,stExcMEThist.Chisquare(f1_exc3Normalized,"R"))
	    print "The chi2 for %s is %s" %( f2_exc3Normalized.GetName()	,stExcMEThist.Chisquare(f2_exc3Normalized,"R"))
	    print "The chi2 for %s is %s" %( f3_exc3Normalized.GetName()	,stExcMEThist.Chisquare(f3_exc3Normalized,"R"))

        fLow.SetLineColor(kBlue)
        fLow.SetLineStyle(6)
        fLow.Draw("SAME")

    if argv[4] == "useMHT" :
        stExcMHThist.GetXaxis().SetRangeUser(1000, 7000)
        stExcMHThist.GetXaxis().SetTitle("S_{T} (GeV)")
        stExcMHThist.GetYaxis().SetTitle("Events")
        stExcMHThist.SetMarkerColor(kBlack)
        stExcMHThist.SetMarkerStyle(8)
        stExcMHThist.SetMarkerSize(0.7)
        stExcMHThist.Draw("EP")
        lowerNormBin = stExcMHThist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("inc%i"%j)))
        upperNormBin = stExcMHThist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("inc%i"%j)))
        lowerNormEdge = stExcMHThist.GetXaxis().GetBinLowEdge(lowerNormBin)
        upperNormEdge = stExcMHThist.GetXaxis().GetBinLowEdge(upperNormBin)
        normBinTotal = 0;
        for normbin in range(lowerNormBin, upperNormBin):
            normBinTotal+=stExcMHThist.GetBinContent(normbin)
        normfactorMHT =  (normBinTotal/f1MHT.Integral(lowerNormEdge, upperNormEdge))*stExcMHThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
        normfactorMHT_exc3 =  (normBinTotal/f1MHT_exc3.Integral(lowerNormEdge, upperNormEdge))*stExcMHThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
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
        f1MHT_exc3Normalized = f1MHT_exc3.Clone()
        f1MHT_exc3Normalized.SetParameter(0, f1MHT_exc3.GetParameter(0)*normfactorMHT_exc3)
        f1MHT_exc3Normalized.Draw("SAME")
        f2MHT_exc3Normalized = f2MHT_exc3.Clone()
        f2MHT_exc3Normalized.SetParameter(0, f2MHT_exc3.GetParameter(0)*normfactorMHT_exc3)
        f2MHT_exc3Normalized.Draw("SAME")
        f3MHT_exc3Normalized = f3MHT_exc3.Clone()
        f3MHT_exc3Normalized.SetParameter(0, f3MHT_exc3.GetParameter(0)*normfactorMHT_exc3)
        f3MHT_exc3Normalized.Draw("SAME")


    legend = TLegend(0.35, 0.68, 0.58, 0.85, "Exclusive ST distributions for n=%i"%j, "brNDC")
    legend.SetTextSize(0.02);
    legend.SetLineColor(1);
    legend.SetLineWidth(1);
    legend.SetFillStyle(1001);
    legend.SetFillColor(10);
    if argv[4] == "useMET" :
        legend.AddEntry(stExcMEThist,"ST using MET","lp");
    if argv[4] == "useMHT" :
        legend.AddEntry(stExcMHThist,"ST using MHT","lp");
    legend.Draw()

    STexcComparisons[j-2].cd()
    lowerExcPads.append(TPad("Exc%02iratiopad"%j, "ratiopad1", 0, 0.04, 1, 0.3))
    lowerExcPads[j-2].SetTopMargin(5)
    lowerExcPads[j-2].SetBottomMargin(0.2)
    lowerExcPads[j-2].Draw()
    lowerExcPads[j-2].cd()

    if argv[4] == "useMET" :
        stExcMETRatio = stExcMEThist.Clone("stExcMETRatio")
        stExcMETRatio.GetXaxis().SetLabelSize(.08)
        stExcMETRatio.GetXaxis().SetTitle("")
        stExcMETRatio.GetYaxis().SetTitle("Ratio to n=2")
        stExcMETRatio.GetYaxis().SetLabelSize(.075)
        stExcMETRatio.GetYaxis().SetTitleSize(.1)
        stExcMETRatio.GetYaxis().SetTitleOffset(.3)
        stExcMETRatio.Divide(stExc2METhist)
        stExcMETRatio.SetTitle("")
        stExcMETRatio.SetMarkerColor(kBlack)
        stExcMETRatio.SetMarkerStyle(8)
        stExcMETRatio.SetMarkerSize(0.7)
        stExcMETRatio.Draw("EP")
        stExcMETRatio.SetStats(0)
    if argv[4] == "useMHT" :
        stExcMHTRatio = stExcMHThist.Clone("stExcMHTRatio")
        stExcMHTRatio.GetXaxis().SetLabelSize(.08)
        stExcMHTRatio.GetXaxis().SetTitle("")
        stExcMHTRatio.GetYaxis().SetTitle("Ratio to n=2")
        stExcMHTRatio.GetYaxis().SetLabelSize(.075)
        stExcMHTRatio.GetYaxis().SetTitleSize(.1)
        stExcMHTRatio.GetYaxis().SetTitleOffset(.3)
        stExcMHTRatio.Divide(stExc2MHThist)
        stExcMHTRatio.SetTitle("")
        stExcMHTRatio.SetMarkerColor(kBlack)
        stExcMHTRatio.SetMarkerStyle(8)
        stExcMHTRatio.SetMarkerSize(0.7)
        stExcMHTRatio.Draw("EP")
        stExcMHTRatio.SetStats(0)
    STexcComparisons[j-2].Write()

    STincComparisons.append(TCanvas("stInc%02iCanvas"%j, "ST, N>=%i"%j, 800, 800))
    STincComparisons[j-2].cd()
    upperIncPads.append(TPad("Inc%02ipad"%j, "pad1", 0, 0.3, 1, 1.0))
    upperIncPads[j-2].SetBottomMargin(2)
    upperIncPads[j-2].Draw()
    upperIncPads[j-2].cd()
    upperIncPads[j-2].SetLogy()
    if argv[4] == "useMET" :
        stIncMEThist.GetXaxis().SetRangeUser(1000, 7000)
        stIncMEThist.GetXaxis().SetTitle("S_{T} (GeV)")
        stIncMEThist.GetYaxis().SetTitle("Events")
        stIncMEThist.SetMarkerColor(kBlack)
        stIncMEThist.SetMarkerStyle(8)
        stIncMEThist.SetMarkerSize(0.7)
        stIncMEThist.Draw("EP")
        lowerNormBin = stIncMEThist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("inc%i"%j)))
        upperNormBin = stIncMEThist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("inc%i"%j)))
        lowerNormEdge = stIncMEThist.GetXaxis().GetBinLowEdge(lowerNormBin)
        upperNormEdge = stIncMEThist.GetXaxis().GetBinLowEdge(upperNormBin)
        normBinTotal = 0;
        for normbin in range(lowerNormBin, upperNormBin):
            normBinTotal+=stIncMEThist.GetBinContent(normbin)
        normfactorMET =  (normBinTotal/f2.Integral(lowerNormEdge, upperNormEdge))*stIncMEThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
        normfactorMET_exc3 =  (normBinTotal/f2_exc3.Integral(lowerNormEdge, upperNormEdge))*stIncMEThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
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
        f1_exc3Normalized = f1_exc3.Clone()
        f1_exc3Normalized.SetParameter(0, f1_exc3.GetParameter(0)*normfactorMET_exc3)
        f1_exc3Normalized.Draw("SAME")
        f2_exc3Normalized = f2_exc3.Clone()
        f2_exc3Normalized.SetParameter(0, f2_exc3.GetParameter(0)*normfactorMET_exc3)
        f2_exc3Normalized.Draw("SAME")
        f3_exc3Normalized = f3_exc3.Clone()
        f3_exc3Normalized.SetParameter(0, f3_exc3.GetParameter(0)*normfactorMET_exc3)
        f3_exc3Normalized.Draw("SAME")
        fLow=TF1("fLow", "(2*([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x))))-([4]/([5]*x)**[6])", 1000, 7000)

        fLow.SetParameters(f2.GetParameter(0)*normfactorMET, f2.GetParameter(1), f2.GetParameter(2), f2.GetParameter(3), f1.GetParameter(0)*normfactorMET, f1.GetParameter(1), f1.GetParameter(2) )
        #print "fLow has parameters" + " " + str(fLow.GetParameter(0)) + " " + str(fLow.GetParameter(1)) + " " + str(fLow.GetParameter(2)) + " " + str(fLow.GetParameter(3)) + " " + str(fLow.GetParameter(4)) + " " + str(fLow.GetParameter(5))
        fLow.SetLineColor(kBlue)
        fLow.SetLineStyle(6)
        fLow.Draw("SAME")
        outputForLimits = open("output/%s_Inclusive%i.txt"%(argv[2],j), "w")
        outputForLimits.write(" STMin    ::   Observed Data   ::   Expected Bkg   ::  Shape Unc  \n")
        for stmin in range(20, 75):
            observed=0
            startbin=stIncMEThist.GetXaxis().FindBin(float(stmin*100))
            for stbin in range (startbin, stIncMEThist.GetXaxis().GetNbins()):
                observed+=stIncMEThist.GetBinContent(stbin)
            expected = f2Normalized.Integral(stmin*100, 9999999)/100
            shapeUnc = abs(f1Normalized.Integral(stmin*100, 999999)/100-expected)/expected +1
                    #max(abs(f2Normalized.Integral(stmin*100, 999999)/100-expected),  #the 100's here are the bin width in GeV
                    #abs(f3Normalized.Integral(stmin*100, 999999)/100-expected)
                    #abs(f2_exc3Normalized.Integral(stmin*100, 999999)/100-expected),
                    #abs(f3_exc3Normalized.Integral(stmin*100, 999999)/100-expected)
                    #) /expected + 1
            if not ((j>5 and stmin<23) or (j>8 and stmin<25) or (j>10 and stmin<26)):
                outputForLimits.write("%i :: %i :: %f :: %f\n" % (stmin*100, observed, expected, shapeUnc))
        outputForLimits.close()

    if argv[4] == "useMHT" :
        stIncMHThist.GetXaxis().SetRangeUser(1000, 7000)
        stIncMHThist.GetXaxis().SetTitle("S_{T} (GeV)")
        stIncMHThist.GetYaxis().SetTitle("Events")
        stIncMHThist.SetMarkerColor(kBlack)
        stIncMHThist.SetMarkerStyle(8)
        stIncMHThist.SetMarkerSize(0.7)
        stIncMHThist.Draw("EP")
        lowerNormBin = stIncMHThist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("inc%i"%j)))
        upperNormBin = stIncMHThist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("inc%i"%j)))
        lowerNormEdge = stIncMHThist.GetXaxis().GetBinLowEdge(lowerNormBin)
        upperNormEdge = stIncMHThist.GetXaxis().GetBinLowEdge(upperNormBin)
        normBinTotal = 0;
        for normbin in range(lowerNormBin, upperNormBin):
            normBinTotal+=stIncMHThist.GetBinContent(normbin)
        normfactorMHT =  (normBinTotal/f1MHT.Integral(lowerNormEdge, upperNormEdge))*stIncMHThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
        normfactorMHT_exc3 =  (normBinTotal/f1MHT_exc3.Integral(lowerNormEdge, upperNormEdge))*stIncMHThist.GetXaxis().GetBinWidth(upperNormBin) # this assumes all the bins have the same width.
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
        f1MHT_exc3Normalized = f1MHT_exc3.Clone()
        f1MHT_exc3Normalized.SetParameter(0, f1MHT_exc3.GetParameter(0)*normfactorMHT_exc3)
        f1MHT_exc3Normalized.Draw("SAME")
        f2MHT_exc3Normalized = f2MHT_exc3.Clone()
        f2MHT_exc3Normalized.SetParameter(0, f2MHT_exc3.GetParameter(0)*normfactorMHT_exc3)
        f2MHT_exc3Normalized.Draw("SAME")
        f3MHT_exc3Normalized = f3MHT_exc3.Clone()
        f3MHT_exc3Normalized.SetParameter(0, f3MHT_exc3.GetParameter(0)*normfactorMHT_exc3)
        f3MHT_exc3Normalized.Draw("SAME")


        outputForLimits = open("output/%s_Inclusive%i.txt"%(argv[2],j), "w")
        outputForLimits.write(" STMin    ::   Observed Data   ::   Expected Bkg   ::  Shape Unc  \n")
        for stmin in range(20, 75):
            observed=0
            startbin=stIncMHThist.GetXaxis().FindBin(float(stmin*100))
            for stbin in range (startbin, stIncMHThist.GetXaxis().GetNbins()):
                observed+=stIncMHThist.GetBinContent(stbin)
            expected = f1MHTNormalized.Integral(stmin*100, 9999999)/100
            shapeUnc = max( abs(f2MHTNormalized.Integral(stmin*100, 999999)/100-expected),  abs(f3MHTNormalized.Integral(stmin*100, 999999)/100-expected)  )/expected + 1
            if not ((j>5 and stmin<23) or (j>8 and stmin<25) or (j>10 and stmin<26)):
                outputForLimits.write("%i :: %i :: %f :: %f\n" % (stmin*100, observed, expected, shapeUnc))
        outputForLimits.close()

    legend = TLegend(0.35, 0.68, 0.58, 0.85, "Inclusive ST distributions for n>=%i"%j, "brNDC")
    legend.SetTextSize(0.02);
    legend.SetLineColor(1);
    legend.SetLineStyle(1);
    legend.SetLineWidth(1);
    legend.SetFillStyle(1001);
    legend.SetFillColor(10);
    if argv[4] == "useMET" :
        legend.AddEntry(stExcMEThist,"ST using MET","lp");
    if argv[4] == "useMHT" :
        legend.AddEntry(stExcMHThist,"ST using MHT","lp");
    legend.Draw()

    STincComparisons[j-2].cd()
    lowerIncPads.append(TPad("Inc%02iratiopad"%j, "incratiopad1", 0, 0.04, 1, 0.3))
    lowerIncPads[j-2].SetTopMargin(5)
    lowerIncPads[j-2].SetBottomMargin(0.2)
    lowerIncPads[j-2].Draw()
    lowerIncPads[j-2].cd()

    if argv[4] == "useMET" :
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
        stIncMETRatio.SetMarkerColor(kBlack)
        stIncMETRatio.SetMarkerStyle(8)
        stIncMETRatio.SetMarkerSize(0.7)
        stIncMETRatio.Draw("EP")

    if argv[4] == "useMHT" :
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
        stIncMHTRatio.SetMarkerColor(kBlack)
        stIncMHTRatio.SetMarkerStyle(8)
        stIncMHTRatio.SetMarkerSize(0.7)
        stIncMHTRatio.Draw("EP")
    STincComparisons[j-2].Write()


#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]
