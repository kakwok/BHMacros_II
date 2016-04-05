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
chi2graph_f1 = TGraph()
chi2graph_f2 = TGraph()
chi2graph_f3 = TGraph()
chi2graph_f1_exc3 = TGraph()
chi2graph_f2_exc3 = TGraph()
chi2graph_f3_exc3 = TGraph()
fitNormRanges = FitAndNormRange(argv[3])
fitNormRanges.showFitRanges()
fitNormRanges.showNormRanges()

def symmetrize(f1,f2):
    formula = "2*"+f1.GetExpFormula("p").Data() + "-" + f2.GetExpFormula("p").Data()
    fname = "2*"+f1.GetName()+"_minus_+"+f2.GetName()
    f1minusf2 = TF1(fname,formula,1000,7000)
    return f1minusf2

# Return the symmetrized funtion w.r.t. bestfit among functions in the range xlow, xup
def getSymmetrizedFuntion(bestfit, functions, xlow, xup):
    bestfit_pos = 0
    fsym_pos = 0
    i=0
    Integrate_bestfit = bestfit.Integral(xlow,xup)
    diffs             = []
    abs_diffs         = []
    for f in functions:
        difference =  Integrate_bestfit  -  f.Integral(xlow,xup) 
        if abs(difference) ==0:
            bestfit_pos = i
        diffs.append( difference )
        abs_diffs.append( abs(difference) )
        i+=1
   # print "The index of best fit is ", bestfit_pos 
    print "Differences with best fit: ",diffs
   # print "Largest differece is ",max(abs_diffs)
    fsym_pos = abs_diffs.index(max(abs_diffs))  
   # print "The index of fsym is ", fsym_pos
    print " %s is symmetrized with respect to %s" % ( functions[bestfit_pos].GetName(), functions[fsym_pos].GetName() )
    return symmetrize( functions[bestfit_pos], functions[fsym_pos])
    #if( diffs[fsym_pos] > 0):
    #    print "Best fit is greater then the outlier function"
    #    return symmetrize( functions[bestfit_pos], functions[fsym_pos])
    #else:
    #    print "Best fit is smaller then outlier function"
    #    return symmetrize( functions[bestfit_pos], functions[fsym_pos])
        
def getNormalizedFuntion(f, hist, xlowbin, xupbin, xlowedge, xupedge, binwidth):
    normBinTotal = 0;
    for normbin in range(xlowbin, xupbin):
        normBinTotal+=hist.GetBinContent(normbin)  
    normfactor =  (normBinTotal/f.Integral(xlowedge, xupedge))*binwidth # this assumes all the bins have the same width.
    fNormalized = f.Clone()
    fNormalized.SetRange(xlowedge, 14000)
    fNormalized.SetParameter(0, f.GetParameter(0)*normfactor)
    return fNormalized

if argv[4] == "useMET" :
    f1 = TF1("f1", "[0]/([1]+x)**[2]", 1000, 7000)
    f2 = TF1("f2", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
    f3 = TF1("f3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
    f1_exc3 = TF1("f1_exc3", "[0]/([1]+x)**[2]", 1000, 7000)
    f2_exc3 = TF1("f2_exc3", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
    f3_exc3 = TF1("f3_exc3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
    fLow=TF1("fLow", "[0]*(2*[1]/([2]*x)**[3]-[4]/([5] + [6]*x + x**2)**[7])", 1000, 7000)
    Funcs ={0:f1,1:f2,2:f3,3:f1_exc3,4:f2_exc3,5:f3_exc3,6:fLow}
for i in range(2,4):
    if argv[4] == "useMET" :
        stExcMEThist=PlotsDir.Get("stExc%02iHist"%i)
        stExc2METhist=PlotsDir.Get("stExc02Hist")
        stIncMEThist=PlotsDir.Get("stInc%02iHist"%i)
    if (i == 2 or i==3 ):
        if argv[4] == "useMET" :
            if i==2:
                f1.SetParameters(2.08e9, 4.28e-3, 6.7)
                for j in range(0, 20):
                    stExcMEThist.Fit("f1", "Q0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                #stExcMEThist.Fit("f1", "E0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                stExcMEThist.Fit("f1", "Q0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                print "done fitting f1"
                print "f1 has parameters %f, %f, %f)"%(f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2))
		Chi2List.append(f1.GetChisquare())

                f2.SetParameters(1.7e12, 2.1, .16, 0.62)
                for j in range(0, 20):
                    stExcMEThist.Fit("f2", "Q0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                #stExcMEThist.Fit("f2", "E0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                stExcMEThist.Fit("f2", "Q0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                print "done fitting f2"
                print "f2 has parameters %f, %f, %f, %f)"%(f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2), f2.GetParameter(3))
		Chi2List.append(f2.GetChisquare())

                f3.SetParameters(2e27, 4e5, -4e2, 3.6)
                for j in range(0, 20):
                    stExcMEThist.Fit("f3", "Q0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                #stExcMEThist.Fit("f3", "E0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
                stExcMEThist.Fit("f3", "Q0L", "", fitNormRanges.getLowerFitBound("exc2"), fitNormRanges.getUpperFitBound("exc2") )
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
                #stExcMEThist.Fit("f1_exc3", "E0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                stExcMEThist.Fit("f1_exc3", "Q0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                print "done fitting f1_exc3"
                print "f1_exc3 has parameters %f, %f, %f)"%(f1_exc3.GetParameter(0), f1_exc3.GetParameter(1), f1_exc3.GetParameter(2))
		Chi2List.append(f1_exc3.GetChisquare())

                f2_exc3.SetParameters(1.7e12, 2.1, .16, 0.62)
                for j in range(0, 20):
                    stExcMEThist.Fit("f2_exc3", "Q0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                #stExcMEThist.Fit("f2_exc3", "E0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                stExcMEThist.Fit("f2_exc3", "Q0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                print "done fitting f2_exc3"
                print "f2_exc3 has parameters %f, %f, %f, %f)"%(f2_exc3.GetParameter(0), f2_exc3.GetParameter(1), f2_exc3.GetParameter(2), f2_exc3.GetParameter(3))
		Chi2List.append(f2_exc3.GetChisquare())

                f3_exc3.SetParameters(2e27, 4e5, -4e2, 3.6)
                for j in range(0, 20):
                    stExcMEThist.Fit("f3_exc3", "Q0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
		#stExcMEThist.Fit("f3_exc3", "E0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
		stExcMEThist.Fit("f3_exc3", "Q0L", "", fitNormRanges.getLowerFitBound("exc3"), fitNormRanges.getUpperFitBound("exc3") )
                print "done fitting f3_exc3"
                print "f3_exc3 has parameters %f, %f, %f, %f)"%(f3_exc3.GetParameter(0), f3_exc3.GetParameter(1), f3_exc3.GetParameter(2), f3_exc3.GetParameter(3))
        	
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
    print "Printing list of chi2"
for j in range(0,len(Chi2List)):
	print Chi2List[j]
print "The minimum chi2 is %s %.3f" % (Chi2List.index(min(Chi2List)),min(Chi2List))
print "This should be the same chi2 %.3f" % (Funcs[Chi2List.index(min(Chi2List))].GetChisquare() )

for j in range(2,11):
    if argv[4] == "useMET" :
        stExcMEThist=PlotsDir.Get("stExc%02iHist"%j)
        stIncMEThist=PlotsDir.Get("stInc%02iHist"%j)
    
    # For exclusive STs
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
        binwidth      = stExcMEThist.GetXaxis().GetBinWidth(upperNormBin)

        f1Normalized = getNormalizedFuntion( f1, stExcMEThist,lowerNormBin, upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f1Normalized has parameters" + " " + str(f1Normalized.GetParameter(0)) + " " + str(f1Normalized.GetParameter(1)) + " " + str(f1Normalized.GetParameter(2))
        f1Normalized.Draw("SAME")
        f2Normalized = getNormalizedFuntion( f2, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f2Normalized has parameters" + " " + str(f2Normalized.GetParameter(0)) + " " + str(f2Normalized.GetParameter(1)) + " " + str(f2Normalized.GetParameter(2))
        f2Normalized.Draw("SAME")
        f3Normalized = getNormalizedFuntion( f3, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f3Normalized has parameters" + " " + str(f3Normalized.GetParameter(0)) + " " + str(f3Normalized.GetParameter(1)) + " " + str(f3Normalized.GetParameter(2))
        f3Normalized.Draw("SAME")
        f1_exc3Normalized = getNormalizedFuntion( f1_exc3, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f1_exc3Normalized has parameters" + " " + str(f1_exc3Normalized.GetParameter(0)) + " " + str(f1_exc3Normalized.GetParameter(1)) + " " + str(f1_exc3Normalized.GetParameter(2))
        f1_exc3Normalized.Draw("SAME")
        f2_exc3Normalized = getNormalizedFuntion( f2_exc3, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f2_exc3Normalized has parameters" + " " + str(f2_exc3Normalized.GetParameter(0)) + " " + str(f2_exc3Normalized.GetParameter(1)) + " " + str(f2_exc3Normalized.GetParameter(2))
        f2_exc3Normalized.Draw("SAME")
        f3_exc3Normalized = getNormalizedFuntion( f3_exc3, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f3_exc3Normalized has parameters" + " " + str(f3_exc3Normalized.GetParameter(0)) + " " + str(f3_exc3Normalized.GetParameter(1)) + " " + str(f3_exc3Normalized.GetParameter(2))
        f3_exc3Normalized.Draw("SAME")

        chi2_f1      = stExcMEThist.Chisquare(f1Normalized,"R")/ f1Normalized.GetNDF()
        chi2_f2      = stExcMEThist.Chisquare(f2Normalized,"R")/ f2Normalized.GetNDF()
        chi2_f3      = stExcMEThist.Chisquare(f3Normalized,"R")/ f3Normalized.GetNDF()
        chi2_f1_exc3 = stExcMEThist.Chisquare(f1_exc3Normalized,"R")/f1_exc3Normalized.GetNDF()
        chi2_f2_exc3 = stExcMEThist.Chisquare(f2_exc3Normalized,"R")/f2_exc3Normalized.GetNDF()
        chi2_f3_exc3 = stExcMEThist.Chisquare(f3_exc3Normalized,"R")/f3_exc3Normalized.GetNDF()
       
        #print "Printing info of Exclusive %i fittings:" % j  
        #print "The chi2/Ndof for %s is %s" %( f1Normalized.GetName()	        , chi2_f1    ) 
	#print "The chi2/Ndof for %s is %s" %( f2Normalized.GetName()	        , chi2_f2    )
	#print "The chi2/Ndof for %s is %s" %( f3Normalized.GetName()	        , chi2_f3    )
	#print "The chi2/Ndof for %s is %s" %( f1_exc3Normalized.GetName()	, chi2_f1_exc3)
	#print "The chi2/Ndof for %s is %s" %( f2_exc3Normalized.GetName()	, chi2_f2_exc3)
	#print "The chi2/Ndof for %s is %s" %( f3_exc3Normalized.GetName()	, chi2_f3_exc3)

        chi2_list = [chi2_f1, chi2_f2, chi2_f3, chi2_f1_exc3, chi2_f2_exc3, chi2_f3_exc3]
        functions = [f1Normalized,f2Normalized,f3Normalized,f1_exc3Normalized,f2_exc3Normalized,f3_exc3Normalized ]
        fbest    = f2Normalized
        fLow     = getSymmetrizedFuntion( fbest, functions, upperNormEdge, 14000)
 
        fLow.SetLineColor(kBlue)
        fLow.SetLineStyle(6)
        fLow.Draw("SAME")

   
    legend = TLegend(0.35, 0.68, 0.7, 0.85, "Exclusive ST distributions for n=%i"%j, "brNDC")
    legend.SetTextSize(0.04);
    legend.SetLineColor(1);
    legend.SetLineWidth(1);
    legend.SetFillStyle(1001);
    legend.SetFillColor(10);
    if argv[4] == "useMET" :
        legend.AddEntry(stExcMEThist,"ST using MET","lp");
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
    STexcComparisons[j-2].Write()

    STincComparisons.append(TCanvas("stInc%02iCanvas"%j, "ST, N>=%i"%j, 800, 800))
    STincComparisons[j-2].cd()
    upperIncPads.append(TPad("Inc%02ipad"%j, "pad1", 0, 0.3, 1, 1.0))
    upperIncPads[j-2].SetBottomMargin(2)
    upperIncPads[j-2].Draw()
    upperIncPads[j-2].cd()
    upperIncPads[j-2].SetLogy()

    # For Inclusive STs
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

        f1Normalized = getNormalizedFuntion( f1, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f1Normalized has parameters" + " " + str(f1Normalized.GetParameter(0)) + " " + str(f1Normalized.GetParameter(1)) + " " + str(f1Normalized.GetParameter(2))
        f1Normalized.Draw("SAME")
        f2Normalized = getNormalizedFuntion( f2, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f2Normalized has parameters" + " " + str(f2Normalized.GetParameter(0)) + " " + str(f2Normalized.GetParameter(1)) + " " + str(f2Normalized.GetParameter(2))
        f2Normalized.Draw("SAME")
        f3Normalized = getNormalizedFuntion( f3, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f3Normalized has parameters" + " " + str(f3Normalized.GetParameter(0)) + " " + str(f3Normalized.GetParameter(1)) + " " + str(f3Normalized.GetParameter(2))
        f3Normalized.Draw("SAME")
        f1_exc3Normalized = getNormalizedFuntion( f1_exc3, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f1_exc3Normalized has parameters" + " " + str(f1_exc3Normalized.GetParameter(0)) + " " + str(f1_exc3Normalized.GetParameter(1)) + " " + str(f1_exc3Normalized.GetParameter(2))
        f1_exc3Normalized.Draw("SAME")
        f2_exc3Normalized = getNormalizedFuntion( f2_exc3, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f2_exc3Normalized has parameters" + " " + str(f2_exc3Normalized.GetParameter(0)) + " " + str(f2_exc3Normalized.GetParameter(1)) + " " + str(f2_exc3Normalized.GetParameter(2))
        f2_exc3Normalized.Draw("SAME")
        f3_exc3Normalized = getNormalizedFuntion( f3_exc3, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #print "f3_exc3Normalized has parameters" + " " + str(f3_exc3Normalized.GetParameter(0)) + " " + str(f3_exc3Normalized.GetParameter(1)) + " " + str(f3_exc3Normalized.GetParameter(2))
        f3_exc3Normalized.Draw("SAME")

        chi2_f1      = stIncMEThist.Chisquare(f1Normalized,"R")/ f1Normalized.GetNDF()
        chi2_f2      = stIncMEThist.Chisquare(f2Normalized,"R")/ f2Normalized.GetNDF()
        chi2_f3      = stIncMEThist.Chisquare(f3Normalized,"R")/ f3Normalized.GetNDF()
        chi2_f1_exc3 = stIncMEThist.Chisquare(f1_exc3Normalized,"R")/f1_exc3Normalized.GetNDF()
        chi2_f2_exc3 = stIncMEThist.Chisquare(f2_exc3Normalized,"R")/f2_exc3Normalized.GetNDF()
        chi2_f3_exc3 = stIncMEThist.Chisquare(f3_exc3Normalized,"R")/f3_exc3Normalized.GetNDF()
        print "Printing info of Inclusive %i fittings:" % j  
        print "The chi2/Ndof for %s is %s" %( f1Normalized.GetName()	        , chi2_f1    ) 
	print "The chi2/Ndof for %s is %s" %( f2Normalized.GetName()	        , chi2_f2    )
	print "The chi2/Ndof for %s is %s" %( f3Normalized.GetName()	        , chi2_f3    )
	print "The chi2/Ndof for %s is %s" %( f1_exc3Normalized.GetName()	, chi2_f1_exc3)
	print "The chi2/Ndof for %s is %s" %( f2_exc3Normalized.GetName()	, chi2_f2_exc3)
	print "The chi2/Ndof for %s is %s" %( f3_exc3Normalized.GetName()	, chi2_f3_exc3)

        chi2graph_f1.SetPoint( chi2graph_f1.GetN(), j , chi2_f1)
        chi2graph_f2.SetPoint( chi2graph_f2.GetN(), j , chi2_f2)
        chi2graph_f3.SetPoint( chi2graph_f3.GetN(), j , chi2_f3)
        chi2graph_f1_exc3.SetPoint( chi2graph_f1_exc3.GetN(), j , chi2_f1_exc3)
        chi2graph_f2_exc3.SetPoint( chi2graph_f2_exc3.GetN(), j , chi2_f2_exc3)
        chi2graph_f3_exc3.SetPoint( chi2graph_f3_exc3.GetN(), j , chi2_f3_exc3)

        chi2_list = [chi2_f1, chi2_f2, chi2_f3, chi2_f1_exc3, chi2_f2_exc3, chi2_f3_exc3]
        functions = [f1Normalized,f2Normalized,f3Normalized,f1_exc3Normalized,f2_exc3Normalized,f3_exc3Normalized ]
        fbest    = f2Normalized
        #fbest    = f2_exc3Normalized
        #fLow     = getSymmetrizedFuntion( fbest, functions, upperNormEdge, 14000)
        fLow     = getSymmetrizedFuntion( fbest, functions, 3000, 14000)

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
            expected = fbest.Integral(stmin*100, 9999999)/100
                    #expected = f2Normalized.Integral(stmin*100, 9999999)/100
            shapeUnc = abs(fLow.Integral(stmin*100, 999999)/100-expected)/expected +1
                    #shapeUnc = abs(f1Normalized.Integral(stmin*100, 999999)/100-expected)/expected +1
                    #max(abs(f2Normalized.Integral(stmin*100, 999999)/100-expected),  #the 100's here are the bin width in GeV
                    #abs(f3Normalized.Integral(stmin*100, 999999)/100-expected)
                    #abs(f2_exc3Normalized.Integral(stmin*100, 999999)/100-expected),
                    #abs(f3_exc3Normalized.Integral(stmin*100, 999999)/100-expected)
                    #) /expected + 1
            if not ((j>5 and stmin<23) or (j>8 and stmin<25) or (j>10 and stmin<26)):
                outputForLimits.write("%i :: %i :: %f :: %f\n" % (stmin*100, observed, expected, shapeUnc))
        outputForLimits.close()

    legend = TLegend(0.35, 0.68, 0.7, 0.85, "Inclusive ST distributions for n>=%i"%j, "brNDC")
    legend.SetTextSize(0.04);
    legend.SetLineColor(1);
    legend.SetLineStyle(1);
    legend.SetLineWidth(1);
    legend.SetFillStyle(1001);
    legend.SetFillColor(10);
    if argv[4] == "useMET" :
        legend.AddEntry(stExcMEThist,"ST using MET","lp");
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

    STincComparisons[j-2].Write()

c1= TCanvas("chi2graph","Chi2 vs N", 800,600)
c1.cd()
chi2graph_f1.SetLineColor(kRed)
chi2graph_f2.SetLineColor(kBlack)
chi2graph_f3.SetLineColor(kViolet)
chi2graph_f1_exc3.SetLineColor(kGreen)
chi2graph_f2_exc3.SetLineColor(kCyan)
chi2graph_f3_exc3.SetLineColor(kMagenta)

chi2graph_f1.Draw()
chi2graph_f2.Draw("SAME")
chi2graph_f3.Draw("SAME")
chi2graph_f1_exc3.Draw("SAME")
chi2graph_f2_exc3.Draw("SAME")
chi2graph_f3_exc3.Draw("SAME")
chi2graph_f1.SetTitle("Chi2/Ndof for different fit functions")
chi2graph_f1.GetXaxis().SetTitle("Inclusive Multiplicity")
chi2graph_f1.GetYaxis().SetTitle("Chi2/Ndof")
leg = TLegend(0.7,0.7,0.9,0.9)
leg.AddEntry(chi2graph_f1, "f1","L")
leg.AddEntry(chi2graph_f2, "f2","L")
leg.AddEntry(chi2graph_f3, "f3","L")
leg.AddEntry(chi2graph_f1_exc3, "f1_exc3","L")
leg.AddEntry(chi2graph_f2_exc3, "f2_exc3","L")
leg.AddEntry(chi2graph_f3_exc3, "f3_exc3","L")
leg.SetFillStyle(1001);
leg.SetFillColor(0);
leg.Draw()
c1.Write()

#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]
