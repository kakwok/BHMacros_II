# Script to fit ST distributions from the BHflatTuplizer, make them pretty, and spit out text files for the combine tool. John Hakala 12/1/2015
# this guy takes four arguments: the input BHflatTuple filename, the output filename, the name of the file defining the fit ranges, and either useMET or useMHT
# python fitSThists.py myBHflatTuple.root myOutputFile.root myFitNormRangesFile.txt useMHT
from ROOT import *
from fitAndNormRanges import *
from sys import argv
import CMS_lumi

##################################################################
CMS_lumi.lumi_13TeV = "2.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 33
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4
H_ref = 600;
W_ref = 800;
W = W_ref
H  = H_ref
# references for T, B, L, R
T = 0.08*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.04*W_ref

##################################################################
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
chi2graph_f4 = TGraph()
chi2graph_f5 = TGraph()
chi2graph_f1_exc3 = TGraph()
chi2graph_f2_exc3 = TGraph()
chi2graph_f3_exc3 = TGraph()
chi2graph_f4_exc3 = TGraph()
chi2graph_f5_exc3 = TGraph()
fitNormRanges = FitAndNormRange(argv[3])
fitNormRanges.showFitRanges()
fitNormRanges.showNormRanges()
rebin         = False 
#f_outlier     = null


STup = 8000
def getratio(f1,f2):
    formula = "("+f1.GetExpFormula("p").Data() + ")/(" + f2.GetExpFormula("p").Data()+")-1"
    fname = f1.GetName()+"_over_+"+f2.GetName()
    f1overf2 = TF1(fname,formula,1000,STup)
    return f1overf2

def symmetrize(f1,f2):
    formula = "2*"+f1.GetExpFormula("p").Data() + "-" + f2.GetExpFormula("p").Data()
    fname = "2*"+f1.GetName()+"_minus_+"+f2.GetName()
    f1minusf2 = TF1(fname,formula,1000,STup)
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
   # f_outlier= functions[fsym_pos]
   # print "The index of fsym is ", fsym_pos
    print " %s is symmetrized with respect to %s" % ( functions[bestfit_pos].GetName(), functions[fsym_pos].GetName() )
    return symmetrize( functions[bestfit_pos], functions[fsym_pos])
        
def getNormalizedFuntion(f, hist, xlowbin, xupbin, xlowedge, xupedge, binwidth):
    normBinTotal = 0;
    for normbin in range(xlowbin, xupbin):
        normBinTotal+=hist.GetBinContent(normbin)  
    normfactor =  (normBinTotal/f.Integral(xlowedge, xupedge))*binwidth # this assumes all the bins have the same width.
    print " The normfactor for %s is %.3f    bin sum=%s    integral = %s" % ( f.GetName(), normfactor, normBinTotal, f.Integral(xlowedge,xupedge))
    fNormalized = f.Clone()
    fNormalized.SetRange(xlowedge, 14000)
    fNormalized.SetParameter(0, f.GetParameter(0)*normfactor)
    return fNormalized

def customfit(f, Sthist, norm):
#    print "Start fitting %s ...." % f.GetName()
    for j in range(0, 20):
        Sthist.Fit(f.GetName(), "Q0LR", "", fitNormRanges.getLowerFitBound(norm), fitNormRanges.getUpperFitBound(norm) )
    Sthist.Fit(f.GetName(), "Q0LR", "", fitNormRanges.getLowerFitBound(norm), fitNormRanges.getUpperFitBound(norm) )
    pars=[]
    for i in range(0,f.GetNpar()):
        pars.append(f.GetParameter(i))
    print "Done fitting %s, %s has parameters :"%( f.GetName(), f.GetName()) , pars
    Chi2List.append(f.GetChisquare())

def ratioplot(fbest, sthist,xlow,xup):
    h = sthist.Clone("h_ratio")
   

if argv[4] == "useMET" :
    # John fits
    #f1 = TF1("f1", "[0]/([1]+x)**[2]", 1000, STup)
    #f2 = TF1("f2", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, STup)
    #f3 = TF1("f3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, STup)
    #f4 = TF1("f4", "([0]*(1+x)^[1])/(x**([2]*TMath::Log(x)))", 1000, STup)
    #f5 = TF1("f5", "([0]*(1-x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, STup)
    #f1_exc3 = TF1("f1_exc3", "[0]/([1]+x)**[2]", 1000, STup)
    #f2_exc3 = TF1("f2_exc3", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, STup)
    #f3_exc3 = TF1("f3_exc3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, STup)
    #f4_exc3 = TF1("f4_exc3", "([0]*(1+x)^[1])/(x**([2]*TMath::Log(x)))", 1000, STup)
    #f5_exc3 = TF1("f5_exc3", "([0]*(1-x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, STup)

    f1 = TF1("f1", "[0]/([1]+0.001*x)**[2]", 1000, STup)
    f2 = TF1("f2", "([0]*(1+x*0.001)^[1])/((0.001*x)**([2]+[3]*TMath::Log(x*0.001)))", 1000, STup)
    f3 = TF1("f3", "[0]/([1] + [2]*x*0.001 + (0.001*x)**2)**[3]", 1000, STup)
    f4 = TF1("f4", "([0]*(1+0.001*x)^[1])/((0.001*x)**([2]*TMath::Log(x*0.001)))", 1000, STup)
    f5 = TF1("f5", "([0]*(1-0.001*x)^[1])/((0.001*x)**([2]+[3]*TMath::Log(x*0.001)))", 1000, STup)
    f1_exc3 = TF1("f1_exc3", "[0]/([1]+0.001*x)**[2]", 1000, STup)
    f2_exc3 = TF1("f2_exc3", "([0]*(1+x*0.001)^[1])/((0.001*x)**([2]+[3]*TMath::Log(x*0.001)))", 1000, STup)
    f3_exc3 = TF1("f3_exc3", "[0]/([1] + [2]*x*0.001 + (0.001*x)**2)**[3]", 1000, STup)
    f4_exc3 = TF1("f4_exc3", "([0]*(1+0.001*x)^[1])/((0.001*x)**([2]*TMath::Log(x*0.001)))", 1000, STup)
    f5_exc3 = TF1("f5_exc3", "([0]*(1-0.001*x)^[1])/((0.001*x)**([2]+[3]*TMath::Log(x*0.001)))", 1000, STup)

    fLow=TF1("fLow", "[0]*(2*[1]/([2]*x)**[3]-[4]/([5] + [6]*x + x**2)**[7])", 1000, STup)
    Funcs ={0:f1,1:f2,2:f3,3:f1_exc3,4:f2_exc3,5:f3_exc3,6:fLow}
for i in range(2,4):
    if argv[4] == "useMET" :
        stExcMEThist=PlotsDir.Get("stExc%02iHist"%i)
        stExc2METhist=PlotsDir.Get("stExc02Hist")
        stIncMEThist=PlotsDir.Get("stInc%02iHist"%i)
    if (i == 2 or i==3 ):
        if argv[4] == "useMET" :
            if i==2:
		# John's parameter
                #f1.SetParameters(2.08e9, 4.28e-3, 6.7)
                #f2.SetParameters(1.7e12, 2.1, .16, 0.62)
                #f3.SetParameters(2e27, 4e5, -4e2, 3.6)
                #f4.SetParameters(8e11, -0.86,  0.62)
                #f5.SetParameters(1.7e12, -1, -3, 0.62)
		# Tutanon's function parameters
                f1.SetParameters(8e6, 0.5, 9)
                f2.SetParameters(2.4e6, -3, 4.7, 0.4)
                f3.SetParameters(6e5, 0.4, -0.1, 4)
                f4.SetParameters(1.5e9, -12,  -0.7)
                f5.SetParameters(2e5, 0, 6, 0.7)

                customfit(f1,stExcMEThist,"exc2")
                customfit(f2,stExcMEThist,"exc2")
                customfit(f3,stExcMEThist,"exc2")
                customfit(f4,stExcMEThist,"exc2")
                #customfit(f5,stExcMEThist,"exc2")

                f1.SetLineColor(kRed)
                f1.SetLineStyle(6)
                f2.SetLineColor(kBlack)
                f2.SetLineStyle(6)
                f3.SetLineColor(kViolet)
                f3.SetLineStyle(6)
                f4.SetLineColor(kPink)
                f4.SetLineStyle(6)
                f5.SetLineColor(kYellow)
                f5.SetLineStyle(6)

            if i==3:
		# John's parameter
                #f1_exc3.SetParameters(2.08e9, 4.28e-3, 6.7)
                #f2_exc3.SetParameters(1.7e12, 2.1, .16, 0.62)
                #f3_exc3.SetParameters(2e27, 4e5, -4e2, 3.6)
                #f4_exc3.SetParameters(8e11, -0.86,  0.62)
                #f5_exc3.SetParameters(1.7e12, -1, -3, 0.62)
		# Tutanon's function parameter
                f1_exc3.SetParameters(8e6, 0.5, 9)
                f2_exc3.SetParameters(2.4e6, -3, 4.7, 0.4)
                f3_exc3.SetParameters(6e5, 0.4, -0.1, 4)
                f4_exc3.SetParameters(1.5e9, -12,  -0.7)
                f5_exc3.SetParameters(2e5, 0, 6, 0.7)


                customfit(f1_exc3,stExcMEThist,"exc3")
                customfit(f2_exc3,stExcMEThist,"exc3")
                customfit(f3_exc3,stExcMEThist,"exc3")
                customfit(f4_exc3,stExcMEThist,"exc3")
                #customfit(f5_exc3,stExcMEThist,"exc3")
                f1_exc3.SetLineColor(kGreen)
                f1_exc3.SetLineStyle(4)
                f2_exc3.SetLineColor(kCyan)
                f2_exc3.SetLineStyle(4)
                f3_exc3.SetLineColor(kMagenta)
                f3_exc3.SetLineStyle(4)
                f4_exc3.SetLineColor(kOrange)
                f4_exc3.SetLineStyle(4)
                f5_exc3.SetLineColor(kSpring+10)
                f5_exc3.SetLineStyle(4)

                #print "f1 has value " + str(f1.Eval(6500)) + " at ST=6500"
                #print "f2 has value " + str(f2.Eval(6500)) + " at ST=6500"
                #print "f3 has value " + str(f3.Eval(6500)) + " at ST=6500"
                #print "f1_exc3 has value " + str(f1_exc3.Eval(6500)) + " at ST=6500"
                #print "f2_exc3 has value " + str(f2_exc3.Eval(6500)) + " at ST=6500"
                #print "f3_exc3 has value " + str(f3_exc3.Eval(6500)) + " at ST=6500"
#    print "Printing list of chi2"
#for j in range(0,len(Chi2List)):
#	print Chi2List[j]
print "The minimum chi2 is %s %.3f" % (Chi2List.index(min(Chi2List)),min(Chi2List))
print "This should be the same chi2 %.3f" % (Funcs[Chi2List.index(min(Chi2List))].GetChisquare() )

for j in range(2,11):
    if argv[4] == "useMET" :
        stExcMEThist=PlotsDir.Get("stExc%02iHist"%j)
        stIncMEThist=PlotsDir.Get("stInc%02iHist"%j)
    
    if rebin :
        stExcMEThist.Rebin(10)
        stIncMEThist.Rebin(10)
    # For exclusive STs
    STexcComparisons.append(TCanvas("stExc%02iCanvas"%j, "ST, N=%i"%j, 800, 800))
    STexcComparisons[j-2].cd()
    upperExcPads.append(TPad("Exc%02ipad"%j, "pad1", 0, 0.3, 1, 1.0))
    upperExcPads[j-2].SetBottomMargin(2)
    upperExcPads[j-2].Draw()
    upperExcPads[j-2].cd()
    upperExcPads[j-2].SetLogy()
    #if argv[4] == "useMET" :
    if argv[4] == "useMET" and (j==2 or j==3):
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
        stExcMEThist.GetXaxis().SetRangeUser(lowerNormEdge, STup)
        #stExcMEThist.GetXaxis().SetRangeUser(1400, STup)
        stExcMEThist.SetMinimum(1e-1)

        print "in N=%i, upperNormEdge = %s, lowerNormEdge = %s" % (j,upperNormEdge, lowerNormEdge)
        f1Normalized = getNormalizedFuntion( f1, stExcMEThist,lowerNormBin, upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f1Normalized.Draw("SAME")
        f2Normalized = getNormalizedFuntion( f2, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f2Normalized.Draw("SAME")
        f3Normalized = getNormalizedFuntion( f3, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f3Normalized.Draw("SAME")
        f4Normalized = getNormalizedFuntion( f4, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f4Normalized.Draw("SAME")
        #f5Normalized = getNormalizedFuntion( f5, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #f5Normalized.Draw("SAME")
        f1_exc3Normalized = getNormalizedFuntion( f1_exc3, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f1_exc3Normalized.Draw("SAME")
        f2_exc3Normalized = getNormalizedFuntion( f2_exc3, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f2_exc3Normalized.Draw("SAME")
        f3_exc3Normalized = getNormalizedFuntion( f3_exc3, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f3_exc3Normalized.Draw("SAME")
        f4_exc3Normalized = getNormalizedFuntion( f4_exc3, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f4_exc3Normalized.Draw("SAME")
        #f5_exc3Normalized = getNormalizedFuntion( f5_exc3, stExcMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #f5_exc3Normalized.Draw("SAME")

        chi2_f1      = stExcMEThist.Chisquare(f1Normalized,"R")/ f1Normalized.GetNDF()
        chi2_f2      = stExcMEThist.Chisquare(f2Normalized,"R")/ f2Normalized.GetNDF()
        chi2_f3      = stExcMEThist.Chisquare(f3Normalized,"R")/ f3Normalized.GetNDF()
        chi2_f4      = stExcMEThist.Chisquare(f4Normalized,"R")/ f4Normalized.GetNDF()
        #chi2_f5      = stExcMEThist.Chisquare(f5Normalized,"R")/ f5Normalized.GetNDF()
        chi2_f1_exc3 = stExcMEThist.Chisquare(f1_exc3Normalized,"R")/f1_exc3Normalized.GetNDF()
        chi2_f2_exc3 = stExcMEThist.Chisquare(f2_exc3Normalized,"R")/f2_exc3Normalized.GetNDF()
        chi2_f3_exc3 = stExcMEThist.Chisquare(f3_exc3Normalized,"R")/f3_exc3Normalized.GetNDF()
        chi2_f4_exc3 = stExcMEThist.Chisquare(f4_exc3Normalized,"R")/f4_exc3Normalized.GetNDF()
        #chi2_f5_exc3 = stExcMEThist.Chisquare(f5_exc3Normalized,"R")/f5_exc3Normalized.GetNDF()
       
        #print "Printing info of Exclusive %i fittings:" % j  
        #print "The chi2/Ndof for %s is %s" %( f1Normalized.GetName()	        , chi2_f1    ) 
	#print "The chi2/Ndof for %s is %s" %( f2Normalized.GetName()	        , chi2_f2    )
	#print "The chi2/Ndof for %s is %s" %( f3Normalized.GetName()	        , chi2_f3    )
	#print "The chi2/Ndof for %s is %s" %( f1_exc3Normalized.GetName()	, chi2_f1_exc3)
	#print "The chi2/Ndof for %s is %s" %( f2_exc3Normalized.GetName()	, chi2_f2_exc3)
	#print "The chi2/Ndof for %s is %s" %( f3_exc3Normalized.GetName()	, chi2_f3_exc3)

        #chi2_list = [chi2_f1, chi2_f2, chi2_f3,chi2_f4,chi2_f5, chi2_f1_exc3, chi2_f2_exc3, chi2_f3_exc3,chi2_f4_exc3,chi2_f5_exc3]
        #functions = [f1Normalized,f2Normalized,f3Normalized,f4Normalized,f5Normalized,f1_exc3Normalized,f2_exc3Normalized,f3_exc3Normalized ,f4_exc3Normalized ,f5_exc3Normalized ]
        chi2_list = [chi2_f1, chi2_f2, chi2_f3,chi2_f4, chi2_f1_exc3, chi2_f2_exc3, chi2_f3_exc3,chi2_f4_exc3]
        functions = [f1Normalized,f2Normalized,f3Normalized,f4Normalized,f1_exc3Normalized,f2_exc3Normalized,f3_exc3Normalized ,f4_exc3Normalized ]
        fbest    = f2Normalized
        #fLow     = getSymmetrizedFuntion( fbest, functions, upperNormEdge, 14000)
        fLow     = getSymmetrizedFuntion( fbest, functions, 3000, 14000)
 
        fLow.SetLineColor(kBlue)
        fLow.SetLineStyle(6)
        fLow.Draw("SAME")

   
    legend = TLegend(0.35, 0.68, 0.8, 0.85, "Exclusive ST distributions for n=%i"%j, "brNDC")
    legend.SetTextSize(0.04);
    legend.SetLineWidth(1);
    legend.SetBorderSize(0);
    legend.SetFillStyle(1001);
    legend.SetFillColor(10);
    if argv[4] == "useMET" :
        #legend.AddEntry(stExcMEThist,"ST using MET","lp");
        legend.AddEntry(stExcMEThist,"data","lp");
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
        binwidth      = stIncMEThist.GetXaxis().GetBinWidth(upperNormBin)
        stIncMEThist.GetXaxis().SetRangeUser(lowerNormEdge, STup)
        stIncMEThist.SetMinimum(1e-1)
        
        print "in N=%i, upperNormEdge = %s, lowerNormEdge = %s" % (j,upperNormEdge, lowerNormEdge)

        print "\nPrinting info of Inclusive %i fittings:" % j  
        f1Normalized = getNormalizedFuntion( f1, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f1Normalized.Draw("SAME")
        f2Normalized = getNormalizedFuntion( f2, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f2Normalized.Draw("SAME")
        f3Normalized = getNormalizedFuntion( f3, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f3Normalized.Draw("SAME")
        f4Normalized = getNormalizedFuntion( f4, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f4Normalized.Draw("SAME")
        #f5Normalized = getNormalizedFuntion( f5, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #f5Normalized.Draw("SAME")
        f1_exc3Normalized = getNormalizedFuntion( f1_exc3, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f1_exc3Normalized.Draw("SAME")
        f2_exc3Normalized = getNormalizedFuntion( f2_exc3, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f2_exc3Normalized.Draw("SAME")
        f3_exc3Normalized = getNormalizedFuntion( f3_exc3, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f3_exc3Normalized.Draw("SAME")
        f4_exc3Normalized = getNormalizedFuntion( f4_exc3, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        f4_exc3Normalized.Draw("SAME")
        #f5_exc3Normalized = getNormalizedFuntion( f5_exc3, stIncMEThist, lowerNormBin,upperNormBin, lowerNormEdge, upperNormEdge, binwidth )
        #f5_exc3Normalized.Draw("SAME")

        chi2_f1      = stIncMEThist.Chisquare(f1Normalized,"R")/ f1Normalized.GetNDF()
        chi2_f2      = stIncMEThist.Chisquare(f2Normalized,"R")/ f2Normalized.GetNDF()
        chi2_f3      = stIncMEThist.Chisquare(f3Normalized,"R")/ f3Normalized.GetNDF()
        chi2_f4      = stIncMEThist.Chisquare(f4Normalized,"R")/ f4Normalized.GetNDF()
        #chi2_f5      = stIncMEThist.Chisquare(f5Normalized,"R")/ f5Normalized.GetNDF()
        chi2_f1_exc3 = stIncMEThist.Chisquare(f1_exc3Normalized,"R")/f1_exc3Normalized.GetNDF()
        chi2_f2_exc3 = stIncMEThist.Chisquare(f2_exc3Normalized,"R")/f2_exc3Normalized.GetNDF()
        chi2_f3_exc3 = stIncMEThist.Chisquare(f3_exc3Normalized,"R")/f3_exc3Normalized.GetNDF()
        chi2_f4_exc3 = stIncMEThist.Chisquare(f4_exc3Normalized,"R")/f4_exc3Normalized.GetNDF()
        #chi2_f5_exc3 = stIncMEThist.Chisquare(f5_exc3Normalized,"R")/f5_exc3Normalized.GetNDF()
 
        #print "The chi2/Ndof for %s is %s" %( f1Normalized.GetName()	        , chi2_f1    ) 
	#print "The chi2/Ndof for %s is %s" %( f2Normalized.GetName()	        , chi2_f2    )
	#print "The chi2/Ndof for %s is %s" %( f3Normalized.GetName()	        , chi2_f3    )
	#print "The chi2/Ndof for %s is %s" %( f1_exc3Normalized.GetName()	, chi2_f1_exc3)
	#print "The chi2/Ndof for %s is %s" %( f2_exc3Normalized.GetName()	, chi2_f2_exc3)
	#print "The chi2/Ndof for %s is %s" %( f3_exc3Normalized.GetName()	, chi2_f3_exc3)

        chi2graph_f1.SetPoint( chi2graph_f1.GetN(), j , chi2_f1)
        chi2graph_f2.SetPoint( chi2graph_f2.GetN(), j , chi2_f2)
        chi2graph_f3.SetPoint( chi2graph_f3.GetN(), j , chi2_f3)
        chi2graph_f4.SetPoint( chi2graph_f4.GetN(), j , chi2_f4)
        #chi2graph_f5.SetPoint( chi2graph_f5.GetN(), j , chi2_f5)
        chi2graph_f1_exc3.SetPoint( chi2graph_f1_exc3.GetN(), j , chi2_f1_exc3)
        chi2graph_f2_exc3.SetPoint( chi2graph_f2_exc3.GetN(), j , chi2_f2_exc3)
        chi2graph_f3_exc3.SetPoint( chi2graph_f3_exc3.GetN(), j , chi2_f3_exc3)
        chi2graph_f4_exc3.SetPoint( chi2graph_f4_exc3.GetN(), j , chi2_f4_exc3)
        #chi2graph_f5_exc3.SetPoint( chi2graph_f5_exc3.GetN(), j , chi2_f5_exc3)

        #chi2_list = [chi2_f1, chi2_f2, chi2_f3,chi2_f4,chi2_f5, chi2_f1_exc3, chi2_f2_exc3, chi2_f3_exc3,chi2_f4_exc3,chi2_f5_exc3]
        #functions = [f1Normalized,f2Normalized,f3Normalized,f4Normalized,f5Normalized,f1_exc3Normalized,f2_exc3Normalized,f3_exc3Normalized ,f4_exc3Normalized ,f5_exc3Normalized ]
        chi2_list = [chi2_f1, chi2_f2, chi2_f3,chi2_f4, chi2_f1_exc3, chi2_f2_exc3, chi2_f3_exc3,chi2_f4_exc3]
        chi2_devlist =[abs(chi2_f1-1),abs(chi2_f2-1),abs(chi2_f3-1),abs(chi2_f4-1),abs(chi2_f1_exc3-1),abs(chi2_f2_exc3-1),abs(chi2_f3_exc3-1),abs(chi2_f4_exc3-1)]
        functions = [f1Normalized,f2Normalized,f3Normalized,f4Normalized,f1_exc3Normalized,f2_exc3Normalized,f3_exc3Normalized ,f4_exc3Normalized]
        #fbest    = f2Normalized
        #fbest    = f2_exc3Normalized
        fbest     = functions[ chi2_devlist.index( min(chi2_devlist) ) ]
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
            expected = fbest.Integral(stmin*100, 14000)/binwidth
            #expected = fbest.Integral(stmin*100, 9999999)/100
                    #expected = f2Normalized.Integral(stmin*100, 9999999)/100
            shapeUnc = abs(fLow.Integral(stmin*100, 14000)/binwidth-expected)/expected +1
                    #shapeUnc = abs(fLow.Integral(stmin*100, 999999)/100-expected)/expected +1
                    #shapeUnc = abs(f1Normalized.Integral(stmin*100, 999999)/100-expected)/expected +1
                    #max(abs(f2Normalized.Integral(stmin*100, 999999)/100-expected),  #the 100's here are the bin width in GeV
                    #abs(f3Normalized.Integral(stmin*100, 999999)/100-expected)
                    #abs(f2_exc3Normalized.Integral(stmin*100, 999999)/100-expected),
                    #abs(f3_exc3Normalized.Integral(stmin*100, 999999)/100-expected)
                    #) /expected + 1
            if not ((j>5 and stmin<23) or (j>8 and stmin<25) or (j>10 and stmin<26)):
                outputForLimits.write("%i :: %i :: %f :: %f\n" % (stmin*100, observed, expected, shapeUnc))
        outputForLimits.close()

    legend = TLegend(0.35, 0.68, 0.8, 0.85, "Inclusive ST distributions for n>=%i"%j, "brNDC")
    legend.SetTextSize(0.04);
    legend.SetLineColor(1);
    legend.SetBorderSize(0);
    legend.SetLineStyle(1);
    legend.SetLineWidth(1);
    legend.SetFillStyle(1001);
    legend.SetFillColor(10);
    if argv[4] == "useMET" :
        #legend.AddEntry(stExcMEThist,"ST using MET","lp");
        legend.AddEntry(stExcMEThist,"data","lp");
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
        #stIncMETRatio.Divide(stExc2METhist)
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
chi2graph_f4.SetLineColor(kPink)
#chi2graph_f5.SetLineColor(kYellow)
chi2graph_f1_exc3.SetLineColor(kGreen)
chi2graph_f2_exc3.SetLineColor(kCyan)
chi2graph_f3_exc3.SetLineColor(kMagenta)
chi2graph_f4_exc3.SetLineColor(kOrange)
#chi2graph_f5_exc3.SetLineColor(kSpring+10)

chi2graph_f1.Draw()
chi2graph_f2.Draw("SAME")
chi2graph_f3.Draw("SAME")
chi2graph_f4.Draw("SAME")
#chi2graph_f5.Draw("SAME")
chi2graph_f1_exc3.Draw("SAME")
chi2graph_f2_exc3.Draw("SAME")
chi2graph_f3_exc3.Draw("SAME")
chi2graph_f4_exc3.Draw("SAME")
#chi2graph_f5_exc3.Draw("SAME")
chi2graph_f1.SetTitle("Chi2/Ndof for different fit functions")
chi2graph_f1.GetXaxis().SetTitle("Inclusive Multiplicity")
chi2graph_f1.GetYaxis().SetTitle("Chi2/Ndof")
chi2graph_f1.GetYaxis().SetRangeUser(0,50)
leg = TLegend(0.7,0.7,0.9,0.9)
leg.AddEntry(chi2graph_f1, "f1","L")
leg.AddEntry(chi2graph_f2, "f2","L")
leg.AddEntry(chi2graph_f3, "f3","L")
leg.AddEntry(chi2graph_f4, "f4","L")
#leg.AddEntry(chi2graph_f5, "f5","L")
leg.AddEntry(chi2graph_f1_exc3, "f1_exc3","L")
leg.AddEntry(chi2graph_f2_exc3, "f2_exc3","L")
leg.AddEntry(chi2graph_f3_exc3, "f3_exc3","L")
leg.AddEntry(chi2graph_f4_exc3, "f4_exc3","L")
#leg.AddEntry(chi2graph_f5_exc3, "f5_exc3","L")
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
