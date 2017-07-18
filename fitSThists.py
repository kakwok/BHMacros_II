# Script to fit ST distributions from the BHflatTuplizer, make them pretty, and spit out text files for the combine tool. John Hakala 12/1/2015
# this guy takes four arguments: the input BHflatTuple filename, the output filename, the name of the file defining the fit ranges, and either useMET or useMHT
# python fitSThists.py myBHflatTuple.root myOutputFile.root myFitNormRangesFile.txt useMHT
from ROOT import *
from fitAndNormRanges import *
from sys import argv
import CMS_lumi
import numpy as np

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
STcomparisons = {}   #  Dictionary of canvas
upperPads     = {}   #  Dictionary of upperPads in Canvas
lowerPads     = {}   #  Dictionary of lowerPads in Canvas
Chi2List     = []
f2_list      = {}     #  functions for fitted to ST n=2  
f3_list      = {}     #  functions for fitted to ST n=3
f4_list      = {}     #  functions for fitted to ST n=3
f23_list      = {}    #  functions for fitted to ST n<=3
f2_norm_list = {}     #  functions for fitted to ST n=2 after normalization
f3_norm_list = {}     #  functions for fitted to ST n=2 after normalization
f23_norm_list = {}    #  functions for fitted to ST n<=3 after normalization
PlotsFname = argv[1]
PlotsFile = TFile(argv[1])
PlotsDir = PlotsFile.Get("ST")
#PlotsDir = PlotsFile.Get("ST_tight")
OutFile = TFile("output/%s"%argv[2], "RECREATE")
chi2graphs   = {}
fitNormRanges = FitAndNormRange(argv[3])
fitNormRanges.showFitRanges()
fitNormRanges.showNormRanges()
rebin          = False 
WriteDataCards = False 
DrawUncertainty= True 
#f_outlier     = null


STup = 9000
def getratio(f1,f2):
    formula = "("+f1.GetExpFormula("p").Data() + ")/(" + f2.GetExpFormula("p").Data()+")-1"
    fname = f1.GetName()+"_over_+"+f2.GetName()
    f1overf2 = TF1(fname,formula,1000,STup)
    return f1overf2

#returns: f1 - (f2-f1), f1 is the best fit
def symmetrize(f1,f2):
    formula = "2*"+f1.GetExpFormula("p").Data() + "-" + f2.GetExpFormula("p").Data()
    fname = "2*"+f1.GetName()+"_minus_+"+f2.GetName()
    f1minusf2 = TF1(fname,formula,1000,STup)
    return f1minusf2
def symmetrizeFormula(f1,f2):
    formula = "2*"+f1.GetExpFormula("p").Data() + "-" + f2.GetExpFormula("p").Data()
    return formula


# Return the symmetrized funtion w.r.t. bestfit among functions in the range xlow, xup
def getSymmetrizedFunction(bestfit, functions, xlow, xup):
    bestfit_pos = 0
    fsym_pos    = -1
    #get the position of the best fit
    i=0
    fnames= []
    for f in functions:
        difference =  bestfit.Eval(xup)  -  f.Eval(xup) 
        if abs(difference) ==0:
            bestfit_pos = i
        fnames.append( f.GetName() )
        i+=1
    SymRange_xlow     = xlow    #range to symmetrize function
    SymRange_xup      = xup     #range to symmetrize function
    #[xlow,xlow+100, ... ,xup]
    fLowFormula = ""
    fUpFormula  = ""

    diffs             = []
    SignChanged       = False
    for x in np.arange(xlow,xup+100,100):
        last_diffs        = diffs
        diffs             = []
        abs_diffs         = []
        for f in functions:
            difference =  bestfit.Eval(x)  -  f.Eval(x) 
            diffs.append( difference )
            abs_diffs.append( abs(difference) )
        #print "At x= ",x
        #print "fnames = ",fnames
        #print "diffs  = ",diffs
        iMaxDiff = abs_diffs.index(max(abs_diffs)) 
        if (len(last_diffs)>0):
            if( diffs[iMaxDiff] * last_diffs[fsym_pos] < 0):
                SignChanged = True
            else:
                SignChanged = False

        SymRange_xup = x
        # Find do we need to change symmetrize function
        if((not fsym_pos == iMaxDiff) or (SignChanged)):
            if not(fsym_pos==-1):
                rangeString = "(x>="+str(SymRange_xlow)+" && x<"+str(SymRange_xup)+")*("
                print "Symmetrize function changed from %s to %s at %s" % (fnames[fsym_pos],fnames[iMaxDiff],x)
                if(not fLowFormula==""): fLowFormula += "+"
                if(not fUpFormula=="" ): fUpFormula  += "+"
               
                # make sure fLow picks up the smaller function 
                if( ( bestfit.Eval(SymRange_xlow) - functions[fsym_pos].Eval(SymRange_xlow))<0 ):
                    fLowFormula  += rangeString + symmetrizeFormula(bestfit, functions[fsym_pos]) +")"
                    fUpFormula   += rangeString + functions[fsym_pos].GetExpFormula("p").Data() +")"
                else:
                    fLowFormula  += rangeString +  functions[fsym_pos].GetExpFormula("p").Data()  +")"
                    fUpFormula   += rangeString + symmetrizeFormula(bestfit, functions[fsym_pos]) +")"
            #Mark the first pass
            fsym_pos      = iMaxDiff 
            SymRange_xlow = x
    #Fill the formula up to xup
    rangeString = "(x>="+str(SymRange_xlow)+" && x<"+str(SymRange_xup)+")*("
    if(not fLowFormula==""): fLowFormula += "+"
    if(not fUpFormula=="" ): fUpFormula  += "+"
    if(( bestfit.Eval(SymRange_xlow) - functions[fsym_pos].Eval(SymRange_xlow))<0):
        fLowFormula  += rangeString + symmetrizeFormula(bestfit, functions[fsym_pos]) +")"
        fUpFormula   += rangeString + functions[fsym_pos].GetExpFormula("p").Data() +")"
    else:
        fLowFormula  += rangeString +  functions[fsym_pos].GetExpFormula("p").Data()  +")"
        fUpFormula   += rangeString + symmetrizeFormula(bestfit, functions[fsym_pos]) +")"

    #print "Final fLow = ",fLowFormula
    #print "Final fUp  = ",fUpFormula
    
    flow = TF1("fLow_symmetrized",fLowFormula,xlow,xup)
    fup  = TF1("fUp_symmetrized",fUpFormula,xlow,xup)
    return (flow,fup)

# Return the TGraph bounded by fLow and fUp for drawing
def getFillGraph(fLow, fUp):
    g = TGraph()
    #for x in np.arange(fLow.GetXmin(),fLow.GetXmax(),130):
    for x in np.arange(fLow.GetXmin(),fLow.GetXmax(),100):
        g.SetPoint(g.GetN(), x, min(fLow.Eval(x),fUp.Eval(x)) )
    for x in np.arange(fUp.GetXmax(),fLow.GetXmin(),-100):
        g.SetPoint(g.GetN(), x, max(fUp.Eval(x),fLow.Eval(x)) )
    return g
 
# Return the TGraph bounded by fLow, fUp and normalized by fbest for drawing
def getRatioFillGraph(fLow, fUp, fbest):
    gUp = TGraph()
    gDown=TGraph()
    gbest=TGraph()
    gFill=TGraph()
    for x in np.arange(fLow.GetXmin(),fLow.GetXmax(),50):
        gFill.SetPoint(gFill.GetN(), x, (fLow.Eval(x)-fbest.Eval(x))/fbest.Eval(x))
        gDown.SetPoint(gDown.GetN(), x, (fLow.Eval(x)-fbest.Eval(x))/fbest.Eval(x))
        gbest.SetPoint(gbest.GetN(), x, 0)
    for x in np.arange(fUp.GetXmax(),fUp.GetXmin(),-50):
        gFill.SetPoint(gFill.GetN(), x, (fUp.Eval(x)-fbest.Eval(x))/fbest.Eval(x))
        gUp.SetPoint(  gUp.GetN(), x, (fUp.Eval(x)-fbest.Eval(x))/fbest.Eval(x))
    gDict={"gUp":gUp,"gDown":gDown,"gFill":gFill,"gbest":gbest}
    return gDict

# Return the mean value of bin content
def getMeanBinContent(hist):
    sumY=0
    nFilledBin = 0
    for ibin in range(1,hist.GetNbinsX()):
        if( hist.GetBinContent(ibin)!=0):
            sumY += hist.GetBinContent(ibin)
            nFilledBin +=1
    sumY = sumY /nFilledBin  
    return sumY
       
def getNormalizedFunction(f, hist, xlowbin, xupbin, refHist, xlowedge, xupedge, binwidth):
    debug = False;
    normBinTotal = 0;
    refHistTotal = 0;

    for normbin in range(xlowbin, xupbin):
        normBinTotal+=hist.GetBinContent(normbin)

    # this assumes all the bins have the same width and the fit well describe the normalization range, which is OK, because it's close to fit range
    normfactor =  (normBinTotal/f.Integral(xlowedge, xupedge))*binwidth 
    if debug:
        print " The normfactor for %s is %.3f  bin sum(numerator)=%s integral(denorminator) = %s" % ( f.GetName(), normfactor, normBinTotal, f.Integral(xlowedge,xupedge))
    fNormalized = f.Clone()
    fNormalized.SetRange(xlowedge, 14000)
    fNormalized.SetParameter(0, f.GetParameter(0)*normfactor)
    return fNormalized

def customfit(f, Sthist, norm):
    print "------------------------------------"
    print "Start fitting %s ...." % f.GetName()
    for j in range(0, 30):
        Sthist.Fit(f.GetName(), "Q0LR", "", fitNormRanges.getLowerFitBound(norm), fitNormRanges.getUpperFitBound(norm) )
    r = Sthist.Fit(f.GetName(), "0LR", "", fitNormRanges.getLowerFitBound(norm), fitNormRanges.getUpperFitBound(norm) )
    pars=[]
    for i in range(0,f.GetNpar()):
        pars.append(f.GetParameter(i))
    print "Done fitting %s, result = %s, %s has parameters :"%( f.GetName(),int(r), f.GetName()) , pars
    Chi2List.append(f.GetChisquare())

def ratioplot(fbest, sthist,xlow,xup):
    h = sthist.Clone("h_ratio")
 
###################################
##  This function takes stHist, get normalized fit function, symmetrize it and draw it on a canvas
##
##  stHist    = Inclusive/Exclusive ST histogram from flatuple  
##  j         = multiplicity of the histogram
## ExcOrInc   = "Exc" or "Inc"
## stRefHist  = Reference histogram (N=2/N=3/N<=3) for drawing ratio plot 
###################################
def FitAndDrawST(stHist,j,ExcOrInc,stRefHist,WriteCanvas):
    if rebin :
        stHist.Rebin(10)
    # For exclusive STs
    if ("Exc02" in stRefHist.GetName()):
	    canvasName = "st%s%02iCanvas_Exc02"%(ExcOrInc,j)
    if ("Exc03" in stRefHist.GetName()):
	    canvasName = "st%s%02iCanvas_Exc03"%(ExcOrInc,j)
    if ("Exc0203" in stRefHist.GetName()):
	    canvasName = "st%s%02iCanvas_Exc0203"%(ExcOrInc,j)
    UpperPadName    = "%s%02ipad"%(ExcOrInc,j)
    STcomparisons[canvasName] = TCanvas(canvasName, "ST, N=%i"%j, 700, 600)
    STcomparisons[canvasName].cd()
    upperPads[UpperPadName]        = (TPad("%s%02ipad"%(ExcOrInc,j), "pad1", 0, 0.3, 1, 1.0))
    upperPads[UpperPadName].SetBottomMargin(0.01)
    upperPads[UpperPadName].Draw()
    upperPads[UpperPadName].cd()
    upperPads[UpperPadName].SetLogy(1)
    stHist.SetTitle("")
    stHist.GetYaxis().SetTitle("Events/100 GeV")
    stHist.SetMarkerColor(kBlack)
    stHist.SetMarkerStyle(8)
    stHist.SetMarkerSize(0.7)
    stHist.Draw("EP")
    if (ExcOrInc=="Exc"):
        lowerNormBin  = stHist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("exc%i"%j)))
        upperNormBin  = stHist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("exc%i"%j)))
    if (ExcOrInc=="Inc"):
        lowerNormBin  = stHist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("inc%i"%j)))
        upperNormBin  = stHist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("inc%i"%j)))

    lowerNormEdge = stHist.GetXaxis().GetBinLowEdge(lowerNormBin)
    upperNormEdge = stHist.GetXaxis().GetBinLowEdge(upperNormBin)
    binwidth      = stHist.GetXaxis().GetBinWidth(upperNormBin)
    #stHist.GetXaxis().SetRangeUser(lowerNormEdge, STup)
    if (ExcOrInc=="Exc"):
        stHist.GetXaxis().SetRangeUser(fitNormRanges.getLowerPlotRange("exc%i"%j),fitNormRanges.getUpperPlotRange("exc%i"%j) )
    if (ExcOrInc=="Inc"):
        stHist.GetXaxis().SetRangeUser(fitNormRanges.getLowerPlotRange("inc%i"%j),fitNormRanges.getUpperPlotRange("inc%i"%j) )
    
    stHist.GetXaxis().SetLabelSize(0)
    stHist.SetMinimum(1e-1)

    print "in N=%i, upperNormEdge = %s, lowerNormEdge = %s" % (j,upperNormEdge, lowerNormEdge)
    chi2_list     = []
    #chi2 deviation list
    chi2_devlist  = []
    functions     = []
    # Normalize all fit functions
    for flist in AllFitList:
        for fname in flist:
	    fname_norm = fname +"_norm"
    	    fnorm      = getNormalizedFunction( flist[fname], stHist, lowerNormBin, upperNormBin, stRefHist,lowerNormEdge, upperNormEdge, binwidth)
            chi2 = stHist.Chisquare( fnorm ,"R")/ fnorm.GetNDF()
    	    chi2_list.append(chi2)
    	    chi2_devlist.append(abs(chi2-1))
    	    functions.append( fnorm )
    	    if(ExcOrInc=="Inc"):
                chi2graphs[fname].SetPoint( chi2graphs[fname].GetN(), j , chi2)
    	    if not DrawUncertainty:
                fnorm.Draw("SAME")
    
    #fbest    = f2Normalized
    #fLow     = getSymmetrizedFunction( fbest, functions, upperNormEdge, 14000)
    #fbest    = f2_norm_list["f2_norm"]
    fbest     = functions[ chi2_devlist.index( min(chi2_devlist) ) ]
    print "-----------------------------------------"
    print "In N=%i, fbest is chosen to be %s\n"%(j,fbest.GetName())
    fLow,fUp  = getSymmetrizedFunction( fbest, functions, 2500, 14000)
    fillGraph= getFillGraph( fLow, fUp )
    
    if DrawUncertainty:
        fillGraph.SetFillColorAlpha(kGray,0.35)
        fillGraph.SetLineColor(kBlue)
        fillGraph.Draw("sameF")
    	stHist.Draw("sameEP")
        fUp.SetLineColor(kBlue)
        fUp.SetLineStyle(1)
        fUp.SetLineWidth(1)
        fLow.SetLineWidth(1)
        fUp.Draw("SAME")
        fbest.SetLineColor(kBlue)
        fbest.SetLineStyle(1)
        fbest.SetLineWidth(2)
        fbest.Draw("SAME")
    else:
        # Draw a legend for all functions
        leg2 = TLegend(0.6,0.5, 0.85, 0.7,"", "brNDC")
        leg2.SetNColumns(2)
        leg2.SetBorderSize(0)
        for fnorm in functions:
            if(fnorm.GetName()==fbest.GetName()):
                leg2.AddEntry(fnorm,fnorm.GetName()+"(best-fit)","l")
            else:
                leg2.AddEntry(fnorm,fnorm.GetName(),"l")
        leg2.SetTextSize(0.03)
        leg2.Draw("SAME")
	fbest.Draw("SAME")
    fLow.SetLineColor(kBlue)
    fLow.Draw("SAME")

    chi2_list    = []
    chi2_devlist = []
    functions    = []

    legend = TLegend(0.6, 0.7, 0.8, 0.85,"", "brNDC")
    legend.SetTextSize(0.04);
    legend.SetLineWidth(1);
    legend.SetBorderSize(0);
    legend.SetFillStyle(1001);
    legend.SetFillColor(10);
    if(ExcOrInc=="Exc"):
        if("QCD" in PlotsFname):
            legend.AddEntry(stHist,"QCD: multiplicity =%i"%j,"ep");
        else:
            legend.AddEntry(stHist,"Data: multiplicity =%i"%j,"ep");
    if(ExcOrInc=="Inc"):
        if("QCD" in PlotsFname):
            legend.AddEntry(stHist,"QCD: multiplicity >=%i"%j,"ep");
        else:
            legend.AddEntry(stHist,"Data: multiplicity >=%i"%j,"ep");
    if DrawUncertainty:
	legend.AddEntry(fillGraph,"Background from fit","fl");
    legend.Draw()

    #CMS_lumi.CMS_lumi(STcomparisons[canvasName], iPeriod, iPos)

    #STcomparisons[canvasName].cd()
    #STcomparisons[canvasName].Update()
    #STcomparisons[canvasName].RedrawAxis()
    #frame = STcomparisons[canvasName].GetFrame()
    #frame.Draw()

    LowerPadName = "%s%02iratiopad"%(ExcOrInc,j)
    STcomparisons[canvasName].cd()
    lowerPads[LowerPadName] = (TPad(LowerPadName, "ratiopad1", 0, 0.04, 1, 0.3))
    lowerPads[LowerPadName].SetTopMargin(0.1)
    lowerPads[LowerPadName].SetBottomMargin(0.25)
    lowerPads[LowerPadName].Draw()
    lowerPads[LowerPadName].cd()

    # Save ST ratio
    stExcRatio = stHist.Clone("st%s%02i_RatioToExc%s"%(ExcOrInc,j,stRefHist.GetName()[6]))
    stExcRatio.Sumw2()
    if(ExcOrInc=="Exc"):
        stExcRatio.GetYaxis().SetTitle("Ratio of n=%i to n=%s"%(j,stRefHist.GetName()[6]))
    if(ExcOrInc=="Inc"):
        stExcRatio.GetYaxis().SetTitle("Ratio of n>=%i to n=%s"%(j,stRefHist.GetName()[6]))
    stExcRatio.Divide(stRefHist)
    print "%s has meanY =%s" %(stExcRatio.GetName(),getMeanBinContent(stExcRatio))
    stExcRatio.GetYaxis().SetRangeUser(0,getMeanBinContent(stExcRatio)*2)
    stExcRatio.Write()

    # Reuse StExcRatio to draw lower panel of data
    if DrawUncertainty:
        stExcRatio = stHist.Clone("st%s%02i_fitPanel"%(ExcOrInc,j))
        stExcRatio.Sumw2()
        stExcRatio.GetYaxis().SetTitle("(Data-Fit)/Fit")
        stExcRatio.Add(fbest,-1)    # Subtract best fit
        stExcRatio.Divide(fbest,1)  # Divide by best fit    
        stExcRatio.GetYaxis().SetRangeUser(-1,1)
    stExcRatio.GetXaxis().SetLabelSize(0.1)
    stExcRatio.GetXaxis().SetTitleSize(0.1)
    stExcRatio.GetXaxis().SetTitle("S_{T} (GeV)")
    stExcRatio.GetYaxis().SetLabelSize(.075)
    stExcRatio.GetYaxis().SetTitleSize(.1)
    stExcRatio.GetYaxis().SetTitleOffset(.3)
    stExcRatio.SetTitle("")
    stExcRatio.SetMarkerColor(kBlack)
    stExcRatio.SetMarkerStyle(8)
    stExcRatio.SetMarkerSize(0.7)
    stExcRatio.Draw("EP")
    stExcRatio.SetStats(0)

    #Draw Fit uncertainty
    if DrawUncertainty:
    	RatioFillGraphs = getRatioFillGraph( fLow, fUp, fbest )
    	RatioFillGraphs["gFill"].SetFillColorAlpha(kGray,0.35)
    	RatioFillGraphs["gFill"].Draw("sameF")
    	RatioFillGraphs["gUp"].SetLineColor(kBlue)
    	RatioFillGraphs["gUp"].Draw("sameC")
    	RatioFillGraphs["gDown"].SetLineColor(kBlue)
    	RatioFillGraphs["gDown"].Draw("sameC")
    	RatioFillGraphs["gbest"].SetLineColor(kBlack)
    	RatioFillGraphs["gbest"].Draw("sameC")
    	stExcRatio.Draw("sameEP")
    if(WriteCanvas):
        STcomparisons[canvasName].Write()

    if (WriteDataCards and ExcOrInc=="Inc"):
        outputForLimits = open("output/%s_Inclusive%i.txt"%(argv[2].replace(".root",""),j), "w")
        outputForLimits.write(" STMin    ::   Observed Data   ::   Expected Bkg   ::  Shape Unc  \n")
        for stmin in range(20, 90):
            observed=0
            startbin=stHist.GetXaxis().FindBin(float(stmin*100))
            for stbin in range (startbin, stHist.GetXaxis().GetNbins()):
                observed+=stHist.GetBinContent(stbin)
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



#################################################################################
#  Main
#################################################################################

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

f1_string = "[0]/([1]+0.001*x)**[2]"
f2_string = "([0]*(1+x*0.001)^[1])/((0.001*x)**([2]+[3]*TMath::Log(x*0.001)))"
f3_string = "[0]/([1] + [2]*x*0.001 + (0.001*x)**2)**[3]"
f4_string = "([0]*(1+0.001*x)^[1])/((0.001*x)**([2]*TMath::Log(x*0.001)))"
f5_string = "([0]*(1-0.001*x)^[1])/((0.001*x)**([2]+[3]*TMath::Log(x*0.001)))"


# Define dictionaries of functions for fitting different histograms
fnames = {"f1":f1_string,"f2":f2_string,"f3":f3_string,"f4":f4_string}
for fname in fnames:
	f2_list[fname+"_exc2"]     = TF1(fname+"_exc2",fnames[fname],1000,STup)
	f3_list[fname+"_exc3"]     = TF1(fname+"_exc3",fnames[fname],1000,STup)
	f4_list[fname+"_exc4"]     = TF1(fname+"_exc4",fnames[fname],1000,STup)

AllFitList=[f3_list,f4_list]
#AllFitList=[f3_list]

for flist in AllFitList:
    for fname in flist:
        chi2graphs[fname]=TGraph()
        if("f1" in fname):
            flist[fname].SetLineColor(kBlack)
            chi2graphs[fname].SetLineColor(kBlack)
        if("f2" in fname):
            flist[fname].SetLineColor(kRed)
            chi2graphs[fname].SetLineColor(kRed)
        if("f3" in fname):
            flist[fname].SetLineColor(kGreen)
            chi2graphs[fname].SetLineColor(kGreen)
        if("f4" in fname):
            flist[fname].SetLineColor(kMagenta)
            chi2graphs[fname].SetLineColor(kMagenta)

#fLow=TF1("fLow", "[0]*(2*[1]/([2]*x)**[3]-[4]/([5] + [6]*x + x**2)**[7])", 1000, STup)
#Funcs ={0:f1,1:f2,2:f3,3:f1_exc3,4:f2_exc3,5:f3_exc3,6:fLow}

##############################################
##  Fit the ST histogram, get the functions
##############################################
if(argv[4]=="useMET"):
    if( "ST_tight" in PlotsDir.GetName()):
        histname02 = ("stExc02Hist_tight")
        histname03 = ("stExc03Hist_tight")
        histname04 = ("stExc04Hist_tight")
    else:
        histname02 = ("stExc02Hist")
        histname03 = ("stExc03Hist")
        histname04 = ("stExc04Hist")
if(argv[4]=="useMHT"):
    if( "ST_tight" in PlotsDir.GetName()):
        histname02 = ("stExc02HistMHT_tight")
        histname03 = ("stExc03HistMHT_tight")
        histname04 = ("stExc04HistMHT_tight")
    else:
        histname02 = ("stExc02HistMHT")
        histname03 = ("stExc03HistMHT")
        histname04 = ("stExc04HistMHT")

stExc2Hist =PlotsDir.Get(histname02)
stExc3Hist =PlotsDir.Get(histname03)
stExc4Hist =PlotsDir.Get(histname04)
#stExc2or3Hist = stExc2Hist.Clone("stExc0203Hist")
#stExc2or3Hist.Add(stExc3Hist)

#####    Fit N=2   ###############
f2_list["f1_exc2"].SetParameters(8e6, 0.5, 9)
f2_list["f2_exc2"].SetParameters(2.4e6, -3, 4.7, 0.4)
f2_list["f3_exc2"].SetParameters(6e5, 0.4, -0.1, 4)
f2_list["f4_exc2"].SetParameters(1.5e9, -12,  -0.7)

#####    Fit N=3   #################
f3_list["f1_exc3"].SetParameters(2e13, 2, 14)
f3_list["f2_exc3"].SetParameters(2.4e6, -3, 4.7, 0.4)
f3_list["f3_exc3"].SetParameters(6e5, 0.4, -0.1, 4)
f3_list["f4_exc3"].SetParameters(1.5e9, -12,  -0.7)

#####    Fit N=4   #################
f4_list["f1_exc4"].SetParameters(8e6, 0.5, 9)
f4_list["f2_exc4"].SetParameters(2.4e6, -3, 4.7, 0.4)
f4_list["f3_exc4"].SetParameters(3e8, 0.3, -1, 4)
f4_list["f4_exc4"].SetParameters(1.5e10, -10,  0.3)

for flist in AllFitList:
    for fname in flist:
        if("exc2" in fname):
            refhist = stExc2Hist
            flist[fname].SetLineStyle(1)
            customfit( flist[fname], refhist,"exc2")
        if("exc3" in fname):
            refhist = stExc3Hist
            flist[fname].SetLineStyle(2)
            chi2graphs[fname].SetLineStyle(2)
            customfit( flist[fname], refhist,"exc3")
        if("exc4" in fname):
            refhist = stExc4Hist
            flist[fname].SetLineStyle(3)
            chi2graphs[fname].SetLineStyle(3)
            customfit( flist[fname], refhist,"exc4")

#    print "Printing list of chi2"
#for j in range(0,len(Chi2List)):
#	print Chi2List[j]
print "The minimum chi2 is %s %.3f" % (Chi2List.index(min(Chi2List)),min(Chi2List))

for j in range(2,12):
    if (argv[4]=="useMET"):
        if("ST_tight" in PlotsDir.GetName()):
           ExcHistName    =("stExc%02iHist_tight"%j)
           IncHistName    =("stInc%02iHist_tight"%j)
           Exc02HistName  =("stExc02Hist_tight")
           Exc03HistName  =("stExc03Hist_tight")
        else:
           ExcHistName    =("stExc%02iHist"%j)
           IncHistName    =("stInc%02iHist"%j)
           Exc02HistName  =("stExc02Hist")
           Exc03HistName  =("stExc03Hist")
          
    if (argv[4]=="useMHT"):
        if("ST_tight" in PlotsDir.GetName()):
           ExcHistName    =("stExc%02iHistMHT_tight"%j)
           IncHistName    =("stInc%02iHistMHT_tight"%j)
           Exc02HistName  =("stExc02HistMHT_tight")
           Exc03HistName  =("stExc03HistMHT_tight")
        else:
           ExcHistName    =("stExc%02iHistMHT"%j)
           IncHistName    =("stInc%02iHistMHT"%j)
           Exc02HistName  =("stExc02HistMHT")
           Exc03HistName  =("stExc03HistMHT")

    stExcHist =PlotsDir.Get(ExcHistName  )
    stIncHist =PlotsDir.Get(IncHistName  )
    stExc2Hist=PlotsDir.Get(Exc02HistName)
    stExc3Hist=PlotsDir.Get(Exc03HistName)
    
    if j==2:
        FitAndDrawST(stExcHist,j,"Exc",stExc3Hist,True)
    if j==3:
        FitAndDrawST(stExcHist,j,"Exc",stExc2Hist,True)
    if j==4:
        FitAndDrawST(stExcHist,j,"Exc",stExc3Hist,True)
        
    FitAndDrawST(stIncHist,j,"Inc",stExc2Hist,True)
    FitAndDrawST(stIncHist,j,"Inc",stExc3Hist,True)
    #FitAndDrawST(stExcHist,j,"Exc",stExc2or3Hist,False)
    #FitAndDrawST(stIncHist,j,"Inc",stExc2or3Hist,False)

c1= TCanvas("chi2graph","Chi2 vs N", 800,600)
OutFile.Append(c1)

#chi2graph_f5.SetLineColor(kYellow)
#chi2graph_f5_exc3.SetLineColor(kSpring+10)

leg = TLegend(0.7,0.7,0.9,0.9)
chi2graphs["f1_exc4"].Draw()
for gname in chi2graphs:
   print "now drawing chi2 graphs %s"%gname
   chi2graphs[gname].Draw("SAME")
   leg.AddEntry(chi2graphs[gname],gname,"L")
chi2graphs["f1_exc4"].SetTitle("Chi2/Ndof for different fit functions")
chi2graphs["f1_exc4"].GetXaxis().SetTitle("Inclusive Multiplicity")
chi2graphs["f1_exc4"].GetYaxis().SetTitle("Chi2/Ndof")
chi2graphs["f1_exc4"].GetYaxis().SetRangeUser(0,50)
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
