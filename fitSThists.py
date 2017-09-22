# Script to fit ST distributions from the BHflatTuplizer, make them pretty, and spit out text files for the combine tool. John Hakala 12/1/2015
# this guy takes four arguments: the input BHflatTuple filename, the output filename, the name of the file defining the fit ranges, and either useMET or useMHT
# python fitSThists.py myBHflatTuple.root myOutputFile.root myFitNormRangesFile.txt useMHT
from ROOT import *
from fitAndNormRanges import *
from sys import argv
from math import sqrt
import CMS_lumi
import numpy as np
from tabulate import tabulate

##################################################################
CMS_lumi.lumi_13TeV = "2.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 13
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4
H_ref = 600
W_ref = 800
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
f_integrals = []
f_chi2 =[] 
chi2graphs_norm   = {}
chi2graphs_fit    = {}
fitNormRanges = FitAndNormRange(argv[3])
fitNormRanges.showFitRanges()
fitNormRanges.showNormRanges()
rebin          = False   # Rebin from 50GeV to 100GeV
WriteDataCards = False 
DrawUncertainty= False
#f_outlier     = null
fbestName = ""

chi2Table          = []
chi2Table_head     = ["Name","Exc3/Exc4","chi2","ndof","chi2/ndof","chi2","full ndof","chi2up/ndof","integral(up)","fitResult"]
chi2Table.append(chi2Table_head)

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
def getEnvelopeFunctions(bestfit, functions, xlow, xup, mode):
    bestfit_pos = 0
    fsym_pos    = -1
    fup_pos     = -1
    flow_pos    = -1
    #get the position of the best fit
    i=0
    fnames= []
    for f in functions:
        difference =  bestfit.Eval(xup)  -  f.Eval(xup) 
        if abs(difference) ==0:
            bestfit_pos = i
            pass
        fnames.append( f.GetName() )
        i+=1
        pass
    SymRange_xlow     = xlow    #range to symmetrize function
    SymRange_xup      = xup     #range to symmetrize function
    #[xlow,xlow+100, ... ,xup]
    fLowFormula = ""
    fUpFormula  = ""

    diffs             = []
    SignChanged       = False
    if (mode=="symmetrize"):
        for x in np.arange(xlow,xup+50,50):
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
            iMinDiff = diffs.index(min(diffs)) # for picking lowest functions
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
    elif (mode=="shade"):
        for x in np.arange(xlow,xup+50,50):
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
            iMaxDiff = diffs.index(min(diffs)) 
            iMinDiff = diffs.index(max(diffs)) 

            fUpSymRange_xup = x
            fLowSymRange_xup = x
            # Find do we need to change symmetrize function
            if((not fup_pos == iMaxDiff) or not flow_pos == iMinDiff):
                if not(fup_pos==-1):
                    fUpRangeString = "(x>="+str(fUpSymRange_xlow)+" && x<"+str(fUpSymRange_xup)+")*("
                    #print "Maximum function changed from %s to %s at %s" % (fnames[fup_pos],fnames[iMaxDiff],x)
                    if(not fUpFormula=="" ): fUpFormula  += "+"
                    fUpFormula   += fUpRangeString +  functions[fup_pos].GetExpFormula("p").Data() +")"
                if not(flow_pos==-1):
                    fLowRangeString = "(x>="+str(fLowSymRange_xlow)+" && x<"+str(fLowSymRange_xup)+")*("
                    #print "Min function changed from %s to %s at %s" % (fnames[flow_pos],fnames[iMinDiff],x)
                    if(not fLowFormula=="" ): fLowFormula  += "+"
                    fLowFormula   += fLowRangeString +  functions[flow_pos].GetExpFormula("p").Data() +")"
                #Mark the first pass
                fup_pos      = iMaxDiff 
                flow_pos     = iMinDiff 
                fUpSymRange_xlow = x
                fLowSymRange_xlow = x
        #Fill the formula up to xup
        fUpRangeString  = "(x>="+str(fUpSymRange_xlow)+" && x<"+str(fUpSymRange_xup)+")*("
        fLowRangeString = "(x>="+str(fLowSymRange_xlow)+" && x<"+str(fLowSymRange_xup)+")*("
        if(not fLowFormula==""): fLowFormula += "+"
        if(not fUpFormula=="" ): fUpFormula  += "+"
        fLowFormula  += fUpRangeString  +  functions[flow_pos].GetExpFormula("p").Data() +")"
        fUpFormula   += fLowRangeString +  functions[fup_pos].GetExpFormula("p").Data() +")"
    else:
        print "ERROR! please use \"shade\" or \"symmetrize\" for mode"

    #print "Final fLow = ",fLowFormula
    #print "Final fUp  = ",fUpFormula
    
    flow = TF1("fLow_symmetrized",fLowFormula,xlow,xup)
    fup  = TF1("fUp_symmetrized",fUpFormula,xlow,xup)
    return (flow,fup)

# Add normalization error to fLow and fUp, return as TGraph. normErr given in fraction 
def AddNormError(fLow,fUp,fbest, normErr):
    fLow_norm = TGraph()
    fUp_norm  = TGraph()
    fUp_norm.SetName( fUp.GetName() + "_norm")
    fLow_norm.SetName( fLow.GetName() + "_norm")
    for x in np.arange(fUp.GetXmin(),fUp.GetXmax(),50):
        shape_err = (fUp.Eval(x) - fbest.Eval(x))
        norm_err  = fbest.Eval(x)* normErr
        delta     = sqrt( norm_err**2 + shape_err**2)
        fUp_norm.SetPoint(fUp_norm.GetN(), x, fbest.Eval(x)+delta )

    for x in np.arange(fLow.GetXmin(),fLow.GetXmax(),50):
        shape_err = (fbest.Eval(x)-fLow.Eval(x) )
        norm_err  = fbest.Eval(x)* normErr
        delta     = sqrt( norm_err**2 + shape_err**2)
        fLow_norm.SetPoint(fLow_norm.GetN(), x, fbest.Eval(x)-delta )
    return fLow_norm, fUp_norm
    

# Return the TGraph bounded by fLow and fUp for drawing
def getFillGraph(fLow, fUp):
    g = TGraph()
    #for x in np.arange(fLow.GetXmin(),fLow.GetXmax(),130):
    for x in np.arange(fLow.GetXmin(),fLow.GetXmax(),50):
        g.SetPoint(g.GetN(), x, min(fLow.Eval(x),fUp.Eval(x)) )
    for x in np.arange(fUp.GetXmax(),fLow.GetXmin(),-50):
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
        pass
    for x in np.arange(fUp.GetXmax(),fUp.GetXmin(),-50):
        gFill.SetPoint(gFill.GetN(), x, (fUp.Eval(x)-fbest.Eval(x))/fbest.Eval(x))
        gUp.SetPoint(  gUp.GetN(), x, (fUp.Eval(x)-fbest.Eval(x))/fbest.Eval(x))
        pass
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
            pass
        pass
    sumY = sumY /nFilledBin  
    return sumY
       
def getNormalizedFunctionWithChi2(f, hist, ExcOrInc, j, STlow=0, STup=0):
    debug = False;
    histBinTotal = 0;
    normBinTotal = 0;
    Total = 0;

    if ("exc3" in f.GetName()): normHist = stExc3Hist
    if ("exc4" in f.GetName()): normHist = stExc4Hist

    #For exclusive multiplicity, calculate chi2 in Fitting region to judge quality of fit
    if (ExcOrInc=="Exc"):
        LowerNormBound    = float(fitNormRanges.getLowerNormBound("exc%i"%j))
        UpperNormBound    = float(fitNormRanges.getUpperNormBound("exc%i"%j))
        lowerNormBin      = hist.GetXaxis().FindBin(LowerNormBound)
        upperNormBin      = hist.GetXaxis().FindBin(UpperNormBound)
        lowerNormBin_ref  = normHist.GetXaxis().FindBin(LowerNormBound)
        upperNormBin_ref  = normHist.GetXaxis().FindBin(UpperNormBound)
        pass
    #For inclusive multiplicity, calculate chi2 in norm region for reference.
    if (ExcOrInc=="Inc"):
        LowerNormBound    = float(fitNormRanges.getLowerNormBound("inc%i"%j))
        UpperNormBound    = float(fitNormRanges.getUpperNormBound("inc%i"%j))
        lowerNormBin      = hist.GetXaxis().FindBin(LowerNormBound)
        upperNormBin      = hist.GetXaxis().FindBin(UpperNormBound)
        lowerNormBin_ref  = normHist.GetXaxis().FindBin(LowerNormBound)
        upperNormBin_ref  = normHist.GetXaxis().FindBin(UpperNormBound)
        pass

    for normbin in range(lowerNormBin, upperNormBin):
        histBinTotal+=hist.GetBinContent(normbin)*hist.GetBinWidth(normbin)
    for normbin in range(lowerNormBin_ref, upperNormBin_ref):
        normBinTotal+=normHist.GetBinContent(normbin)*normHist.GetBinWidth(normbin)

    #normfactor =  (normBinTotal/f.Integral(xlowedge, xupedge))*binwidth 
    normfactor =  histBinTotal/normBinTotal 
    if debug: print " The normfactor for %s is %.3f  | bin sum(numerator)=%s bin sum(denorminator) = %s" % ( f.GetName(), normfactor, histBinTotal, normBinTotal )
    fNormalized = f.Clone()
    if not (STlow==0 and STup==0):
        fNormalized.SetRange(STlow, STup)
    else:
        fNormalized.SetRange(LowerNormBound, UpperNormBound)
    fNormalized.SetParameter(0, f.GetParameter(0)*normfactor)

    chi2sum = 0
    for normbin in range(lowerNormBin, upperNormBin):
        x       = hist.GetBinCenter(normbin)
        y       = hist.GetBinContent(normbin) 
        errY    = hist.GetBinError(normbin) 
        #chi2sum += pow( (y - fNormalized.Eval(x) )/errY  ,2)
        #chi2me = chi2sum / fNormalized.GetNDF()
        #chi2sum += pow( (y - fNormalized.Eval(x) )/errY  ,2)
        pass
    #chi2me = chi2sum / fNormalized.GetNDF()
    #chi2pDOF  = hist.Chisquare( fNormalized ,"R")/ fNormalized.GetNDF()
    chi2pDOF  = hist.Chisquare( fNormalized ,"R")
    fNormalized.SetRange(LowerNormBound, 14000)
    #print "%s   chi2 = %.3f    chi2_me = %.3f" % (fNormalized.GetName(), chi2, chi2me)

    #return fNormalized ,chi2me
    return fNormalized ,chi2pDOF

def customfit(f, Sthist, ExcN):
    print "------------------------------------"
    print "Start fitting %s ...." % f.GetName()
    # Guide the fit with the amplitude of the refhist in the fit region
    f.SetParameter(0, Sthist.Integral( Sthist.FindBin(fitNormRanges.getLowerFitBound(ExcN)), Sthist.FindBin(fitNormRanges.getUpperFitBound(ExcN))))
    f.SetParLimits(0, 0, 1E10)
    for j in range(0, 20):
        Sthist.Fit(f.GetName(), "Q0LRB", "", fitNormRanges.getLowerFitBound(ExcN), fitNormRanges.getUpperFitBound(ExcN) )
    r = Sthist.Fit(f.GetName(), "0LRB", "", fitNormRanges.getLowerFitBound(ExcN), fitNormRanges.getUpperFitBound(ExcN) )
    fClone = f.Clone()
    pars=[]
    for i in range(0,f.GetNpar()): pars.append(f.GetParameter(i))
    Chi2List.append(f.GetChisquare())
    chi2pNDF = f.GetChisquare()/ f.GetNDF()
    #Calculate chi2 for the full range with clone of f
    fClone.SetName(f.GetName()+"_fullRange")
    #fClone.SetRange(fitNormRanges.getUpperFitBound(ExcN),13000)
    fClone.SetRange(4500,13000)
    Sthist.Fit( fClone.GetName(), "Q0LRB", "" , 4500, 13000)
    chi2full     = Sthist.Chisquare( fClone, "R") 
    UpperInt     = fClone.Integral( 5000.0, 7000.0)
    ndf_full=0
    for i in range(Sthist.FindBin(4500),Sthist.FindBin(13000)+1):
        if not(Sthist.GetBinContent(i)==0):  ndf_full+=1
    ndf_full -= f.GetNpar()
    chi2fullpNDF = chi2full / ndf_full
    
    chi2sum =0
    #chi2chklist=[]
    #for iBin in range( Sthist.FindBin( fitNormRanges.getLowerFitBound(ExcN)), Sthist.FindBin(fitNormRanges.getUpperFitBound(ExcN))+1):
    #    x       = Sthist.GetBinCenter(iBin)
    #    y       = Sthist.GetBinContent(iBin) 
    #    errY    = Sthist.GetBinError(iBin) 
    #    chi2sum += pow( (y - f.Eval(x) )/errY  ,2)
    #    chi2chklist.append({"ST":x,"chi2term":pow( (y - f.Eval(x) )/errY  ,2)})
 
    if("exc3" in f.GetName()):
        fname = f.GetName().replace("_exc3","")
        chi2graphs_fit[fname].SetPoint( chi2graphs_fit[fname].GetN(), 3, chi2pNDF)
        f_integrals.append( UpperInt)
        f_chi2.append( chi2full)
        chi2Table_row = [fname,"3", "%.3f"%(f.GetChisquare()),f.GetNDF(), "%.3f"%chi2pNDF,"%.3f"%chi2full, ndf_full, "%.3f"%chi2fullpNDF ,"%.3f"%UpperInt  ,int(r)]
        chi2Table.append(chi2Table_row)
        pass
        
    if("exc4" in f.GetName()):
        fname = f.GetName().replace("_exc4","")
        chi2graphs_fit[fname].SetPoint( chi2graphs_fit[fname].GetN(), 4 ,chi2pNDF)
        #chi2Table_row = [fname,"4", "%.3f"%(f.GetChisquare()), "%.3f"%f.GetNDF(), "%.3f"%chi2pNDF,int(r)]
        chi2Table_row = [fname,"4", "%.3f"%(f.GetChisquare()), f.GetNDF(), "%.3f"%chi2pNDF,"%.3f"%chi2full, ndf_full, "%.3f"%chi2fullpNDF ,"%.3f"%UpperInt  ,int(r)]
        chi2Table.append(chi2Table_row)
        pass
    print "Done fitting %s, result = %s, %s has parameters :"%( f.GetName(),int(r), f.GetName()) , pars
    print "Done fitting %s, %s has chi2perNDF = %.5f :"%(f.GetName(), f.GetName(), chi2pNDF)
    #for x in chi2chklist:
    #    print x["ST"],"%.3f"%x["chi2term"]

def pickBestFit( functions):
    fInts = []
    int2func = {}
    for f in functions:
        fInt = f.Integral(6000,7000)
        fName = f.GetName()
        fInts.append(fInt)
        int2func[fInt] = f
        pass
    fInts = np.sort(np.array(fInts))
    medianInt = np.median(fInts)
    meanInt = np.mean(fInts)
    aveInt = (medianInt+meanInt)/2.0
    bestIntDiff = 999999.9
    for integral in np.nditer(fInts):
        if abs(integral - aveInt) < bestIntDiff: 
            bestIntF = int2func[float(integral)]
            bestIntDiff = abs(integral - aveInt)
            pass
        pass
    return bestIntF.GetName()[:-5]
    #return functions[ chi2_devlist.index( min(chi2_devlist) ) ]
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
def NormAndDrawST(stHist,j,ExcOrInc,stRefHist,WriteCanvas):
    # For exclusive STs
    if ("Exc02" in stRefHist.GetName()):   canvasName = "st%s%02iCanvas_Exc02"%(ExcOrInc,j)
    if ("Exc03" in stRefHist.GetName()):   canvasName = "st%s%02iCanvas_Exc03"%(ExcOrInc,j)
    if ("Exc0203" in stRefHist.GetName()): canvasName = "st%s%02iCanvas_Exc0203"%(ExcOrInc,j)
    UpperPadName    = "%s%02ipad"%(ExcOrInc,j)
    STcomparisons[canvasName] = TCanvas(canvasName, "ST, N=%i"%j, 700, 600)
    STcomparisons[canvasName].cd()
    upperPads[UpperPadName]        = (TPad("%s%02ipad"%(ExcOrInc,j), "pad1", 0, 0.3, 1, 1.0))
    upperPads[UpperPadName].SetBottomMargin(0.01)
    upperPads[UpperPadName].Draw()
    upperPads[UpperPadName].cd()
    upperPads[UpperPadName].SetLogy(1)
    stHist.SetTitle("")
    stHist.GetYaxis().SetTitle("Events/%i GeV"%stHist.GetBinWidth(1))
    stHist.SetMarkerColor(kBlack)
    stHist.SetMarkerStyle(8)
    stHist.SetMarkerSize(0.7)
    #stHist.Sumw2(False)
    stHist.Draw("EP")
    if (ExcOrInc=="Exc"):
        lowerNormBin  = stHist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("exc%i"%j)))
        upperNormBin  = stHist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("exc%i"%j)))
        pass
    if (ExcOrInc=="Inc"):
        lowerNormBin  = stHist.GetXaxis().FindBin(float(fitNormRanges.getLowerNormBound("inc%i"%j)))
        upperNormBin  = stHist.GetXaxis().FindBin(float(fitNormRanges.getUpperNormBound("inc%i"%j)))
        pass

    lowerNormEdge = stHist.GetXaxis().GetBinLowEdge(lowerNormBin)
    upperNormEdge = stHist.GetXaxis().GetBinLowEdge(upperNormBin)
    binwidth      = stHist.GetXaxis().GetBinWidth(upperNormBin)
    #stHist.GetXaxis().SetRangeUser(lowerNormEdge, STup)
    if (ExcOrInc=="Exc"): stHist.GetXaxis().SetRangeUser(fitNormRanges.getLowerPlotRange("exc%i"%j),fitNormRanges.getUpperPlotRange("exc%i"%j) )
    if (ExcOrInc=="Inc"): stHist.GetXaxis().SetRangeUser(fitNormRanges.getLowerPlotRange("inc%i"%j),fitNormRanges.getUpperPlotRange("inc%i"%j) )
    
    stHist.GetXaxis().SetLabelSize(0)
    stHist.SetMinimum(1e-1)

    print "in N=%i, upperNormEdge = %s, lowerNormEdge = %s" % (j,upperNormEdge, lowerNormEdge)
    chi2_list     = []
    #chi2 deviation list
    chi2_devlist  = []
    functions     = []
    fbest = TF1()
    fbest.SetName("failedToFind_fbest")
    # Normalize all fit functions
    for flist in AllFitList:
        for fname in flist:
            #fname_norm = fname +"_norm"
            fnorm,chi2pDOF      = getNormalizedFunctionWithChi2( flist[fname], stHist, ExcOrInc, j )
            if("Exc04" in stHist.GetName() and fbestName in fname): fbest = fnorm
            elif(fbestName+"_exc3" in fname): fbest = fnorm
            chi2_list.append(chi2pDOF)
            chi2_devlist.append(abs(chi2pDOF-1))
            if(ExcOrInc=="Exc"):
                if((ExcOrInc+str(j)).lower() in fname):
                    functions.append( fnorm )
            else:
                functions.append( fnorm )
            if(ExcOrInc=="Inc"):
                chi2graphs_norm[fname].SetPoint( chi2graphs_norm[fname].GetN(), j , chi2pDOF)
            if not DrawUncertainty:
                if (ExcOrInc=="Inc"): fnorm.Draw("SAME")
                elif( (ExcOrInc+str(j)).lower() in fname): fnorm.Draw("SAME")          # draw only the fitted function for exclusive spectrums
                pass
            pass
        pass
    #fbest    = f2Normalized
    #fLow     = getSymmetrizedFunction( fbest, functions, upperNormEdge, 14000)
    #fbest    = f2_norm_list["f2_norm"]
    #fbest     = functions[ chi2_devlist.index( min(chi2_devlist) ) ]
    #fbest = pickBestFit( functions, chi2_devlist )
    print "-----------------------------------------"
    print "In N=%i, fbest is chosen to be %s\n"%(j,fbest.GetName())
    fLow,fUp  = getEnvelopeFunctions( fbest, functions, 2500, 7500, "shade")
    fillGraph= getFillGraph( fLow, fUp )
    fLow_norm,fUp_norm= AddNormError( fLow, fUp, fbest, 0.1 )
    
    if DrawUncertainty:
        fillGraph.SetFillColorAlpha(kGray,0.35)
        fillGraph.SetLineColor(kBlue)
        fillGraph.Draw("sameF")
        stHist.Draw("sameEP")
        fUp.SetLineColor(kBlue)
        fUp.SetLineStyle(1)
        fUp.SetLineWidth(1)
        fLow.SetLineColor(kCyan)
        fLow.SetLineWidth(1)
        fLow.Draw("SAME")
        fUp.Draw("SAME")
        fUp_norm.SetLineColor(kRed)
        fUp_norm.SetLineStyle(1)
        fUp_norm.SetLineWidth(1)
        fUp_norm.Draw("LSAME")
        fLow_norm.SetLineColor(kGreen)
        fLow_norm.SetLineStyle(1)
        fLow_norm.SetLineWidth(1)
        fLow_norm.Draw("LSAME")

        fbest.SetLineColor(kBlue)
        fbest.SetLineStyle(1)
        fbest.SetLineWidth(2)
        fbest.Draw("SAME")
        pass
    else:
        # Draw a legend for all functions
        leg2 = TLegend(0.5,0.5, 0.85, 0.7,"", "brNDC")
        leg2.SetNColumns(2)
        leg2.SetBorderSize(0)
        for fnorm in functions:
            if(fnorm.GetName()==fbest.GetName()): leg2.AddEntry(fnorm,fnorm.GetName()+"(best-fit)","l")
            else: leg2.AddEntry(fnorm,fnorm.GetName(),"l")
        leg2.SetTextSize(0.03)
        leg2.Draw("SAME")

	fbest.Draw("SAME")

    legend = TLegend(0.6, 0.7, 0.8, 0.85,"", "brNDC")
    legend.SetTextSize(0.04)
    legend.SetLineWidth(1)
    legend.SetBorderSize(0)
    legend.SetFillStyle(1001)
    legend.SetFillColor(10)
    if(ExcOrInc=="Exc"):
        if("qcd" in PlotsFname or "QCD" in PlotsFname): legend.AddEntry(stHist,"QCD: multiplicity =%i"%j,"ep")
        else: legend.AddEntry(stHist,"Data: multiplicity =%i"%j,"ep")
        pass
    if(ExcOrInc=="Inc"):
        if("qcd" in PlotsFname or "QCD" in PlotsFname): legend.AddEntry(stHist,"QCD: multiplicity >=%i"%j,"ep")
        else: legend.AddEntry(stHist,"Data: multiplicity >=%i"%j,"ep")
        pass
    if DrawUncertainty: legend.AddEntry(fillGraph,"Background from fit","fl")
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
    if(ExcOrInc=="Exc"): stExcRatio.GetYaxis().SetTitle("Ratio of n=%i to n=%s"%(j,stRefHist.GetName()[6]))
    if(ExcOrInc=="Inc"): stExcRatio.GetYaxis().SetTitle("Ratio of n>=%i to n=%s"%(j,stRefHist.GetName()[6]))
    stExcRatio.Divide(stRefHist)
    print "%s has meanY =%s" %(stExcRatio.GetName(),getMeanBinContent(stExcRatio))
    stExcRatio.GetYaxis().SetRangeUser(0,getMeanBinContent(stExcRatio)*2)
    stExcRatio.Write()

    #Draw pulls in the lower panel
    stExcRatio = stHist.Clone("st%s%02i_fitPanel"%(ExcOrInc,j))
    stExcRatio.Sumw2()
    stExcRatio.GetYaxis().SetTitle("(QCD-Fit)/Fit")
    stExcRatio.Add(fbest,-1)    # Subtract best fit
    stExcRatio.Divide(fbest,1)  # Divide by best fit    
    stExcRatio.GetYaxis().SetRangeUser(-1,1)
    #stExcRatio.GetYaxis().SetRangeUser(-0.2,0.2)
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
        pass
    else:
        fpulls = []
        for fnorm in functions:
            fdiff_formula = "("+fnorm.GetExpFormula("p").Data() + "-" + fbest.GetExpFormula("p").Data() +")/("+fbest.GetExpFormula("p").Data()+")"
            f_diff = TF1(fnorm.GetName()+"_pull_"+fbest.GetName(),fdiff_formula,1000,STup)
            #print "fnorm.Eval(4000)=%s, fbest.Eval(4000)=%s, fdiff.Eval(4000)=%s"%(fnorm.Eval(4000),fbest.Eval(4000),f_diff.Eval(4000))
            #print "f_diff name=%s fdiff.Eval(4000)=%s"%(f_diff.GetName(),f_diff.Eval(4000))
            f_diff.SetLineColor( fnorm.GetLineColor() )
            f_diff.SetLineStyle( fnorm.GetLineStyle() )
            f_diff.SetLineWidth( fnorm.GetLineWidth() )
            fpulls.append(f_diff)
            pass
        for f in fpulls: f.Draw("same")
        pass
        
    if(WriteCanvas): STcomparisons[canvasName].Write()

    if (WriteDataCards and ExcOrInc=="Inc"):
        outputForLimits = open("output/%s_Inclusive%i.txt"%(argv[2].replace(".root",""),j), "w")
        outputForLimits.write(" STMin    ::   Observed Data   ::   Expected Bkg   ::  Shape Unc  \n")
        for stmin in range(20, 90):
            observed=0
            startbin=stHist.GetXaxis().FindBin(float(stmin*100))
            for stbin in range (startbin, stHist.GetXaxis().GetNbins()): observed+=stHist.GetBinContent(stbin)
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
            if not ((j>5 and stmin<23) or (j>8 and stmin<25) or (j>10 and stmin<26)): outputForLimits.write("%i :: %i :: %f :: %f\n" % (stmin*100, observed, expected, shapeUnc))
        outputForLimits.close()
        pass
    pass


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

#f1_string       ="[0]/([1]+x/13000)**[2]"                      # Cannot fit with 2500-4000
f3_string       ="[0]/([1] + [2]*x*0.001 + (0.001*x)**2)**[3]" # Failed to fit 
CMSBH1_string   ="([0]*(1+x/13000)^[1])/((x/13000)**([2]*TMath::Log(x/13000)))"
CMSBH2_string   ="([0]*(1+x/13000)^[1])/((x/13000)**([2]+[3]*TMath::Log(x/13000)))"
ATLAS1_string   ="([0]*(1-(x/13000)^(1/3))^[1])/((x/13000)^[2])"
ATLAS2_string   ="([0]*(1-(x/13000)^(1/3))^[1])/((x/13000)^([2]+[3]*(TMath::Log(x/13000))^2))"
dijet1_string   ="([0]*(1-x/13000)^[1])/((x/13000)**([2]))" 
dijet2_string   ="([0]*(1-x/13000)^[1])/((x/13000)**([2]+[3]*TMath::Log(x/13000)))" 
#dijet3_string   ="([0]*(1-x/13000)^[1])/((x/13000)**([2]+[3]*TMath::Log(x/13000)+[4]*TMath::Log(x/13000)^2))"
UA21_string     ="[0]*(x/13000)^[1]*TMath::Exp(1)^([2]*(x/13000))"
UA22_string     ="[0]*(x/13000)^([1])*TMath::Exp(1)^([2]*(x/13000)+[3]*(x/13000)^2)"
#UA23_string     ="[0]*(x/13000)^([1])*TMath::Exp(1)^([2]*(x/13000)+[3]*(x/13000)^2+[4]*(x/13000)^3)"
dijetMod        ="([0]*(1-x/13000)^[1])*(1+[2]*(x/13000))/((x/13000)**([3]+[4]*TMath::Log(x/13000)))"

ATLASBH1_string="([0]*(1-x/13000)^[1])*((x/13000)**([2]*TMath::Log(x/13000)))"           # standard dijet without constant exponent 
ATLASBH2_string="([0]*(1-x/13000)^[1])*((1+x/13000)**([2]*TMath::Log(x/13000)))"         # replace x by 1+x
ATLASBH3_string="([0]*(1-x/13000)^[1])*((TMath::Exp(1))**([2]*(x/13000)^2))"             # replace x by exp
ATLASBH4_string="([0]*(1-(x/13000)^(1/3))^[1])*((x/13000)^([2]*(TMath::Log(x/13000))))"  # replace (1-x) by (1-x^1/3)
ATLASBH5_string="([0]*(1-(x/13000))^[1])*((x/13000)^([2]*(x/13000)))"                    # replace lnx by x in exponent
ATLASBH6_string="([0]*(1-x/13000)^[1])*((1+x/13000)**([2]*(x/13000)))"                   # replace x by 1+x, replace lnx by x in exponent
#UA21_string     ="[0]*(x/13000)^[1]*TMath::Exp(1)^([2]*(x/13000))"
#UA22_string     ="[0]*(x/13000)^([1])*TMath::Exp(1)^([2]*(x/13000)+[3]*(x/13000)^2)"
#UA23_string     ="[0]*(x/13000)^([1])*TMath::Exp(1)^([2]*(x/13000)+[3]*(x/13000)^2+[4]*(x/13000)^3)"
#dijetMod        ="([0]*(1-x/13000)^[1])*(1+[2]*(x/13000))/((x/13000)**([2]+[3]*TMath::Log(x/13000)))"


# Define dictionaries of functions for fitting different histograms
#fnames = {"f1":f1_string,"f2":f2_string,"f3":f3_string,"f4":f4_string,"f5":f5_string}
#fnames = {"f1":f1_string,"f2":f2_string,"f4":f4_string}
#fnames.update({"VVdijet2":VVdijet2_string,"VVdijet1":VVdijet1_string})
#fnames.update({"ATLAS1":ATLAS1_string,"ATLAS2":ATLAS2_string})
#fnames.update({"dijet2":dijet2_string,"dijet3":dijet3_string})
#fnames.update({"exPow":exPow_string})
fnames = {
#"CMSBH1":CMSBH1_string,
#"CMSBH2":CMSBH2_string,
"dijet1":dijet1_string,
"dijet2":dijet2_string,
#"dijet3":dijet3_string,
"ATLAS1":ATLAS1_string,
"ATLAS2":ATLAS2_string,
#"UA21":UA21_string,
"UA22":UA22_string,
#"UA23":UA23_string,
"dijetMod":dijetMod,
"ATLASBH1":ATLASBH1_string,
"ATLASBH2":ATLASBH2_string,
"ATLASBH3":ATLASBH3_string,
"ATLASBH4":ATLASBH4_string,
"ATLASBH5":ATLASBH5_string,
#"ATLASBH6":ATLASBH6_string
}

#fnames = {"VVdijet2":VVdijet2_string,"VVdijet1":VVdijet1_string}
#fnames.update({"ATLAS1":ATLAS1_string,"ATLAS2":ATLAS2_string})
#fnames = {"dijet2":dijet2_string,"dijet3":dijet3_string}
#fnames.update({"exPow":exPow_string})
for fname in fnames:
    f2_list[fname+"_exc2"]     = TF1(fname+"_exc2",fnames[fname],1000,STup)
    f3_list[fname+"_exc3"]     = TF1(fname+"_exc3",fnames[fname],1000,STup)
    f4_list[fname+"_exc4"]     = TF1(fname+"_exc4",fnames[fname],1000,STup)
    chi2graphs_fit[fname] =TGraph()

AllFitList=[f3_list,f4_list]
#AllFitList=[f3_list]

for flist in AllFitList:
    for fname in sorted(flist.iterkeys()):
        chi2graphs_norm[fname]=TGraph()
        if("f1" in fname):
            flist[fname].SetLineColor(kBlack)
            chi2graphs_norm[fname].SetLineColor(kBlack)
            chi2graphs_fit["f1"].SetLineColor(kBlack)
        if("CMSBH2" in fname):
            flist[fname].SetLineColor(kRed)
            chi2graphs_norm[fname].SetLineColor(kRed)
            chi2graphs_fit["CMSBH2"].SetLineColor(kRed)
        if("f3" in fname):
            flist[fname].SetLineColor(kGreen)
            chi2graphs_norm[fname].SetLineColor(kGreen)
            chi2graphs_fit["f3"].SetLineColor(kGreen)
        if("CMSBH1" in fname):
            flist[fname].SetLineColor(kRed)
            flist[fname].SetLineWidth(3)
            chi2graphs_norm[fname].SetLineColor(kRed)
            chi2graphs_fit["CMSBH1"].SetLineColor(kRed)
        if("dijet1" in fname):
            flist[fname].SetLineColor(kOrange)
            chi2graphs_norm[fname].SetLineColor(kOrange)
            chi2graphs_fit["dijet1"].SetLineColor(kOrange)
        if("dijet2" in fname):
            flist[fname].SetLineColor(kOrange)
            flist[fname].SetLineWidth(3)
            chi2graphs_norm[fname].SetLineColor(kOrange)
            chi2graphs_fit["dijet2"].SetLineColor(kOrange)
        if("dijet3" in fname):
            flist[fname].SetLineColor(kOrange)
            flist[fname].SetLineWidth(5)
            chi2graphs_norm[fname].SetLineColor(kOrange)
            chi2graphs_fit["dijet3"].SetLineColor(kOrange)
        if("UA21" in fname):
            flist[fname].SetLineColor(kViolet)
            chi2graphs_norm[fname].SetLineColor(kViolet)
            chi2graphs_fit["UA21"].SetLineColor(kViolet)
        if("UA22" in fname):
            flist[fname].SetLineColor(kViolet)
            flist[fname].SetLineWidth(3)
            chi2graphs_norm[fname].SetLineColor(kViolet)
            chi2graphs_fit["UA22"].SetLineColor(kViolet)
        if("UA23" in fname):
            flist[fname].SetLineColor(kViolet)
            flist[fname].SetLineWidth(5)
            chi2graphs_norm[fname].SetLineColor(kViolet)
            chi2graphs_fit["UA22"].SetLineColor(kViolet)
        if("ATLAS1" in fname):
            flist[fname].SetLineColor(kCyan)
            chi2graphs_norm[fname].SetLineColor(kCyan)
            chi2graphs_fit["ATLAS1"].SetLineColor(kCyan)
        if("ATLAS2" in fname):
            flist[fname].SetLineColor(kCyan)
            flist[fname].SetLineWidth(3)
            chi2graphs_norm[fname].SetLineColor(kCyan)
            chi2graphs_fit["ATLAS2"].SetLineColor(kCyan)
        if("dijetMod" in fname):
            flist[fname].SetLineColor(kGreen)
            chi2graphs_norm[fname].SetLineColor(kGreen)
            chi2graphs_fit["dijetMod"].SetLineColor(kGreen)
        if("ATLASBH" in fname):
            i = int(fname.strip("ATLASTBH")[0])
            flist[fname].SetLineColor(kBlue-3+i)
            chi2graphs_norm[fname].SetLineColor(kBlue-3+i)
            chi2graphs_fit[fname[:-5]].SetLineColor(kBlue-3+i)



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
        pass
    else:
        histname02 = ("stExc02Hist")
        histname03 = ("stExc03Hist")
        histname04 = ("stExc04Hist")
        pass
    pass
if(argv[4]=="useMHT"):
    if( "ST_tight" in PlotsDir.GetName()):
        histname02 = ("stExc02HistMHT_tight")
        histname03 = ("stExc03HistMHT_tight")
        histname04 = ("stExc04HistMHT_tight")
        pass
    else:
        histname02 = ("stExc02HistMHT")
        histname03 = ("stExc03HistMHT")
        histname04 = ("stExc04HistMHT")
        pass
    pass

stExc2Hist =PlotsDir.Get(histname02)
stExc3Hist =PlotsDir.Get(histname03)
stExc4Hist =PlotsDir.Get(histname04)
#stExc2or3Hist = stExc2Hist.Clone("stExc0203Hist")
#stExc2or3Hist.Add(stExc3Hist)

#####    Fit N=3   #################
#if("f1_exc3" in f3_list): f3_list["f1_exc3"].SetParameters(2e13, 0.3, 1)
#if("f2_exc3" in f3_list): f3_list["f2_exc3"].SetParameters(2.4e6, -27, 1.9, -0.48)
#if("f3_exc3" in f3_list): f3_list["f3_exc3"].SetParameters(6e5, 0.4, -0.1, 4)
#if("f4_exc3" in f3_list): f3_list["f4_exc3"].SetParameters(1.5e9, -31,  -0.8)
#if("f5_exc3" in f3_list): f3_list["f5_exc3"].SetParameters(37, 9 ,3,  -0.9)
#
######    Fit N=4   #################
#if("f1_exc4" in f4_list): f4_list["f1_exc4"].SetParameters(8e6, 0.3, 1)
#if("f2_exc4" in f4_list): f4_list["f2_exc4"].SetParameters(2.4e6, -27, 1.9, -0.48)
#if("f3_exc4" in f4_list): f4_list["f3_exc4"].SetParameters(3e8, 0.3, -1, 4)
#if("f4_exc4" in f4_list): f4_list["f4_exc4"].SetParameters(1.5e10, -31,  -0.8)
#if("f5_exc4" in f4_list): f4_list["f5_exc4"].SetParameters(37, 9,3,  -0.9)

for flist in AllFitList:
    for fname in sorted(flist.iterkeys()):
        if("f1"       in fname):  flist[fname].SetParameters(1, 0.3,     1)
        if("CMSBH2"   in fname):  flist[fname].SetParameters(1, -27,   1.9,  -0.48)
        if("f3"       in fname):  flist[fname].SetParameters(1, 0.4,  -0.1,      4)
        if("CMSBH1"   in fname):  flist[fname].SetParameters(1, -31,  -0.8)
        if("dijet1"   in fname):  flist[fname].SetParameters(1,   9,     3,   -0.9)
        if("dijet2"   in fname):  flist[fname].SetParameters(1,  10,     3,   -0.7)
        if("dijet3"   in fname):  flist[fname].SetParameters(1,  15,     3,   2)
        #if("dijet3"   in fname):  flist[fname].SetParameters(1,  10,     3,   -0.2,   1)
        if("ATLAS1"   in fname):  flist[fname].SetParameters(1,   1,     1)
        if("ATLAS2"   in fname):  flist[fname].SetParameters(1,   1,     1,      1)
        if("UA21"     in fname):  flist[fname].SetParameters(1,   1,     -1)
        if("UA22"     in fname):  flist[fname].SetParameters(1,   -1,     -1,      1)
        if("UA23"     in fname):  flist[fname].SetParameters(1,   -1,     -1,      -1,   1)
        if("dijetMod_exc3" in fname):  flist[fname].SetParameters(1,   13,     0,     -1.2, -2)
        if("dijetMod_exc4" in fname):  flist[fname].SetParameters(1,   14.6,     0,      2.3, -0.5)

        if("exc2" in fname):
            refhist = stExc2Hist
            flist[fname].SetLineStyle(1)
            customfit( flist[fname], refhist,"exc2")
            pass
        if("exc3" in fname):
            refhist = stExc3Hist
            flist[fname].SetLineStyle(2)
            chi2graphs_norm[fname].SetLineStyle(2)
            customfit( flist[fname], refhist,"exc3")
            pass
        if("exc4" in fname):
            refhist = stExc4Hist
            flist[fname].SetLineStyle(3)
            chi2graphs_norm[fname].SetLineStyle(3)
            customfit( flist[fname], refhist,"exc4")
            pass
        pass
    pass

funcs4best     = []
# Normalize all fit functions
for flist in AllFitList:
    for fname in flist:
        if("exc3" in fname): funcs4best.append(flist[fname])
        pass
    pass
fbestName = pickBestFit(funcs4best)
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
            Exc04HistName  =("stExc04Hist_tight")
            pass
        else:
            ExcHistName    =("stExc%02iHist"%j)
            IncHistName    =("stInc%02iHist"%j)
            Exc02HistName  =("stExc02Hist")
            Exc03HistName  =("stExc03Hist")
            Exc04HistName  =("stExc04Hist")
            pass
        pass
          
    if (argv[4]=="useMHT"):
        if("ST_tight" in PlotsDir.GetName()):
            ExcHistName    =("stExc%02iHistMHT_tight"%j)
            IncHistName    =("stInc%02iHistMHT_tight"%j)
            Exc02HistName  =("stExc02HistMHT_tight")
            Exc03HistName  =("stExc03HistMHT_tight")
            Exc04HistName  =("stExc04HistMHT_tight")
            pass
        else:
            ExcHistName    =("stExc%02iHistMHT"%j)
            IncHistName    =("stInc%02iHistMHT"%j)
            Exc02HistName  =("stExc02HistMHT")
            Exc03HistName  =("stExc03HistMHT")
            Exc04HistName  =("stExc04HistMHT")
            pass
        pass

    stExcHist =PlotsDir.Get(ExcHistName  )
    stIncHist =PlotsDir.Get(IncHistName  )
    stExc2Hist=PlotsDir.Get(Exc02HistName)
    stExc3Hist=PlotsDir.Get(Exc03HistName)
    stExc4Hist=PlotsDir.Get(Exc04HistName)

    if rebin:
        stExc2Hist.Rebin()
        stExc3Hist.Rebin()
        stExcHist.Rebin()
        stIncHist.Rebin()
        pass
    if j>6:
        stIncHist.Rebin()
    #if j==2:
    #    NormAndDrawST(stExcHist,j,"Exc",stExc3Hist,stExc2Hist,True)
    if j==3: NormAndDrawST(stExcHist,j,"Exc",stExc2Hist,True)
    if j==4: NormAndDrawST(stExcHist,j,"Exc",stExc3Hist,True)

    #NormAndDrawST(stIncHist,j,"Inc",stExc2Hist,True)
    NormAndDrawST(stIncHist,j,"Inc",stExc3Hist,True)
    #NormAndDrawST(stExcHist,j,"Exc",stExc2or3Hist,False)
    #NormAndDrawST(stIncHist,j,"Inc",stExc2or3Hist,False)

c1= TCanvas("chi2graph","Chi2 vs N", 800,600)
OutFile.Append(c1)

chi2norm=[]
chi2norm.append(["Inclusive multiplicity","chi2 after normalization"])

leg = TLegend(0.7,0.7,0.9,0.9)
chi2graphs_norm[chi2graphs_norm.keys()[0]].Draw()
for gname in chi2graphs_norm:
   print "now drawing chi2 graphs %s"%gname
   chi2graphs_norm[gname].Draw("SAME")
   leg.AddEntry(chi2graphs_norm[gname],gname,"L")

for i in range(chi2graphs_norm["ATLAS1_exc3"].GetN()):
    x=Double(0.0)
    y=Double(0.0)
    chi2graphs_norm["ATLAS1_exc3"].GetPoint(i,x,y)
    print x,y
    

chi2graphs_norm[chi2graphs_norm.keys()[0]].SetTitle("Chi2/Ndof for different fit functions after normalization")
chi2graphs_norm[chi2graphs_norm.keys()[0]].GetXaxis().SetTitle("Inclusive Multiplicity")
chi2graphs_norm[chi2graphs_norm.keys()[0]].GetYaxis().SetTitle("Chi2/Ndof")
chi2graphs_norm[chi2graphs_norm.keys()[0]].GetYaxis().SetRangeUser(0,1)
leg.SetFillStyle(1001);
leg.SetFillColor(0);
leg.Draw()
c1.Write()
# Write the chi2pNDF graphs for fitting
c1.Clear()
leg.Clear()
c1.SetName("chi2graph_fit")
chi2graphs_fit[chi2graphs_fit.keys()[0]].Draw()

#for gname in chi2graphs_fit:
   #print "now drawing chi2 graphs %s"%gname
chi2graphs_fit[chi2graphs_fit.keys()[0]].GetXaxis().SetTitle("Exclusive Multiplicity")
chi2graphs_fit[chi2graphs_fit.keys()[0]].GetYaxis().SetTitle("Chi2/Ndof")
chi2graphs_fit[chi2graphs_fit.keys()[0]].GetYaxis().SetRangeUser(0,2)
leg.SetFillStyle(1001)
leg.SetFillColor(0)
leg.Draw()
print tabulate(chi2Table,"firstrow")
print "Median Integral: %f"%np.median(np.array(f_integrals)) 
print "Mean Integral: %f"%np.mean(np.array(f_integrals)) 
print "(min+max)/2 Integral: %f"%( (np.amin(f_integrals)+np.amax(np.array(f_integrals)))/2.0 )
print np.sort(np.array(f_integrals)) 
print "Median chi2: %f"%np.median(np.array(f_chi2))
#for row in chi2Table:
#    print row
c1.Write()

#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]
