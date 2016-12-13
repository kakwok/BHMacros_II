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
STcomparisons = {}   #  Dictionary of canvas
upperPads     = {}   #  Dictionary of upperPads in Canvas
lowerPads     = {}   #  Dictionary of lowerPads in Canvas
Chi2List     = []
f2_list      = {}     #  functions for fitted to ST n=2  
f3_list      = {}     #  functions for fitted to ST n=3
f2_norm_list = {}     #  functions for fitted to ST n=2 after normalization
f3_norm_list = {}     #  functions for fitted to ST n=2 after normalization
PlotsFile = TFile(argv[1])
#PlotsDir = PlotsFile.Get("ST")
PlotsDir = PlotsFile.Get("ST_tight")
OutFile = TFile("output/%s"%argv[2], "RECREATE")
chi2graphs   = {}
fitNormRanges = FitAndNormRange(argv[3])
fitNormRanges.showFitRanges()
fitNormRanges.showNormRanges()
rebin          = False 
WriteDataCards = False
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
        
def getNormalizedFunction(f, hist, xlowbin, xupbin, xlowedge, xupedge, binwidth):
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
 
###################################
##  This function takes stHist, get normalized fit function, symmetrize it and draw it on a canvas
##
##  stHist    = Inclusive/Exclusive ST histogram from flatuple  
##  j         = multiplicity of the histogram
## ExcOrInc   = "Exc" or "Inc"
## stExc2Hist = histogram of N=2 for drawing ratio plot 
###################################
def FitAndDrawST(stHist,j,ExcOrInc,stExc2Hist):
    if rebin :
        stHist.Rebin(10)
    # For exclusive STs
    canvasName = "st%s%02iCanvas"%(ExcOrInc,j)
    UpperPadName    = "%s%02ipad"%(ExcOrInc,j)
    STcomparisons[canvasName] = TCanvas(canvasName, "ST, N=%i"%j, 800, 800)
    STcomparisons[canvasName].cd()
    upperPads[UpperPadName]        = (TPad("%s%02ipad"%(ExcOrInc,j), "pad1", 0, 0.3, 1, 1.0))
    upperPads[UpperPadName].SetBottomMargin(2)
    upperPads[UpperPadName].Draw()
    upperPads[UpperPadName].cd()
    upperPads[UpperPadName].SetLogy()
    stHist.GetXaxis().SetTitle("S_{T} (GeV)")
    stHist.GetYaxis().SetTitle("Events")
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
    stHist.GetXaxis().SetRangeUser(lowerNormEdge, STup)
    #stHist.GetXaxis().SetRangeUser(1400, STup)
    stHist.SetMinimum(1e-1)

    print "in N=%i, upperNormEdge = %s, lowerNormEdge = %s" % (j,upperNormEdge, lowerNormEdge)
    chi2_list     = []
    #chi2 deviation list
    chi2_devlist  = []
    functions     = []
    # Normalize the N=2 fit functions
    for fname in f2_list:
    	fname_norm = fname +"_norm"
    	f2_norm_list[fname_norm] = getNormalizedFunction( f2_list[fname], stHist, lowerNormBin, upperNormBin, lowerNormEdge, upperNormEdge, binwidth)
        chi2 = stHist.Chisquare( f2_norm_list[fname_norm] ,"R")/ f2_norm_list[fname_norm].GetNDF()
    	chi2_list.append(chi2)
    	chi2_devlist.append(abs(chi2-1))
    	functions.append( f2_norm_list[fname_norm] )
    	if(ExcOrInc=="Inc"):
            chi2graphs[fname].SetPoint( chi2graphs[fname].GetN(), j , chi2)

    # Draw the normalized N=2 fit functions
    for fname in f2_norm_list:
    	f2_norm_list[fname].Draw("SAME")

    # Normalize the N=3 fit functions
    for fname in f3_list:
    	fname_norm = fname +"_norm"
    	f3_norm_list[fname_norm] = getNormalizedFunction( f3_list[fname], stHist, lowerNormBin, upperNormBin, lowerNormEdge, upperNormEdge, binwidth)
        chi2 = stHist.Chisquare( f3_norm_list[fname_norm] ,"R")/ f3_norm_list[fname_norm].GetNDF()
    	chi2_list.append(chi2)
    	chi2_devlist.append(abs(chi2-1))
    	functions.append( f3_norm_list[fname_norm] )
    	# Fill chi2 graph
    	if(ExcOrInc=="Inc"):
            chi2graphs[fname].SetPoint( chi2graphs[fname].GetN(), j , chi2)

    # Draw the normalized N=3 fit functions
    for fname in f3_norm_list:
    	f3_norm_list[fname].Draw("SAME")
	
    
    #fbest    = f2Normalized
    #fLow     = getSymmetrizedFuntion( fbest, functions, upperNormEdge, 14000)
    #fbest    = f2_norm_list["f2_norm"]
    fbest     = functions[ chi2_devlist.index( min(chi2_devlist) ) ]
    print "fbest is chosen to be %s\n"%fbest.GetName()
    fLow     = getSymmetrizedFuntion( fbest, functions, 3000, 14000)
    
    fLow.SetLineColor(kBlue)
    fLow.SetLineStyle(6)
    fLow.Draw("SAME")
    chi2_list    = []
    chi2_devlist = []
    functions    = []

    if(ExcOrInc=="Exc"):
        legend = TLegend(0.35, 0.68, 0.8, 0.85, "Exclusive ST distributions for n=%i"%j, "brNDC")
    if(ExcOrInc=="Inc"):
        legend = TLegend(0.35, 0.68, 0.8, 0.85, "Inclusive ST distributions for n>=%i"%j, "brNDC")
   
    legend.SetTextSize(0.04);
    legend.SetLineWidth(1);
    legend.SetBorderSize(0);
    legend.SetFillStyle(1001);
    legend.SetFillColor(10);
    legend.AddEntry(stHist,"data","lp");
    legend.Draw()

    LowerPadName = "%s%02iratiopad"%(ExcOrInc,j)
    STcomparisons[canvasName].cd()
    lowerPads[LowerPadName] = (TPad(LowerPadName, "ratiopad1", 0, 0.04, 1, 0.3))
    lowerPads[LowerPadName].SetTopMargin(5)
    lowerPads[LowerPadName].SetBottomMargin(0.2)
    lowerPads[LowerPadName].Draw()
    lowerPads[LowerPadName].cd()

    stExcRatio = stHist.Clone("stExcRatio")
    stExcRatio.GetXaxis().SetLabelSize(.08)
    stExcRatio.GetXaxis().SetTitle("")
    stExcRatio.GetYaxis().SetTitle("Ratio to n=2")
    stExcRatio.GetYaxis().SetLabelSize(.075)
    stExcRatio.GetYaxis().SetTitleSize(.1)
    stExcRatio.GetYaxis().SetTitleOffset(.3)
    stExcRatio.Divide(stExc2Hist)
    stExcRatio.SetTitle("")
    stExcRatio.SetMarkerColor(kBlack)
    stExcRatio.SetMarkerStyle(8)
    stExcRatio.SetMarkerSize(0.7)
    stExcRatio.Draw("EP")
    stExcRatio.SetStats(0)
    STcomparisons[canvasName].Write()

    if (WriteDataCards and ExcOrInc=="Inc"):
        outputForLimits = open("output/%s_Inclusive%i.txt"%(argv[2],j), "w")
        outputForLimits.write(" STMin    ::   Observed Data   ::   Expected Bkg   ::  Shape Unc  \n")
        for stmin in range(20, 75):
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

f2_list["f1"]=f1
f2_list["f2"]=f2
f2_list["f3"]=f3
f2_list["f4"]=f4
#f2_list["f5"]=f5

f3_list["f1_exc3"]=f1_exc3
f3_list["f2_exc3"]=f2_exc3
f3_list["f3_exc3"]=f3_exc3
f3_list["f4_exc3"]=f4_exc3
#f3_list["f5_exc3"]=f5_exc3

chi2graphs["f1"]=TGraph()
chi2graphs["f2"]=TGraph()
chi2graphs["f3"]=TGraph()
chi2graphs["f4"]=TGraph()
chi2graphs["f1_exc3"]=TGraph()
chi2graphs["f2_exc3"]=TGraph()
chi2graphs["f3_exc3"]=TGraph()
chi2graphs["f4_exc3"]=TGraph()


fLow=TF1("fLow", "[0]*(2*[1]/([2]*x)**[3]-[4]/([5] + [6]*x + x**2)**[7])", 1000, STup)
Funcs ={0:f1,1:f2,2:f3,3:f1_exc3,4:f2_exc3,5:f3_exc3,6:fLow}

##############################################
##  Fit the ST histogram, get the functions
##############################################
for i in range(2,4):
    if(argv[4]=="useMET"):
        if( "ST_tight" in PlotsDir.GetName()):
            histname = ("stExc%02iHist_tight"%i)
        else:
            histname = ("stExc%02iHist"%i)
    if(argv[4]=="useMHT"):
        if( "ST_tight" in PlotsDir.GetName()):
            histname = ("stExc%02iHistMHT_tight"%i)
        else:
            histname = ("stExc%02iHistMHT"%i)

    stExcHist =PlotsDir.Get(histname)

    if (i == 2 or i==3 ):
            if i==2:
		# John's parameter
                #f1.SetParameters(2.08e9, 4.28e-3, 6.7)
                #f2.SetParameters(1.7e12, 2.1, .16, 0.62)
                #f3.SetParameters(2e27, 4e5, -4e2, 3.6)
                #f4.SetParameters(8e11, -0.86,  0.62)
                #f5.SetParameters(1.7e12, -1, -3, 0.62)
		# Tutanon's function parameters
                f2_list["f1"].SetParameters(8e6, 0.5, 9)
                f2_list["f2"].SetParameters(2.4e6, -3, 4.7, 0.4)
                f2_list["f3"].SetParameters(6e5, 0.4, -0.1, 4)
                f2_list["f4"].SetParameters(1.5e9, -12,  -0.7)
                #f2_list["f5"].SetParameters(2e5, 0, 6, 0.7)

		for fname in f2_list:
			print stExcHist.GetName()
			customfit(f2_list[fname],stExcHist,"exc2")
			f2_list[fname].SetLineStyle(6)

                f2_list["f1"].SetLineColor(kRed)
                f2_list["f2"].SetLineColor(kBlack)
                f2_list["f3"].SetLineColor(kViolet)
                f2_list["f4"].SetLineColor(kPink)
                #f2_list["f5"].SetLineColor(kYellow)

            if i==3:
		# John's parameter
                #f1_exc3.SetParameters(2.08e9, 4.28e-3, 6.7)
                #f2_exc3.SetParameters(1.7e12, 2.1, .16, 0.62)
                #f3_exc3.SetParameters(2e27, 4e5, -4e2, 3.6)
                #f4_exc3.SetParameters(8e11, -0.86,  0.62)
                #f5_exc3.SetParameters(1.7e12, -1, -3, 0.62)
		# Tutanon's function parameter
                f3_list["f1_exc3"].SetParameters(8e6, 0.5, 9)
                f3_list["f2_exc3"].SetParameters(2.4e6, -3, 4.7, 0.4)
                f3_list["f3_exc3"].SetParameters(6e5, 0.4, -0.1, 4)
                f3_list["f4_exc3"].SetParameters(1.5e9, -12,  -0.7)
                #f3_list["f5_exc3"].SetParameters(2e5, 0, 6, 0.7)

		for fname in f3_list:
			customfit(f3_list[fname],stExcHist,"exc3")
			f3_list[fname].SetLineStyle(4)
                f3_list["f1_exc3"].SetLineColor(kGreen)
                f3_list["f2_exc3"].SetLineColor(kCyan)
                f3_list["f3_exc3"].SetLineColor(kMagenta)
                f3_list["f4_exc3"].SetLineColor(kOrange)
                #f3_list["f5_exc3"].SetLineColor(kSpring+10)

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

for j in range(2,11):
    if (argv[4]=="useMET"):
        if("ST_tight" in PlotsDir.GetName()):
           ExcHistName    =("stExc%02iHist_tight"%j)
           IncHistName    =("stInc%02iHist_tight"%j)
           Exc02HistName  =("stExc02Hist_tight")
        else:
           ExcHistName    =("stExc%02iHist"%j)
           IncHistName    =("stInc%02iHist"%j)
           Exc02HistName  =("stExc02Hist")
          
    if (argv[4]=="useMHT"):
        if("ST_tight" in PlotsDir.GetName()):
           ExcHistName    =("stExc%02iHistMHT_tight"%j)
           IncHistName    =("stInc%02iHistMHT_tight"%j)
           Exc02HistName  =("stExc02HistMHT_tight")
        else:
           ExcHistName    =("stExc%02iHistMHT"%j)
           IncHistName    =("stInc%02iHistMHT"%j)
           Exc02HistName  =("stExc02HistMHT")

    stExcHist =PlotsDir.Get(ExcHistName  )
    stIncHist =PlotsDir.Get(IncHistName  )
    stExc2Hist=PlotsDir.Get(Exc02HistName)
    
    FitAndDrawST(stExcHist,j,"Exc",stExc2Hist)
    FitAndDrawST(stIncHist,j,"Inc",stExc2Hist)

c1= TCanvas("chi2graph","Chi2 vs N", 800,600)
OutFile.Append(c1)

chi2graphs["f1"].SetLineColor(kRed)
chi2graphs["f2"].SetLineColor(kBlack)
chi2graphs["f3"].SetLineColor(kViolet)
chi2graphs["f4"].SetLineColor(kPink)
chi2graphs["f1_exc3"].SetLineColor(kGreen)
chi2graphs["f2_exc3"].SetLineColor(kCyan)
chi2graphs["f3_exc3"].SetLineColor(kMagenta)
chi2graphs["f4_exc3"].SetLineColor(kOrange)
#chi2graph_f5.SetLineColor(kYellow)
#chi2graph_f5_exc3.SetLineColor(kSpring+10)

leg = TLegend(0.7,0.7,0.9,0.9)
chi2graphs["f1"].Draw()
for gname in chi2graphs:
   print "now drawing chi2 graphs %s"%gname
   chi2graphs[gname].Draw("SAME")
   leg.AddEntry(chi2graphs[gname],gname,"L")
chi2graphs["f1"].SetTitle("Chi2/Ndof for different fit functions")
chi2graphs["f1"].GetXaxis().SetTitle("Inclusive Multiplicity")
chi2graphs["f1"].GetYaxis().SetTitle("Chi2/Ndof")
chi2graphs["f1"].GetYaxis().SetRangeUser(0,50)
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
