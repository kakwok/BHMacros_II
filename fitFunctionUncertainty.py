###################
# FCN=15.239 FROM MINOS     STATUS=SUCCESSFUL    138 CALLS        3060 TOTAL
#                     EDM=7.12799e-11    STRATEGY= 1      ERROR MATRIX ACCURATE
#  EXT PARAMETER                                   STEP         FIRST
#  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
#   1  p0           3.07101e+00   1.94152e-01  -9.80186e-05  -1.39547e+00
#   2  p1          -8.10968e-01   1.83003e-01  -1.00651e-04   1.35899e+00
#   3  p2           7.47502e-01   4.91735e-03   4.91735e-03   3.14992e+05
#
#errors for parameters are 0.183003, 0.004917, and 0.000000

from ROOT import *
# Create one array for each combination of +'s and -'s for each parameter
fs=[]

#Parameters from MINOS fitting
p0     =  1.59439e+09
p1     =  3.07101e+00
p2     = -8.10968e-01
p3     =  7.47502e-01
o1     =  3.07101e+00
o2     = -8.10968e-01
o3     =  7.47502e-01
sigma0 =  9.26517e+07
sigma1 =  1.94152e-01
sigma2 =  1.83003e-01
sigma3 =  4.91735e-03

canvas = TCanvas("canvas","canvas",800,800)
canvas.cd()
canvas.SetLogy()
f1 =TF1("f1", "[0]/([1]*x)**[2]", 1000, 7000)
f1.SetParameters(2.64090e+09, 3.51875e-03, 7.26789e+00)
f2 = TF1("f2", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
f2.SetParameters(p0, p1, p2, p3)
f2normalization=f2.Integral(1400, 2400)


for i in range(0, 3):
    for j in range(0, 2):
        fs.append(f2.Clone())
        if i==0:
            print "i is 0"
            p1=o1+((float(j)-0.5)*2)*sigma1
            p2=o2
            p3=o3
        if i==1:
            print "i is 1"
            p1=o1
            p2=o2+((float(j)-0.5)*2)*sigma2
            p3=o3
        if i==2:
            print "i is 2"
            p1=o1
            p2=o2
            p3=o3+((float(j)-0.5)*2)*sigma3
        fs[i].SetParameters(p0, p1, p2, p3)
        print "This scenario corresponds to p1 + %1.1f sigma1, p2 + %1.1f sigma2, p3 + %1.1f sigma3"%((p1-o1)/sigma1, (p2-o2)/sigma2, (p3-o3)/sigma3)

        normalization=fs[i].Integral(1400,2400)
        fs[i].SetParameter(0, p0 * f2normalization/normalization)
        print "function %i has parameters %f %f %f %f"%(i, fs[i].GetParameters()[0], fs[i].GetParameters()[1], fs[i].GetParameters()[2], fs[i].GetParameters()[2])
        fs[i].SetLineColorAlpha(kBlue, 0.5)
        if i==0:
            fs[i].Draw()
            fs[i].GetYaxis().SetRangeUser(0.5, 10e5)
            fs[i].GetXaxis().SetTitle("S_{T} (GeV)")
            fs[i].GetYaxis().SetTitleOffset(1.2)
            fs[i].GetYaxis().SetTitle("Predicted events from background fit")
        else:
            fs[i].Draw("SAME")
        #if __name__ == '__main__':
        #   rep = ''
        #   while not rep in [ 'n', 'N' ]:
        #      rep = raw_input( 'enter "n" to plot next function: ' )
        #      if 1 < len(rep):
        #         rep = rep[0]

f2.SetLineColor(kBlack)
f2.Draw("SAME")

f1.SetLineColor(kBlack)
f1.SetLineStyle(6)
f1.Draw("SAME")
#f3 = TF1("f3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
#f3.SetParameters(8214546978325997612131418112.000000, 478497.649273, -456.304620, 3.762602)
#f3.SetLineColor(kGreen)
#f3.SetLineStyle(6)
#f3.Draw("SAME")
#fLow=TF1("fLow", "(2*([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))-[4]/([5] + [6]*x + x**2)**[7])", 1000, 7000)
#fLow.SetParameters(f2.GetParameter(0),f2.GetParameter(1),f2.GetParameter(2), f2.GetParameter(3), f3.GetParameter(0), f3.GetParameter(1), f3.GetParameter(2),f3.GetParameter(3))
fLow=TF1("fLow", "(2*([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))-[4]/([5]*x)**[6])", 1000, 7000)
fLow.SetParameters(f2.GetParameter(0),f2.GetParameter(1),f2.GetParameter(2), f2.GetParameter(3), f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2))
fLow.SetLineColor(kBlack)
fLow.SetLineStyle(6)
fLow.Draw("SAME")

canvas.Print("fitfunctions.pdf")
canvas.Print("fitfunctions.png")
