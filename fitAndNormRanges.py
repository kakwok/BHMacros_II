# Class for dealing with fit and normalization ranges
# Initialize the class by calling the constructor with the name of the textfile with the fit/normrange structure
# The fit and norm ranges specified in the textfile must have the structure:
# FitRange, exc, 2,  1300, 2800
# NormRange, exc, 11,  2100, 2300
# There are simple getters to get the ranges and their bounds, for example
# FitAndNormRange.getLowerNormBound("inc8")
# John Hakala, 12/2/2015
from collections import defaultdict

class FitAndNormRange:
    # Loads the fit and normalization ranges for a given set of ST spectra from a text file
    # Makes them available for the fitSThists pyroot macro
    def __init__(self, rangesFile):
        self.rangesFile = rangesFile
        fitRangesArray = []
        normRangesArray = []
        plotRangesArray = []
        self.FitRanges = defaultdict()
        self.NormRanges = defaultdict()
        self.PlotRanges = defaultdict()
        with open(self.rangesFile) as rangesTextFile:
            for line in rangesTextFile.readlines():
                rangeArray = line.split(",")
                key='{}{}'.format(rangeArray[1].strip(), rangeArray[2].strip())
                value=[float(rangeArray[3].strip()), float(rangeArray[4].strip())]
                if rangeArray[0].strip() == "FitRange" :
                    fitRangesArray.append([key, value])
                elif rangeArray[0].strip() == "NormRange" :
                    normRangesArray.append([key, value])
                elif rangeArray[0].strip() == "PlotRange" :
                    plotRangesArray.append([key, value])

        for dist, fitrange in fitRangesArray:
            self.FitRanges[dist]=fitrange
        for dist, normrange in normRangesArray:
            self.NormRanges[dist]=normrange
        for dist, plotrange in plotRangesArray:
            self.PlotRanges[dist]=plotrange

    def showFitRanges(self):
        print "The fit ranges are:"
        for key in self.FitRanges:
            print key,
            print ": ",
            print self.FitRanges[key]

    def showNormRanges(self):
        print "The normalization ranges are:"
        for key in self.NormRanges:
            print key,
            print ": ",
            print self.NormRanges[key]

    def getFitRange(self, key):
        return self.FitRanges[key]

    def getLowerFitBound(self, key):
        return self.FitRanges[key][0]

    def getUpperFitBound(self, key):
        return self.FitRanges[key][1]

    def getNormRange(self, key):
        return self.NormRanges[key]

    def getLowerNormBound(self, key):
        return self.NormRanges[key][0]

    def getUpperNormBound(self, key):
        return self.NormRanges[key][1]
    def getLowerPlotRange(self, key):
        return self.PlotRanges[key][0]
    def getUpperPlotRange(self, key):
        return self.PlotRanges[key][1]
