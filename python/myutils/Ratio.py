import ROOT
def renewHist(hist,reference,min,max):
    theHist = hist.Clone()
    theReference = reference.Clone()
    theHist.SetLineWidth(1)
    theHist.SetMarkerSize(1)
    return theHist, theReference

def getRatio(hist,reference,min,max,yTitle="Data/MC",maxUncertainty = 1000.000,restrict=True, yRange=[0.5,1.75]):
    from ROOT import gROOT,gSystem
    theHist, theReference = renewHist(hist,reference,min,max)
    try:
        from ROOT import coolRatio
    except ImportError:
        gROOT.ProcessLine('.L $CMSSW_BASE/src/Xbb/python/myutils/Ratio.C')
        from ROOT import coolRatio
    thePlotter = coolRatio()
    theRatio = thePlotter.make_rebinned_ratios(theHist,theReference,maxUncertainty,False,0)
    refError = thePlotter.make_rebinned_ratios(theHist,theReference,maxUncertainty,False,1)
    theRatio.GetXaxis().SetRangeUser(min,max)
    print "theRation: range",min," ",max
    #theRatio.GetXaxis().SetRangeUser(0,1)
    if restrict:
        theRatio.SetMinimum(yRange[0])
        theRatio.SetMaximum(yRange[1])
    else:
        theRatio.SetMinimum(int(theRatio.GetMinimum()))
        theRatio.SetMaximum(int(theRatio.GetMaximum()*1.5))
    #theRatio.GetYaxis().SetNdivisions(104)
    theRatio.GetYaxis().SetNdivisions(505)
    theRatio.GetYaxis().SetTitle("Ratio")
    theRatio.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*2.2)
    theRatio.GetYaxis().SetTitleOffset(0.6)
    theRatio.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize() * 2.2)
    theRatio.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*2.2)
    theRatio.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize() * 2.2)
    theRatio.GetYaxis().SetTitleOffset(0.4)
    theRatio.GetYaxis().CenterTitle(ROOT.kTRUE)
    theRatio.GetYaxis().SetDrawOption("M")
    theRatio.SetXTitle("")
    theRatio.SetYTitle(yTitle)
    return theRatio, refError

