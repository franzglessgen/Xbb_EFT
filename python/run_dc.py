#!/usr/bin/env python
from __future__ import print_function
import sys, ROOT, warnings
sys.path.append("/work/fglessge/EFT/CMSSW_10_1_0/src/Xbb/python/myutils")
ROOT.gROOT.SetBatch(True)
ROOT.v5.TFormula.SetMaxima(10000)
#suppres the EvalInstace conversion warning bug
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )
from optparse import OptionParser
from myutils.Datacard import Datacard

# ------------------------------------------------------------------------------
# script to produce datacards from cached tree
# recommended to run it for one sample identifier per job
# ------------------------------------------------------------------------------
class RunDatacards(object):

    def __init__(self, config, region, chunkNumber=-1, useSampleIdentifiers=None, verbose=False, forceRedo=True, runEFTcomponent = False, EFTcomponent = None, EFTcomponentname = None):
        self.verbose = verbose
        self.forceRedo = forceRedo
        self.config = config
        self.region = region
        self.chunkNumber = chunkNumber
        self.useSampleIdentifiers = useSampleIdentifiers
        self.dcMaker = Datacard(config=self.config, region=region, runEFTcomponent = runEFTcomponent, EFTcomponent = EFTcomponent, EFTcomponentname = EFTcomponentname)

    def run(self):
        if self.dcMaker:
            if not self.dcMaker.splitFilesExist(useSampleIdentifier=self.useSampleIdentifiers, chunkNumber=self.chunkNumber) or self.forceRedo:
                self.dcMaker.run(useSampleIdentifiers=self.useSampleIdentifiers, chunkNumber=self.chunkNumber)
            else:
                print ("nothing to do.")
        else:
            print ("\x1b[31mERROR: could not produce datacards\x1b[0m")

if __name__ == "__main__":
    # read arguments
    argv = sys.argv
    parser = OptionParser()
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                              help="Verbose mode.")
    parser.add_option("-C", "--config", dest="config", default=[], action="append",
                          help="configuration file")
    parser.add_option("-t", "--regions", dest="regions", default='',
                          help="cut regions identifier")
    parser.add_option("-S","--sampleIdentifier", dest="sampleIdentifier", default='',
                          help="sample identifier (no subsample!)")
    parser.add_option("-f", "--force", action="store_true", dest="force", default=False,
                          help="force overwriting of already existing datacards")
    parser.add_option("-i", "--chunkNumber", dest="chunkNumber", default='-1',
                          help="number of part to cache")
    parser.add_option("-E", "--runEFTcomponent", dest="runEFTcomponent", default=False,
                          help="Used if the sample has to be run multiple times for different EFT weights")
    parser.add_option("-W", "--EFTcomponent", dest="EFTcomponent", default=None,
                          help="Index of the EFT component to process")
    parser.add_option("-N", "--EFTcomponentname", dest="EFTcomponentname", default=None,
                          help="Datacard name of the EFT process")
    (opts, args) = parser.parse_args(argv)
    if opts.config == "":
            opts.config = "config"

    # Import after configure to get help message
    from myutils import BetterConfigParser

    # load config
    config = BetterConfigParser()
    config.read(opts.config)

    # if no region is given in argument, run it for all of them
    regionsListString = opts.regions if len(opts.regions.strip())>0 else config.get('LimitGeneral', 'List')
    regions = [x.strip() for x in regionsListString.split(',') if len(x.strip()) > 0]
    useSampleIdentifiers = opts.sampleIdentifier.split(',') if len(opts.sampleIdentifier) > 0 else None
    print("regions:", regions)
    for region in regions:
        print("init...", opts.chunkNumber)
        runDC = RunDatacards(config=config, region=region, chunkNumber=int(opts.chunkNumber), useSampleIdentifiers=useSampleIdentifiers, forceRedo=opts.force, runEFTcomponent = opts.runEFTcomponent, EFTcomponent = opts.EFTcomponent, EFTcomponentname = opts.EFTcomponentname)
        print("run...", opts.chunkNumber)
        runDC.run()
