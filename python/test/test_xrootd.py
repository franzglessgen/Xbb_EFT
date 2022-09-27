#!/usr/bin/env python
from __future__ import print_function
import unittest
import sys
import random
sys.path.append('../')
from myutils.sampleTree import SampleTree as SampleTree
from myutils.FileLocator import FileLocator
import os

class TestXrootdMethods(unittest.TestCase):

    scratchDirectory = os.environ['SCRATCH_DIR'] if 'SCRATCH_DIR' in os.environ else '.'

    def getTree(self, path):
        fileNames = [path]
        return SampleTree(fileNames, treeName='Events')

    def test_xrootd(self):
        if 'X509_USER_PROXY' in os.environ and len(os.environ['X509_USER_PROXY'].strip()) > 0:
            path1 = 'root://xrootd-cms.infn.it//store/group/phys_higgs/hbb/ntuples/VHbbPostNano/2017/V11/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17NanoAODv4-PU2017_1282/190510_115113/0000/tree_1.root'
            tree1 = self.getTree(path1)
            print ("ENTRIES:", tree1.GetEntries())
            self.assertEqual(tree1.GetEntries(), 552904L)

            fileLocator = FileLocator()

            path2 = fileLocator.removeRedirector(path1)
            print ("PATH2:", path2)
            self.assertTrue(path2.startswith('/store/group/phys_higgs/'))
            self.assertTrue(path2.endswith('/tree_1.root'))

            path3 = fileLocator.addRedirector(redirector='root://xrootd-cms.infn.it', fileName=path2)
            self.assertEqual(path1, path3)
        else:
            print("INFO: this test is skipped because no X509 proxy certificate is found which is needed to access the files!")


if __name__ == '__main__':
    unittest.main()
