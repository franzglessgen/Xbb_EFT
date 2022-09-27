import os
import json

# A class to apply SF's tabulated in json files
class LeptonSF:
    def __init__(self, lep_json, lep_name, lep_binning, extrapolateFromClosestBin=True) :
        if not os.path.isfile(lep_json):
            self.valid = False
            if lep_json!="":
                pass
                print "[LeptonSF]: Warning: ", lep_json, " is not a valid file. Return."
            else:
                pass
                print "[LeptonSF]: No file has been specified. Return."
        else:
            self.init(lep_json, lep_name, lep_binning, extrapolateFromClosestBin)

    def init(self, lep_json, lep_name, lep_binning, extrapolateFromClosestBin) :
        f = open(lep_json, 'r')             
        #print '[LeptonSF]: Initialize with the following parameters:'
        #print '\tfile:',lep_json
        #print '\titem:', lep_name
        #print '\tbinning:', lep_binning
        results = json.load(f)
        if lep_name not in results.keys():
            self.valid = False
            #print "[LeptonSF]: Warning: ", lep_name , " is not a valid item. Return."
            return False
        self.res = results[lep_name]
        self.lep_name = lep_name
        self.lep_binning = lep_binning
        self.valid = True
        self.extrapolateFromClosestBin = extrapolateFromClosestBin
        f.close()

    def get_1D(self, pt):
        if not self.valid:
            return [1.0, 0.0]

        if self.lep_binning not in self.res.keys():
            return [1.0, 0.0]

        # if no bin is found, search for closest one, and double the uncertainty
        closestPtBin = ""
        closestPt = 9999.

        ptFound = False

        for ptKey, result in sorted(self.res[self.lep_binning].iteritems()) :
            ptL = float(((ptKey[ptKey.find(':')+2:]).rstrip(']').split(',')[0]))
            ptH = float(((ptKey[ptKey.find(':')+2:]).rstrip(']').split(',')[1]))

            #print 'ptL is', ptL
            #print 'ptH is', ptH

            if abs(ptL-pt)<closestPt or abs(ptH-pt)<closestPt and not ptFound:
                closestPt = min(abs(ptL-pt), abs(ptH-pt))
                closestPtBin = ptKey

                if (pt>ptL and pt<ptH):

                    closestPtBin = ptKey
                    ptFound = True

                if ptFound:
                    return [result["value"], result["error"]]

        if self.extrapolateFromClosestBin and not (closestPtBin==""):
            return [self.res[self.lep_binning][closestPtBin]["value"],2*self.res[self.lep_binning][closestPtBin]["error"]]
        else:
            return [1.0, 0.0]
                    

    def get_2D(self, pt, eta):
        if not self.valid:
            return [1.0, 0.0]        

        stripForPt= 4

        stripForEta = 5
        if self.lep_binning not in self.res.keys():
            return [1.0, 0.0]

        if "abseta" in self.lep_binning:
            eta = abs(eta)
            stripForEta = 8

        if "pt_eta_ratio" in self.lep_binning:
            stripForEta = 4
            stripForPt= 5

        #if "abseta_pt_DATA" in self.lep_binning or "abseta_pt_MC" in self.lep_binning:
        #    stripForEta = 8
        #    stripForPt= 4

        # if no bin is found, search for closest one, and double the uncertainty
        closestEtaBin = ""
        closestPtBin = ""
        closestEta = 9999.
        closestPt = 9999.

        etaFound = False
        for etaKey, values in sorted(self.res[self.lep_binning].iteritems()) :
            #print 'etaKey is', etaKey
            #print 'etaKey after strip is', etaKey[stripForEta:]
            #etaL = float(((etaKey[stripForEta:]).rstrip(']').split(',')[0]))
            #etaH = float(((etaKey[stripForEta:]).rstrip(']').split(',')[1]))
            etaL = float(((etaKey[etaKey.find('[')+1:]).rstrip(']').split(',')[0]))
            etaH = float(((etaKey[etaKey.find('[')+1:]).rstrip(']').split(',')[1]))

            ptFound = False

            if abs(etaL-eta)<closestEta or abs(etaH-eta)<closestEta and not etaFound:
                closestEta = min(abs(etaL-eta), abs(etaH-eta))
                closestEtaBin = etaKey

            if (eta>etaL and eta<etaH):
                closestEtaBin = etaKey
                etaFound = True                

            for ptKey, result in sorted(values.iteritems()) :
                #ptL = float(((ptKey[stripForPt:]).rstrip(']').split(',')[0]))
                #ptH = float(((ptKey[stripForPt:]).rstrip(']').split(',')[1]))
                ptL = float(((ptKey[ptKey.find('[')+1:]).rstrip(']').split(',')[0]))
                ptH = float(((ptKey[ptKey.find('[')+1:]).rstrip(']').split(',')[1]))

                if abs(ptL-pt)<closestPt or abs(ptH-pt)<closestPt and not ptFound:
                    closestPt = min(abs(ptL-pt), abs(ptH-pt))
                    closestPtBin = ptKey

                if (pt>ptL and pt<ptH):
                    closestPtBin = ptKey
                    ptFound = True

                if etaFound and ptFound:
                    return [result["value"], result["error"]]

        if self.extrapolateFromClosestBin and not (closestPtBin=="" or closestEtaBin==""):
            return [self.res[self.lep_binning][closestEtaBin][closestPtBin]["value"], 
                    2*self.res[self.lep_binning][closestEtaBin][closestPtBin]["error"]] 
        else:
            return [1.0, 0.0]


##################################################################################################
# EXAMPLE 
#

if __name__ == "__main__":


    #analysis = '2016'
    analysis = '2017'


    wdir = os.environ['CMSSW_BASE']+"/src/Xbb"
    if analysis == 2016:
        jsons = {
            #
            #Muon
            #
            #ID and ISO
            wdir+'/python/json/V25/muon_ID_BCDEFv2.json' : ['MC_NUM_TightID_DEN_genTracks_PAR_pt_eta', 'abseta_pt_ratio'], #eta pt
            wdir+'/python/json/V25/muon_ID_GHv2.json' : ['MC_NUM_TightID_DEN_genTracks_PAR_pt_eta', 'abseta_pt_ratio'],
            ###
            wdir+'/python/json/V25/muon_ISO_BCDEFv2.json' : ['TightISO_TightID_pt_eta', 'abseta_pt_ratio'],
            wdir+'/python/json/V25/muon_ISO_GHv2.json' : ['TightISO_TightID_pt_eta', 'abseta_pt_ratio'],
            #Tracker
            wdir+'/python/json/V25/trk_SF_RunBCDEF.json' : ['Graph', 'ratio_eff_eta3_dr030e030_corr'],
            wdir+'/python/json/V25/trk_SF_RunGH.json' : ['Graph', 'ratio_eff_eta3_dr030e030_corr'],
            #Trigg
            #BCDEF
            wdir+'/python/json/V25/EfficienciesAndSF_RunBtoF.json' : ['IsoMu24_OR_IsoTkMu24_PtEtaBins', 'abseta_pt_ratio'],
            #GH
            wdir+'/python/json/V25/theJSONfile_Period4.json' : ['IsoMu24_OR_IsoTkMu24_PtEtaBins', 'abseta_pt_ratio'],
            ##
            ##Electron
            ##
            ##ID and ISO (grouped as MVAid for electron)
            wdir+'/python/json/V25/EIDISO_WH_out.json' : ['EIDISO_WH', 'eta_pt_ratio'],
            #Tracker
            wdir+'/python/json/V25/ScaleFactor_etracker_80x.json' : ['ScaleFactor_tracker_80x', 'eta_pt_ratio'],
            #Trigg
            wdir+'/python/json/V25/Tight27AfterIDISO_out.json' : ['Tight27AfterIDISO', 'eta_pt_ratio']
            }
    elif analysis == '2017':
        jsons = {
            #
            # Muon
            #
            # ID
            wdir+'/python/json/V25/SingleLepton_2017/RunBCDEF_SF_ID.json' : ['NUM_TightID_DEN_genTracks', 'abseta_pt'],
            #
            # ISO
            wdir+'/python/json/V25/SingleLepton_2017/RunBCDEF_SF_ISO.json' : ['NUM_UltraTightIso4_DEN_TightIDandIPCut', 'abseta_pt'],
            # Tracker
            # no tracker for 2017 (~1)
            #Trigg
            wdir+'/python/json/V25/SingleLepton_2017/RunBCDEF_SF_TRIGG.json' : ['NUM_IsoMu27_DEN_empty', 'abseta_pt'],

            ##
            ##Electron
            ##
            ##ID and ISO (grouped as MVAid for electron)
            wdir+'/python/json/V25/SingleLepton_2017/VHbb1ElectronIDISO2017.json' : ['singleEleIDISO2017', 'eta_pt_ratio'],
            #Tracker
            #-------check----------
            # no tracker for 2017, at least not used in AT
            #Trigg
            wdir+'/python/json/V25/SingleLepton_2017/VHbb1ElectronTrigger2017.json' : ['singleEleTrigger', 'eta_pt_ratio']
            }

    for j, name in jsons.iteritems():

        print 'j is', j
        lepCorr = LeptonSF(j , name[0], name[1])

        #test1
        lepton_pt = 35 
        lepton_eta = 1.7

        ##test2
        #lepton_pt = 75 
        #lepton_eta = 0.9 

        #2-D binned SF
        if not j.find('trk_SF_Run') != -1:
            if 'abseta' in  name[1]:
                weight = lepCorr.get_2D(lepton_pt, abs(lepton_eta))
            else:
                weight = lepCorr.get_2D(lepton_pt, lepton_eta)
        #1-D binned SF
        else:
            weight = lepCorr.get_1D(lepton_eta)
        val = weight[0]
        err = weight[1]
        print 'SF: ',  val, ' +/- ', err
