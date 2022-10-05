import os
import json

class MuonSF:
    def __init__(self, mu_json, mu_name) :
        self.init(mu_json, mu_name)

    def init(self, mu_json, mu_name) :
        f = open(mu_json, 'r')
        results = json.load(f)
        if mu_name not in results.keys():
            return False
        self.res = results[mu_name]
        f.close()

    def get_2D(self, pt, eta):
        eta = abs(eta)
        header = ''
        if "abseta_pt_ratio" in self.res.keys(): header = "abseta_pt_ratio"
        elif "abseta_pt_MC" in self.res.keys(): header = "abseta_pt_MC"
        if header != '':
            for etaKey, values in sorted(self.res[header].iteritems()) :
                etaL = float(((etaKey[8:]).rstrip(']').split(',')[0]))
                etaH = float(((etaKey[8:]).rstrip(']').split(',')[1]))
                if not (eta>etaL and eta<etaH):
                    continue
                #print etaL, etaH
                for ptKey, result in sorted(values.iteritems()) :
                    ptL = float(((ptKey[4:]).rstrip(']').split(',')[0]))
                    ptH = float(((ptKey[4:]).rstrip(']').split(',')[1]))
                    if not (pt>ptL and pt<ptH):
                        continue
                    #print ptL, ptH
                    #print "|eta| bin: %s  pT bin: %s\tdata/MC SF: %f +/- %f" % (etaKey, ptKey, result["value"], result["error"])
                    return [result["value"], result["error"]]
        else:
            return [1.0, 0.0]

        # if nothing was found, return 1 +/- 0
        return [1.0, 0.0]



##################################################################################################
# EXAMPLE
#
#jsons = {
#    'json/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.json' : 'runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins',
#    'json/MuonIso_Z_RunCD_Reco74X_Dec1.json' : 'NUM_LooseRelIso_DEN_LooseID_PAR_pt_spliteta_bin1',
#    'json/MuonID_Z_RunCD_Reco74X_Dec1.json' : 'NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1' ,
#    'json/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1_MC.json' : 'runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins'
#    }
#
## example
#pt = 40.01
#eta = -1.68
#
#for j, name in jsons.iteritems():
#    muonCorr = MuonSF(j , name)
#    weight = muonCorr.get_2D( pt , eta)
#    val = weight[0]
#    err = weight[1]
#    print j, name, ': ',  val, ' +/- ', err
#
