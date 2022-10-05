import os
import json
import math

# A class to apply SF's tabulated in json files
class LeptonSF:
    def __init__(self, lep_json, lep_name, lep_binning, extrapolateFromClosestBin=True) :
        if not os.path.isfile(lep_json):
            self.valid = False
            if lep_json!="":
                pass
                #print "[LeptonSF]: Warning: ", lep_json, " is not a valid file. Return."
            else:
                pass
                #print "[LeptonSF]: No file has been specified. Return."
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

def computeEff_leg(leg):
    #leg is the leg index, can be 0 or 1
    eff_leg = [1.,0.,0.]
    eff_leg[0] = weight[leg][0]
    eff_leg[1] = weight[leg][0]-weight[leg][1]
    eff_leg[2] = weight[leg][0]+weight[leg][1]
    return eff_leg

def computeEventEff_fromleg(effleg1, effleg2):
    '''returns event efficiency and relative uncertainty'''
    eff_event = [1.,0.]
    eff_event[0] = ((effleg1[0][0]**2*effleg2[1][0] + effleg1[1][0]**2*effleg2[0][0])/(effleg1[0][0] + effleg1[1][0]))
    #relative uncertainty down
    uncert_down = (abs(((effleg1[0][1]**2*effleg2[1][1] + effleg1[1][1]**2*effleg2[0][1])/(effleg1[0][1] + effleg1[1][1])) - eff_event[0])/eff_event[0])
    #relative uncertainty up
    uncert_up = (abs(((effleg1[0][2]**2*effleg2[1][2] + effleg1[1][2]**2*effleg2[0][2])/(effleg1[0][2] + effleg1[1][2])) - eff_event[0])/eff_event[0])
    eff_event[1]  = (uncert_down+uncert_up)/2.
    return eff_event



##################################################################################################
# EXAMPLE 
#

#def getDoubleMuTrigger(pt1, eta1, pt2, eta2):
#    '''Return array. [0]: nominal value of the dimuon SF, [1]: down variation (stat), [2]: up variation (stat)'''

if __name__ == "__main__":

    # in this example, we have two muons mu1 and mu2

    #mu1:
    pt1 = 65
    eta1 = 1.5
    #eta1 = 2.3

    #mu2:
    pt2 = 57 
    eta2 = 0.3

    #path to json files
    wdir = '.'
    jsons = {
        #
        #Trigger
        #
        #Leg 17
        wdir+'/RunBCDEF_data_leg17_tag8.json':['NUM_hlt_Mu17Mu8_leg17_DEN_LooseIDnISO', 'abseta_pt'],
        wdir+'/RunBCDEF_mc_leg17_tag8.json':['NUM_hlt_Mu17Mu8_leg17_DEN_LooseIDnISO', 'abseta_pt'],
        #Leg 8
        wdir+'/RunBCDEF_data_leg8.json':['NUM_hlt_Mu17Mu8_leg8_DEN_LooseIDnISO', 'abseta_pt'],
        wdir+'/RunBCDEF_mc_leg8.json':['NUM_hlt_Mu17Mu8_leg8_DEN_LooseIDnISO', 'abseta_pt'],
        #
        # for systematics
        #

        #Leg 17
        wdir+'/RunBCDEF_data_leg17.json':['NUM_hlt_Mu17Mu8_leg17_DEN_LooseIDnISO', 'abseta_pt'],
        wdir+'/RunBCDEF_mc_leg17.json':['NUM_hlt_Mu17Mu8_leg17_DEN_LooseIDnISO', 'abseta_pt'],
        }

    # contrain data and mc efficieny
    # [0] index: mu1. [1] index: mu2
    effData_leg8 = []
    effData_leg17= []
    effMC_leg8= []
    effMC_leg17 = []

    #
    # for systetmatics

    effData_leg17_sys= []
    effMC_leg17_sys = []


    # SF. [nominal, down stat systematic, up stat systematic]
    SF = [1,1,1] 

    # SYS. [nominal, down stat systematic, up stat systematic]
    SYS = [1,1,1] 
 

    # read all per-eg efficiency
    for j, name in jsons.iteritems():
        #print 'j is', j
        weight = []
        lepCorr = LeptonSF(j , name[0], name[1])

        weight.append(lepCorr.get_2D(pt=pt1, eta=eta1))
        weight.append(lepCorr.get_2D(pt=pt2, eta=eta2))

        ###
        # reading  efficiency for each leg
        ###
         
        if j.find('RunBCDEF_data_leg17_tag8.json') != -1:
            #compute the efficiency for both legs
            effData_leg17.append(computeEff_leg(0))
            effData_leg17.append(computeEff_leg(1))
        elif j.find('RunBCDEF_mc_leg17_tag8.json') != -1:
            #compute the efficiency for both legs
            effMC_leg17.append(computeEff_leg(0))
            effMC_leg17.append(computeEff_leg(1))

        elif j.find('RunBCDEF_data_leg8.json') != -1:
            #compute the efficiency for both legs
            effData_leg8.append(computeEff_leg(0))
            effData_leg8.append(computeEff_leg(1))
        elif j.find('RunBCDEF_mc_leg8.json') != -1:
            #compute the efficiency for both legs
            effMC_leg8.append(computeEff_leg(0))
            effMC_leg8.append(computeEff_leg(1))


        ###
        # sys
        ###

        elif j.find('RunBCDEF_data_leg17.json') != -1:
            #compute the efficiency for both legs
            effData_leg17_sys.append(computeEff_leg(0))
            effData_leg17_sys.append(computeEff_leg(1))
        elif j.find('RunBCDEF_mc_leg17.json') != -1:
            #compute the efficiency for both legs
            effMC_leg17_sys.append(computeEff_leg(0))
            effMC_leg17_sys.append(computeEff_leg(1))


    ###
    # get efficiency per event
    ###

    effData_event = computeEventEff_fromleg(effData_leg8,effData_leg17)
    effMC_event = computeEventEff_fromleg(effMC_leg8,effMC_leg17)

    # to compute the systematics
    effData_event_sys = computeEventEff_fromleg(effData_leg8,effData_leg17_sys)
    effMC_event_sys = computeEventEff_fromleg(effMC_leg8,effMC_leg17_sys)
    
    ###
    # compute SF using efficiency ratio
    ###

    SF[0] =  (effData_event[0]/effMC_event[0])
    SF[1] = (1-math.sqrt(effData_event[1]**2+ effMC_event[1]**2))*SF[0]
    SF[2] = (1+math.sqrt(effData_event[1]**2+ effMC_event[1]**2))*SF[0]

    print 'SF is', SF

    ###
    # SYStematics
    ###

    SYS[0] = 1
    SYS[1] = (effData_event_sys[0]/effMC_event_sys[0])/SF[0]
    SYS[2] = 1

    print 'SYS is', SYS
