#!/usr/bin/env python
import ROOT
from ROOT import TVector3, TLorentzVector
from math import pi, sqrt, cos, sin, sinh, log, cosh, acos, acosh
import numpy as np
import array
import os
from BranchTools import Collection
from BranchTools import AddCollectionsModule

class EFT_obs(AddCollectionsModule):

    def __init__(self, branchName='EFT_obs'):
        super(EFT_obs, self).__init__()
        self.branchName = branchName


    def customInit(self, initVars):
        self.sample = initVars['sample']
        self.config = initVars['config']
            
        self.addBranch(self.branchName + '_Theta_e')
        self.addBranch(self.branchName + '_Theta_m')
        self.addBranch(self.branchName + '_Theta_l')
        self.addBranch(self.branchName + '_theta_e')
        self.addBranch(self.branchName + '_theta_m')
        self.addBranch(self.branchName + '_theta_l')
        self.addBranch(self.branchName + '_phi_e')
        self.addBranch(self.branchName + '_phi_m')
        self.addBranch(self.branchName + '_phi_l')
        self.addBranch(self.branchName + '_phi_weight')


        self.addBranch(self.branchName + '_LHE_Theta_e')
        self.addBranch(self.branchName + '_LHE_Theta_m')
        self.addBranch(self.branchName + '_LHE_Theta_l')
        self.addBranch(self.branchName + '_LHE_theta_e')
        self.addBranch(self.branchName + '_LHE_theta_m')
        self.addBranch(self.branchName + '_LHE_theta_l')
        self.addBranch(self.branchName + '_LHE_phi_e')
        self.addBranch(self.branchName + '_LHE_phi_m')
        self.addBranch(self.branchName + '_LHE_phi_l')
        self.addBranch(self.branchName + '_LHE_phi_weight')


        #debugging_only
        self.addBranch(self.branchName + '_LHE_H_mass_from_bb')
        self.addBranch(self.branchName + '_LHE_H_pt_from_bb')
        self.addBranch(self.branchName + '_LHE_H_eta_from_bb')
        self.addBranch(self.branchName + '_LHE_H_phi_from_bb')
        self.addBranch(self.branchName + '_LHE_Vtype')

        self.addBranch(self.branchName + '_VH_mass')
        self.addBranch(self.branchName + '_LHE_VH_mass_from_llbb')

        self.sampleTree = initVars['sampleTree']

    def processEvent(self, tree):
        # if current entry has not been processed yet
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
            self._b(self.branchName + '_Theta_e')[0] = -99997.0
            self._b(self.branchName + '_Theta_m')[0] = -99997.0
            self._b(self.branchName + '_Theta_l')[0] = -99997.0
            self._b(self.branchName + '_theta_e')[0] = -99997.0
            self._b(self.branchName + '_theta_m')[0] = -99997.0
            self._b(self.branchName + '_theta_l')[0] = -99997.0
            self._b(self.branchName + '_phi_e')[0]   = -99997.0
            self._b(self.branchName + '_phi_m')[0]   = -99997.0
            self._b(self.branchName + '_phi_l')[0]   = -99997.0
            self._b(self.branchName + '_phi_weight')[0]   = -99997.0


            self._b(self.branchName + '_LHE_Theta_e')[0] = -99997.0
            self._b(self.branchName + '_LHE_Theta_m')[0] = -99997.0
            self._b(self.branchName + '_LHE_Theta_l')[0] = -99997.0
            self._b(self.branchName + '_LHE_theta_e')[0] = -99997.0
            self._b(self.branchName + '_LHE_theta_m')[0] = -99997.0
            self._b(self.branchName + '_LHE_theta_l')[0] = -99997.0
            self._b(self.branchName + '_LHE_phi_e')[0]   = -99997.0
            self._b(self.branchName + '_LHE_phi_m')[0]   = -99997.0
            self._b(self.branchName + '_LHE_phi_l')[0]   = -99997.0
            self._b(self.branchName + '_LHE_phi_weight')[0]   = -99997.0

            self._b(self.branchName + '_LHE_H_pt_from_bb')[0]   = -99997.0
            self._b(self.branchName + '_LHE_H_eta_from_bb')[0]   = -99997.0
            self._b(self.branchName + '_LHE_H_phi_from_bb')[0]   = -99997.0
            self._b(self.branchName + '_LHE_H_mass_from_bb')[0]   = -99997.0
            self._b(self.branchName + '_LHE_Vtype')[0]   = -99997.0
            
            self._b(self.branchName + '_VH_mass')[0]   = -99997.0
            self._b(self.branchName + '_LHE_VH_mass_from_llbb')[0]   = -99997.0


            self._b(self.branchName + '_LHE_H_pt_from_bb')[0],self._b(self.branchName + '_LHE_H_eta_from_bb')[0],self._b(self.branchName + '_LHE_H_phi_from_bb')[0],self._b(self.branchName + '_LHE_H_mass_from_bb')[0]  = self.getHiggsPtEtaPhiMFrombb(tree)
            self._b(self.branchName + '_LHE_Vtype')[0] = self.getLHELeptons(tree)[2]

            self._b(self.branchName + '_Theta_e')[0], self._b(self.branchName + '_Theta_m')[0], self._b(self.branchName + '_Theta_l')[0] = self.getTheta(tree,False)
            self._b(self.branchName + '_theta_e')[0], self._b(self.branchName + '_theta_m')[0], self._b(self.branchName + '_theta_l')[0] = self.gettheta(tree,False)
            self._b(self.branchName + '_phi_e')[0], self._b(self.branchName + '_phi_m')[0], self._b(self.branchName + '_phi_l')[0] = self.getphi(tree,False)

            self._b(self.branchName + '_LHE_Theta_e')[0], self._b(self.branchName + '_LHE_Theta_m')[0], self._b(self.branchName + '_LHE_Theta_l')[0] = self.getTheta(tree,True)
            self._b(self.branchName + '_LHE_theta_e')[0], self._b(self.branchName + '_LHE_theta_m')[0], self._b(self.branchName + '_LHE_theta_l')[0] = self.gettheta(tree,True)
            self._b(self.branchName + '_LHE_phi_e')[0], self._b(self.branchName + '_LHE_phi_m')[0], self._b(self.branchName + '_LHE_phi_l')[0] = self.getphi(tree,True)

            self._b(self.branchName + '_VH_mass')[0] = self.getVH(tree)[3]
            self._b(self.branchName + '_LHE_VH_mass_from_llbb')[0] = self.getLHEVH(tree)[3]

            self._b(self.branchName + '_phi_weight')[0] = self.getphiweight(tree,False)
            self._b(self.branchName + '_LHE_phi_weight')[0] = self.getphiweight(tree,True)


    def getLeptons(self, tree):
        isZee = tree.isZee
        isZmm = tree.isZmm
        vLidx = tree.vLidx
        Vtype = tree.Vtype

        #basic check for consistency
        if (isZee+isZmm != 1): print "isZee and isZmm inconsistent"
        if (isZee == 1 and Vtype != 1): print "isZee and Vtype inconsistent"
        if (isZmm == 1 and Vtype != 0): print "isZmm and Vtype inconsistent"

        Electron_pt = tree.Electron_pt
        Electron_eta = tree.Electron_eta
        Electron_phi = tree.Electron_phi
        Electron_mass = tree.Electron_mass

        Muon_pt = tree.Muon_pt
        Muon_eta = tree.Muon_eta
        Muon_phi = tree.Muon_phi
        Muon_mass = tree.Muon_mass
        lep1, lep2 = TLorentzVector(), TLorentzVector()

        if (Vtype == 0): #Zmm
            try:
                lep1.SetPtEtaPhiM(Muon_pt[vLidx[0]],Muon_eta[vLidx[0]],Muon_phi[vLidx[0]],Muon_mass[vLidx[0]])
            except:
                lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
            try:
                lep2.SetPtEtaPhiM(Muon_pt[vLidx[1]],Muon_eta[vLidx[1]],Muon_phi[vLidx[1]],Muon_mass[vLidx[1]])
            except:
                lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)

        if (Vtype == 1): #Zee
            try:
                lep1.SetPtEtaPhiM(Electron_pt[vLidx[0]],Electron_eta[vLidx[0]],Electron_phi[vLidx[0]],Electron_mass[vLidx[0]])
            except:
                lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
            try:
                lep2.SetPtEtaPhiM(Electron_pt[vLidx[1]],Electron_eta[vLidx[1]],Electron_phi[vLidx[1]],Electron_mass[vLidx[1]])
            except:
                lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)

        return lep1, lep2, Vtype

    def getLHELeptons(self,tree):
        if tree.isData == 1:
            lep1,lep2 = TLorentzVector(),TLorentzVector()
            lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
            lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
            LHE_Vtype = -1
            return lep1, lep2, LHE_Vtype


        LHE_pdgId = list(tree.LHEPart_pdgId)
        LHE_pt = list(tree.LHEPart_pt)
        LHE_eta = list(tree.LHEPart_eta) 
        LHE_phi = list(tree.LHEPart_phi)
        LHE_mass = list(tree.LHEPart_mass)
        LHE_isZee = False
        LHE_isZmm = False
	LHE_Vtype = -1

        lep1, lep2 = TLorentzVector(), TLorentzVector()


        if ((11 in LHE_pdgId) and (-11 in LHE_pdgId)): LHE_isZee = True
        if ((13 in LHE_pdgId) and (-13 in LHE_pdgId)): LHE_isZmm = True

        #if abs(LHE_pdgId) == 11: LHE_isZee = True
        #if abs(LHE_pdgId) == 13: LHE_isZmm = True

	if (LHE_isZee ^ LHE_isZmm):
            if LHE_isZee: 
                LHE_Vtype = 1
            else:
                pass
            if LHE_isZmm: 
                LHE_Vtype = 0
            else:
                pass
        elif (LHE_isZee and LHE_isZmm): 
            LHE_Vtype = -9 
        else:
            LHE_Vtype = -99


	if (LHE_Vtype == 1):
            order = LHE_pt[LHE_pdgId.index(11)] > LHE_pt[LHE_pdgId.index(-11)]
            try:
                if order:
                    lep1.SetPtEtaPhiM(LHE_pt[LHE_pdgId.index(11)],LHE_eta[LHE_pdgId.index(11)],LHE_phi[LHE_pdgId.index(11)],LHE_mass[LHE_pdgId.index(11)])
                else:
                    lep2.SetPtEtaPhiM(LHE_pt[LHE_pdgId.index(11)],LHE_eta[LHE_pdgId.index(11)],LHE_phi[LHE_pdgId.index(11)],LHE_mass[LHE_pdgId.index(11)])
            except:
                if order:
                    lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
                else:
                    lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)

            try:
                if order:
                    lep2.SetPtEtaPhiM(LHE_pt[LHE_pdgId.index(-11)],LHE_eta[LHE_pdgId.index(-11)],LHE_phi[LHE_pdgId.index(-11)],LHE_mass[LHE_pdgId.index(-11)])
                else:
                    lep1.SetPtEtaPhiM(LHE_pt[LHE_pdgId.index(-11)],LHE_eta[LHE_pdgId.index(-11)],LHE_phi[LHE_pdgId.index(-11)],LHE_mass[LHE_pdgId.index(-11)])
            except:
                if order:
                    lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
                else:
                    lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)


        if (LHE_Vtype == 0):
            order = LHE_pt[LHE_pdgId.index(13)] > LHE_pt[LHE_pdgId.index(-13)]
            try:
                if order:
                    lep1.SetPtEtaPhiM(LHE_pt[LHE_pdgId.index(13)],LHE_eta[LHE_pdgId.index(13)],LHE_phi[LHE_pdgId.index(13)],LHE_mass[LHE_pdgId.index(13)])
                else:
                    lep2.SetPtEtaPhiM(LHE_pt[LHE_pdgId.index(13)],LHE_eta[LHE_pdgId.index(13)],LHE_phi[LHE_pdgId.index(13)],LHE_mass[LHE_pdgId.index(13)])
            except:
                if order:
                    lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
                else:
                    lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)


            try:
                if order:
                    lep2.SetPtEtaPhiM(LHE_pt[LHE_pdgId.index(-13)],LHE_eta[LHE_pdgId.index(-13)],LHE_phi[LHE_pdgId.index(-13)],LHE_mass[LHE_pdgId.index(-13)])
                else:
                    lep1.SetPtEtaPhiM(LHE_pt[LHE_pdgId.index(-13)],LHE_eta[LHE_pdgId.index(-13)],LHE_phi[LHE_pdgId.index(-13)],LHE_mass[LHE_pdgId.index(-13)])
            except:
                if order:
                    lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
                else:
                    lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)


	return lep1, lep2, LHE_Vtype


    



    def getH(self, tree):
        H_pt = tree.H_pt
        H_eta = tree.H_eta
        H_phi = tree.H_phi
        H_mass = tree.H_mass
        H = TLorentzVector()
        H.SetPtEtaPhiM(H_pt, H_eta, H_phi, H_mass)
        return H


    def getLHEH(self,tree): 
        if tree.isData == 1:
            bb = TLorentzVector()
            bb.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
            return bb
            
        LHE_pdgId = list(tree.LHEPart_pdgId)
        LHE_pt = list(tree.LHEPart_pt)
        LHE_eta = list(tree.LHEPart_eta)
        LHE_phi = list(tree.LHEPart_phi)
        LHE_mass = list(tree.LHEPart_mass)

	b1, b2, bb = TLorentzVector(), TLorentzVector(), TLorentzVector()
        bb.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)

	try:
            b1.SetPtEtaPhiM(LHE_pt[LHE_pdgId.index(-5)],LHE_eta[LHE_pdgId.index(-5)],LHE_phi[LHE_pdgId.index(-5)],LHE_mass[LHE_pdgId.index(-5)])
        except:
            b1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
        try:
            b2.SetPtEtaPhiM(LHE_pt[LHE_pdgId.index(5)],LHE_eta[LHE_pdgId.index(5)],LHE_phi[LHE_pdgId.index(5)],LHE_mass[LHE_pdgId.index(5)])
        except:
            b2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)

        bb = b1 + b2
        return bb


    def getLHEVH(self,tree):
        H = self.getLHEH(tree)
        lep1, lep2, LHE_Vtype = self.getLHELeptons(tree)
        V, VH = TLorentzVector(), TLorentzVector()
        V = lep1 + lep2
        VH = V + H
        return VH.Pt(), VH.Eta(), VH.Phi(), VH.M()


    def getVH(self,tree):
        H = self.getH(tree)
        V_pt = tree.V_pt 
        V_eta = tree.V_eta
        V_phi = tree.V_phi
        V_mass = tree.V_mass
        V, VH = TLorentzVector(), TLorentzVector()
        try:
            V.SetPtEtaPhiM(V_pt, V_eta, V_phi, V_mass)
        except:
            V.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
        VH = V + H 
        return VH.Pt(), VH.Eta(), VH.Phi(), VH.M()


    def getHiggsPtEtaPhiMFrombb(self,tree):
        H = self.getLHEH(tree)
        H_mass = H.M()
        H_eta = H.Eta()
        H_phi = H.Phi()
        H_pt = H.Pt()
        return H_pt, H_eta, H_phi, H_mass

    def getTheta(self, tree, LHE_level):
        if LHE_level:
            lep1, lep2, Vtype = self.getLHELeptons(tree)
            H = self.getLHEH(tree)
        else:
            lep1, lep2, Vtype = self.getLeptons(tree)
            H = self.getH(tree)
        
        beam = TLorentzVector()

        tmp_lep1, tmp_lep2, tmp_H = TLorentzVector(), TLorentzVector(), TLorentzVector()

        tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
        tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
        tmp_H.SetPtEtaPhiM(H.Pt(),H.Eta(),H.Phi(),H.M())

        if(lep1.Eta()<-10 or lep2.Eta()<-10 or tmp_H.Eta()<-10):
            return -100, -100, -100

        beam.SetPxPyPzE(0,0,6500,6500)

        V_mom, bVH = TLorentzVector(), TVector3()
        V_mom = tmp_lep1+tmp_lep2
        bVH = (tmp_lep1+tmp_lep2+tmp_H).BoostVector()

        V_mom.Boost(-bVH)

        Theta = float('nan')

        try:
            #Theta  = acos((V_mom.Vect().Unit()).Dot(beam.Vect().Unit()))
            Theta = (V_mom.Vect().Unit()).Angle(beam.Vect().Unit())
        except Exception:
            pass
            Theta = -100

        if (Vtype == 0):
            Theta_l = Theta
            Theta_m = Theta
            Theta_e = -99999

        elif (Vtype == 1):
            Theta_l = Theta
            Theta_e = Theta
            Theta_m = -99999

        else:
            Theta_l = -9999
            Theta_m = -9999
            Theta_e = -9999

        return Theta_e, Theta_m, Theta_l


    def gettheta(self, tree, LHE_level):
        if LHE_level:
            lep1, lep2, Vtype = self.getLHELeptons(tree)
            H = self.getLHEH(tree)
        else:
            lep1, lep2, Vtype = self.getLeptons(tree)
            H = self.getH(tree)

        tmp_lep1, tmp_lep2, tmp_H = TLorentzVector(), TLorentzVector(), TLorentzVector()

        tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
        tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
        tmp_H.SetPtEtaPhiM(H.Pt(),H.Eta(),H.Phi(),H.M())

        if(lep1.Eta()<-10 or lep2.Eta()<-10 or tmp_H.Eta()<-10):
            return -100, -100, -100

        V_mom, bVH, bV = TLorentzVector(), TVector3(), TVector3()

        bVH = (tmp_lep1 + tmp_lep2 + tmp_H).BoostVector()
        V_mom = (tmp_lep1 + tmp_lep2)

        V_mom.Boost(-bVH)
        tmp_lep1.Boost(-bVH)

        bV = V_mom.BoostVector()
        tmp_lep1.Boost(-bV)

        theta = float('nan')
        try:
            theta = (V_mom).Angle(tmp_lep1.Vect())
        except Exception:
            pass
            theta = -100

        if (Vtype == 0):
            theta_l = theta
            theta_m = theta
            theta_e = -99999

        elif (Vtype == 1):
            theta_l = theta
            theta_e = theta
            theta_m = -99999

        else:
            theta_l = -9999
            theta_m = -9999
            theta_e = -9999

        return theta_e, theta_m, theta_l


    def getphi(self, tree, LHE_level):
        if LHE_level:
            lep1, lep2, Vtype = self.getLHELeptons(tree)
            H = self.getLHEH(tree)
        else:
            lep1, lep2, Vtype = self.getLeptons(tree)
            H = self.getH(tree)


        beam = TLorentzVector()

        tmp_lep1, tmp_lep2, tmp_H = TLorentzVector(), TLorentzVector(), TLorentzVector()

        tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
        tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
        tmp_H.SetPtEtaPhiM(H.Pt(),H.Eta(),H.Phi(),H.M())

        if(lep1.Eta()<-10 or lep2.Eta()<-10 or tmp_H.Eta()<-10):
            return -100, -100, -100

        beam.SetPxPyPzE(0,0,6500,6500)

        V_mom, bVH, n_scatter, n_decay = TLorentzVector(), TVector3(), TVector3(), TVector3()
        bVH = (tmp_lep1+tmp_lep2+tmp_H).BoostVector()
        V_mom = tmp_lep1+tmp_lep2

        tmp_lep1.Boost(-bVH)
        tmp_lep2.Boost(-bVH)
        V_mom.Boost(-bVH)
        #beam.Boost(-bVH)

        n_scatter = ((beam.Vect().Unit()).Cross(V_mom.Vect())).Unit()
        n_decay   = (tmp_lep1.Vect().Cross(tmp_lep2.Vect())).Unit()

        sign_flip =  -1 if ( ((n_scatter.Cross(n_decay))*(V_mom.Vect())) < 0 ) else +1

        try:
            phi = sign_flip*acos(n_scatter.Dot(n_decay))
        except Exception:
            pass
            phi = -100

        if (Vtype == 0):
            phi_l = phi
            phi_m = phi
            phi_e = -99999

        elif (Vtype == 1):
            phi_l = phi
            phi_e = phi
            phi_m = -99999

        else:
            phi_l = -9999
            phi_m = -9999
            phi_e = -9999

        return phi_e, phi_m, phi_l


    def getphiweight(self, tree, LHE_level):
        Theta = self.getTheta(tree, LHE_level)[2]
        theta = self.gettheta(tree, LHE_level)[2]
        try:
            weight = np.sin(2*theta) * np.sin(2*Theta)
        except Exception:
            pass 
            weight = -999
        
        return weight


