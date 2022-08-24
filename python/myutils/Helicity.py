#!/usr/bin/env python
import ROOT
from ROOT import TVector3, TLorentzVector
from math import pi, sqrt, cos, sin, sinh, log, cosh, acos, acosh
import numpy as np
import array
import os
from BranchTools import Collection
from BranchTools import AddCollectionsModule

class Helicity(AddCollectionsModule):

    def __init__(self, branchName='Helicity'):
        super(Helicity, self).__init__()
        self.branchName = branchName


    def customInit(self, initVars):
        self.sample = initVars['sample']
        self.config = initVars['config']
        self.NbEvents = 0
           

        #self.bidx0 = self.config.get('General', 'btagidx0')
        #self.bidx1 = self.config.get('General', 'btagidx1')

 
        self.addBranch(self.branchName + '_thrust')
        self.addBranch(self.branchName + '_Zh_beam_angle')


    def processEvent(self, tree):
        # if current entry has not been processed yet
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
            self._b(self.branchName + '_Zh_beam_angle')[0] = self.getZHbeamAngle(tree)
            self._b(self.branchName + '_thrust')[0] = self.getThrust(tree)
            
            self.NbEvents+=1
            print("Event ", self.NbEvents, " out of ", tree.GetEntries())


    def getThrust(self, tree):
        lep1, lep2, Vtype = self.getLeptons(tree)
        
        beam = TLorentzVector()

        tmp_lep1, tmp_lep2, tmp_b1, tmp_b2 = TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector()

        tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
        tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
      
        
        idx0 = tree.hJidx[0]       
        idx1 = tree.hJidx[1]       

        if idx0 < 0 or idx1 < 0:
            return -1
 
 
        tmp_b1.SetPtEtaPhiM(tree.Jet_pt[idx0],tree.Jet_eta[idx0],tree.Jet_phi[idx0],tree.Jet_mass[idx0])
        tmp_b2.SetPtEtaPhiM(tree.Jet_pt[idx1],tree.Jet_eta[idx1],tree.Jet_phi[idx1],tree.Jet_mass[idx1])


        if(lep1.Eta()<-10 or lep2.Eta()<-10 or tree.Jet_eta[idx0]<-10 or tree.Jet_eta[idx1]<-10 ):
            return -1

        v_l1 = tmp_lep1.Vect()
        v_l2 = tmp_lep2.Vect()
        v_b1 = tmp_b1.Vect()
        v_b2 = tmp_b2.Vect()

        v_array = [v_l1, v_l2, v_b1, v_b2]

        
        thrust = float('nan')


        directions = []

        for theta in np.arange(0, np.pi, 0.5):
            for phi in np.arange(0, 2*np.pi, 0.5):
                unit_vec = TVector3(0,0,1)  
                unit_vec.SetMag(1.0)
                unit_vec.SetTheta(theta)
                unit_vec.SetPhi(phi)
                directions.append(unit_vec)

        #directions = np.array(directions) 

        scalar_product = [sum([v.Dot(direction)*np.heaviside(v.Dot(direction), 1.0) for v in v_array]) for direction in directions ]  


        d = max(scalar_product)
        
        pt_tot = sum([v.Mag() for v in v_array])        
        
        if pt_tot > 0:
            d = d/pt_tot
        else:
            return -1

        return d


    def getZHbeamAngle(self, tree):
        lep1, lep2, Vtype = self.getLeptons(tree)
        H = self.getH(tree)
        
        beam = TLorentzVector()

        tmp_lep1, tmp_lep2, tmp_H = TLorentzVector(), TLorentzVector(), TLorentzVector()

        tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
        tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
        tmp_H.SetPtEtaPhiM(H.Pt(),H.Eta(),H.Phi(),H.M())

        if(lep1.Eta()<-10 or lep2.Eta()<-10 or tmp_H.Eta()<-10):
            return -100

        beam.SetPxPyPzE(0,0,6500,6500)

        VH = TLorentzVector()
        VH = tmp_lep1+tmp_lep2+tmp_H

        ZHbeamAngle = float('nan')

        try:
            #Theta  = acos((V_mom.Vect().Unit()).Dot(beam.Vect().Unit()))
            ZHbeamAngle = (VH.Vect().Unit()).Dot(beam.Vect().Unit())
        except Exception:
            #pass
            ZHbeamAngle = -100

        return ZHbeamAngle



    def getLeptons(self, tree):
        isZee = tree.isZee
        isZmm = tree.isZmm
        vLidx = tree.vLidx
        Vtype = tree.Vtype

        #basic check for consistency
        #if (isZee+isZmm != 1): print "isZee and isZmm inconsistent"
        if (isZee == 1 and Vtype != 1): print "isZee and Vtype inconsistent"
        if (isZmm == 1 and Vtype != 0): print "isZmm and Vtype inconsistent"

        Electron_pt = tree.Electron_pt
        Electron_eta = tree.Electron_eta
        Electron_phi = tree.Electron_phi
        Electron_mass = tree.Electron_mass
        Electron_charge = tree.Electron_charge

        Muon_pt = tree.Muon_pt
        Muon_eta = tree.Muon_eta
        Muon_phi = tree.Muon_phi
        Muon_mass = tree.Muon_mass
        Muon_charge = tree.Muon_charge

        lep1, lep2 = TLorentzVector(), TLorentzVector()

        if (Vtype == 0): #Zmm
            order = Muon_charge[vLidx[0]] > Muon_charge[vLidx[1]]
            try:
                if order:
                    lep1.SetPtEtaPhiM(Muon_pt[vLidx[0]],Muon_eta[vLidx[0]],Muon_phi[vLidx[0]],Muon_mass[vLidx[0]])
                else:
                    lep2.SetPtEtaPhiM(Muon_pt[vLidx[0]],Muon_eta[vLidx[0]],Muon_phi[vLidx[0]],Muon_mass[vLidx[0]])
            except:
                if order:
                    lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
                else:
                    lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)                  
            try:
                if order:
                    lep2.SetPtEtaPhiM(Muon_pt[vLidx[1]],Muon_eta[vLidx[1]],Muon_phi[vLidx[1]],Muon_mass[vLidx[1]])
                else:
                    lep1.SetPtEtaPhiM(Muon_pt[vLidx[1]],Muon_eta[vLidx[1]],Muon_phi[vLidx[1]],Muon_mass[vLidx[1]])
            except:
                if order:
                    lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
                else:
                    lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)



        if (Vtype == 1): #Zee
            order = Electron_charge[vLidx[0]] > Electron_charge[vLidx[1]]
            try:
                if order:
                    lep1.SetPtEtaPhiM(Electron_pt[vLidx[0]],Electron_eta[vLidx[0]],Electron_phi[vLidx[0]],Electron_mass[vLidx[0]])
                else:
                    lep2.SetPtEtaPhiM(Electron_pt[vLidx[0]],Electron_eta[vLidx[0]],Electron_phi[vLidx[0]],Electron_mass[vLidx[0]])
            except:
                if order:
                    lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
                else:
                    lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
            try:
                if order:
                    lep2.SetPtEtaPhiM(Electron_pt[vLidx[1]],Electron_eta[vLidx[1]],Electron_phi[vLidx[1]],Electron_mass[vLidx[1]])
                else:
                    lep1.SetPtEtaPhiM(Electron_pt[vLidx[1]],Electron_eta[vLidx[1]],Electron_phi[vLidx[1]],Electron_mass[vLidx[1]])
            except:
                if order:
                    lep2.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)
                else:
                    lep1.SetPtEtaPhiM(-1.0,-1.0,-1.0,-1.0)

        return lep1, lep2, Vtype

    def getH(self, tree):
        H_pt = tree.H_pt
        H_eta = tree.H_eta
        H_phi = tree.H_phi
        H_mass = tree.H_mass
        H = TLorentzVector()
        H.SetPtEtaPhiM(H_pt, H_eta, H_phi, H_mass)
        return H

