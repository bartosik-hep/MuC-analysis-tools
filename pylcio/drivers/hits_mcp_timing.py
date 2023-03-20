import ROOT as R
import numpy as np
from pyLCIO.drivers.Driver import Driver

from pdb import set_trace as br

CONST_C = R.TMath.C()

def get_oldest_mcp_parent(mcp):
    """Recursively looks for the oldest parent of the MCParticle"""
    pars = mcp.getParents()
    if (len(pars) < 1):
        return mcp
    for par in pars:
        if par is mcp:
            continue
        return get_oldest_mcp_parent(par)

class HitsMCPDriver( Driver ):
    """Driver creating histograms of detector hits and their corresponding MCParticles"""

    HIT_COLLECTION_NAMES = {
        'SimTrackerHit': ['VertexBarrelCollection', 'VertexEndcapCollection', 'InnerTrackerBarrelCollection', 'InnerTrackerEndcapCollection', 'OuterTrackerBarrelCollection', 'OuterTrackerEndcapCollection'],
        'SimCalorimeterHit': ['ECalBarrelCollection', 'ECalEndcapCollection', 'HCalBarrelCollection', 'HCalEndcapCollection']
        # 'SimCalorimeterHit': ['ECalBarrelCollection', 'HCalBarrelCollection']
    }
    HIT_PDGS = [22, 2112]
    PDG_IDS = {
        2112: 5,
        2212: 4,
        -2212: -4,
        321: 3,
        -321: -3,
        -211: -2,
        211: 2,
        -111: -1,
        111: 1,
    }
    
    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.histos = {}
        self.output_path = output_path
    
    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""

        # Book histograms for Tracker hits
        name = 'h_hit_trk_time_mt0_vs_mcp_oldest_time'
        self.histos[name] = R.TH2I(name, ';Hit time - T0 [ns];Oldest MCParticle time [ns]', 220,-20,200, 220,-20,200)
        name = 'h_hit_trk_time_mt0'
        self.histos[name] = R.TH1I(name, ';Hit time - T0 [ns];Tracker hits', 550, -50, 500)
        name = 'h_hit_trk_time'
        self.histos[name] = R.TH1I(name, ';Hit time [ns];Tracker hits', 550, -50, 500)
        name = 'h_hit_trk_mcp_oldest_time'
        self.histos[name] = R.TH1I(name, ';MCParticle time [ns];Tracker oldest MCParticles', 550, -50, 500)
        name = 'h_hit_trk_mcp_pdg'
        self.histos[name] = R.TH1I(name, ';Hit PDG ID;Tracker MCParticles', 2600, -400, 2200)
        name = 'h_hit_trk_mcp_oldest_pdg'
        self.histos[name] = R.TH1I(name, ';Hit PDG ID;Tracker oldest MCParticles', 2600, -400, 2200)
        name = 'h_hit_trk_n_time_e'
        self.histos[name] = R.TH2I( name, ';Hit time - T0 [ns];Neutron energy', 310,-20,600, 1000,0,2)
        name = 'h_hit_trk_time_pdg'
        self.histos[name] = R.TH2I( name, ';Hit time - T0 [ns];MCParticle PdgID', 320,-20,300, 50,-25,25)

        # Book histograms for Calorimeter hits
        name = 'h_hit_cal_time_mt0_vs_mcp_oldest_time'
        self.histos[name] = R.TH2I(name, ';Hit time - T0 [ns];Oldest MCParticle time [ns]', 320,-20,300, 320,-20,300)
        name = 'h_hit_cal_time_maxdiff'
        self.histos[name] = R.TH1I(name, ';t_{min} - t_{max} [ns];Calorimeter hits', 500, 0, 1000)
        name = 'h_hit_cal_time_mt0'
        self.histos[name] = R.TH1I(name, ';Hit time - T0 [ns];Calorimeter hits', 550, -50, 500)
        name = 'h_hit_cal_time'
        self.histos[name] = R.TH1I(name, ';Hit time [ns];Calorimeter hits', 550, -50, 500)
        name = 'h_hit_cal_time_ext'
        self.histos[name] = R.TH1I(name, ';Hit time [ns/];Calorimeter hits', 5005, -50, 50000)
        name = 'h_hit_cal_energy'
        self.histos[name] = R.TH1I(name, ';Hit energy [GeV];Calorimeter hits', 1000, 0, 0.001)
        name = 'h_hit_cal_subenergy'
        self.histos[name] = R.TH1I(name, ';Subhit energy [GeV];Calorimeter subhits', 1000, 0, 0.001)
        name = 'h_hit_cal_mcp_oldest_time'
        self.histos[name] = R.TH1I(name, ';MCParticle time [ns];Calorimeter oldest MCParticles', 550, -50, 500)
        name = 'h_hit_cal_mcp_mult'
        self.histos[name] = R.TH1I(name, ';# MCContributions;Calorimeter hits', 100, 0, 100)
        name = 'h_hit_cal_mcp_pdg'
        self.histos[name] = R.TH1I(name, ';Hit PDG ID;Calorimeter MCParticles', 2600, -400, 2200)
        name = 'h_hit_cal_mcp_oldest_pdg'
        self.histos[name] = R.TH1I(name, ';Hit PDG ID;Calorimeter oldest MCParticles', 2600, -400, 2200)
        name = 'h_hit_cal_n_time_e'
        self.histos[name] = R.TH2I( name, ';Hit time - T0 [ns];Neutron energy', 510,-20,1000, 600,0,1.2)
        name = 'h_hit_cal_time_pdg'
        self.histos[name] = R.TH2I( name, ';Hit time - T0 [ns];MCParticle PdgID', 510,-20,1000, 50,-25,25)
        
        # Create histograms for MCParticles
        for pdg in self.HIT_PDGS:
            for suffix in ['hit', 'nohit']:
                name = 'h_mcp_{0:s}_{1:d}_e'.format(suffix, pdg)
                self.histos[name] = R.TH1I( name, ';Energy [GeV];MCParticles', 5500,0,1.1)
                name = 'h_mcp_{0:s}_{1:d}_time'.format(suffix, pdg)
                self.histos[name] = R.TH1I( name, ';Time [ns];MCParticles', 1200,-20,100)
                name = 'h_mcp_{0:s}_{1:d}_zy'.format(suffix, pdg)
                self.histos[name] = R.TH2I( name, ';Vertex Z [mm];Vertex Y [mm]', 1600,-8000,8000, 400,-2000,2000)
        
        for suffix in ['tlow', 'thigh']:
            name = 'h_mcp_pdg_{0:s}'.format(suffix)
            self.histos[name] = R.TH1I( name, ';Pdg number;MCParticles', 50,-25,25)
            name = 'h_mcp_e_{0:s}'.format(suffix)
            self.histos[name] = R.TH1I( name, ';Energy [GeV];MCParticles', 2000,0,2)
    
    def processEvent( self, event ):
        """Called by the event loop for each event"""
        
        # Get the MCParticle collection from the event
        mcParticles = event.getMcParticles()
        histos = self.histos
        hitMCParticles = set()

        # Loop over hits
        print('Event: {0:d}'.format(event.getEventNumber()))
        for col_type in ['SimTrackerHit', 'SimCalorimeterHit']:
        # for col_type in ['SimCalorimeterHit']:
            for iCol, col_name in enumerate(self.HIT_COLLECTION_NAMES[col_type]):
                # print('Event: {0:d} Col: {1:s}'.format(event.getEventNumber(), col_name))
                col = event.getCollection(col_name)
                # print('  N elements: {0:d}'.format(col.getNumberOfElements()))
                # Filling the Tracker hit properties
                if col_type == 'SimTrackerHit':
                    for iHit in range(col.getNumberOfElements()):
                        hit = col.getElementAt(iHit)
                        hit_time = hit.getTime()
                        self.histos['h_hit_trk_time'].Fill(hit_time)
                        # Calculating the T0 based on the hit position in ns
                        t0 = hit.getPositionVec().Mag() / (CONST_C / 1e6)
                        self.histos['h_hit_trk_time_mt0'].Fill(hit_time - t0)
                        # Getting the MCParticle of the hit
                        mcp = hit.getMCParticle()
                        self.histos['h_hit_trk_mcp_pdg'].Fill(mcp.getPDG())
                        # Getting the oldest MCParticle of the hit
                        mcp_o = get_oldest_mcp_parent(mcp)
                        hitMCParticles.add(mcp_o.id())
                        self.histos['h_hit_trk_mcp_oldest_pdg'].Fill(mcp_o.getPDG())
                        self.histos['h_hit_trk_mcp_oldest_time'].Fill(mcp_o.getTime())
                        self.histos['h_hit_trk_time_mt0_vs_mcp_oldest_time'].Fill(hit_time - t0, mcp_o.getTime())
                        pdg = mcp_o.getPDG()
                        if pdg in self.PDG_IDS:
                            pdg = self.PDG_IDS[pdg]
                        hit_time_t0 = hit_time - t0
                        self.histos['h_hit_trk_time_pdg'].Fill(hit_time_t0, pdg)
                        if mcp_o.getPDG() == 2112:
                            lv = mcp_o.getLorentzVec()
                            self.histos['h_hit_trk_n_time_e'].Fill(hit_time_t0, lv.P())

                # Filling the Calorimeter hit properties
                elif col_type == 'SimCalorimeterHit':
                    nMcp_max = 0
                    for iHit in range(col.getNumberOfElements()):
                        hit = col.getElementAt(iHit)
                        # Looping over the MCParticles of the hit
                        nMcp = hit.getNMCContributions()
                        # self.histos['h_hit_cal_mcp_mult'].Fill(nMcp)
                        hit_times = np.zeros(nMcp, dtype=np.float32)
                        # print('E: {0:.3e} #: {1:d}'.format(hit.getEnergy(), nMcp))
                        self.histos['h_hit_cal_energy'].Fill(hit.getEnergy())
                        for iM in range(nMcp):
                            mcp = hit.getParticleCont(iM)
                            hit_time = hit.getTimeCont(iM)
                            # print('{0:d}. cell: {1:d}  ene: {2:.3e}  tim: {3:.1f}  mcp: {4:d}  pdg: {5:d}'.format(iM, hit.getCellID0(), hit.getEnergyCont(iM), hit_time, mcp.id(), mcp.getPDG()))
                            # print(mcp)
                            # print(iM, hit.getCellID0(), hit.getEnergyCont(iM), hit_time, mcp)
                            hit_times[iM] = hit_time
                            self.histos['h_hit_cal_time'].Fill(hit_time)
                            self.histos['h_hit_cal_time_ext'].Fill(hit_time)
                            self.histos['h_hit_cal_subenergy'].Fill(hit.getEnergyCont(iM))
                            # Calculating the T0 based on the hit position in ns
                            t0 = hit.getPositionVec().Mag() / (CONST_C / 1e6)
                            hit_time_t0 = hit_time - t0
                            self.histos['h_hit_cal_time_mt0'].Fill(hit_time_t0)
                            self.histos['h_hit_cal_mcp_pdg'].Fill(mcp.getPDG())
                            # Getting the oldest MCParticle of the hit
                            mcp_o = get_oldest_mcp_parent(mcp)
                            hitMCParticles.add(mcp_o.id())
                            self.histos['h_hit_cal_mcp_oldest_pdg'].Fill(mcp_o.getPDG())
                            self.histos['h_hit_cal_mcp_oldest_time'].Fill(mcp_o.getTime())
                            self.histos['h_hit_cal_time_mt0_vs_mcp_oldest_time'].Fill(hit_time - t0, mcp_o.getTime())
                            # Filling the histograms on the calorimeter hit timing
                            suffix = 'tlow' if hit_time_t0 < 10 else 'thigh'
                            name = 'h_mcp_pdg_{0:s}'.format(suffix)
                            pdg = mcp_o.getPDG()
                            if pdg in self.PDG_IDS:
                                pdg = self.PDG_IDS[pdg]
                            self.histos[name].Fill(pdg)
                            self.histos['h_hit_cal_time_pdg'].Fill(hit_time_t0, pdg)
                            if mcp_o.getPDG() == 2112:
                                name = 'h_mcp_e_{0:s}'.format(suffix)
                                lv = mcp_o.getLorentzVec()
                                self.histos[name].Fill(lv.P())
                                self.histos['h_hit_cal_n_time_e'].Fill(hit_time_t0, lv.P())
                        if nMcp > 1:
                            self.histos['h_hit_cal_time_maxdiff'].Fill(hit_times.max() - hit_times.min())

        # Loop over all gen-level MCParticles
        for mcp in mcParticles:
            if mcp.getGeneratorStatus() != 1:
                continue
            pdg = mcp.getPDG()
            if pdg not in self.HIT_PDGS:
                continue
            # Get the four vector of the MCParticle
            lv = mcp.getLorentzVec()
            suffix = 'hit' if mcp.id() in hitMCParticles else 'nohit'
            hname = 'h_mcp_{0:s}_{1:d}_e'.format(suffix, pdg)
            self.histos[hname].Fill(lv.P())
            vtx = mcp.getVertex()
            hname = 'h_mcp_{0:s}_{1:d}_zy'.format(suffix, pdg)
            self.histos[hname].Fill(vtx[2], vtx[1])
            hname = 'h_mcp_{0:s}_{1:d}_time'.format(suffix, pdg)
            self.histos[hname].Fill(mcp.getTime())
        
        print('  Hits from {0:d} particles').format(len(hitMCParticles))

    def endOfData( self ):
        """Called by the event loop at the end of the loop"""
        
        # Storing histograms to the output ROOT file
        if self.output_path is not None:
            out_file = R.TFile(self.output_path, 'RECREATE')
            for hname, histo in self.histos.iteritems():
                histo.Write(hname)
            out_file.Close()
