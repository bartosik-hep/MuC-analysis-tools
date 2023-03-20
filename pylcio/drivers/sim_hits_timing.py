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

class HitsTimingDriver( Driver ):
    """Driver creating histograms of detector hits timing and energy"""

    HIT_COLLECTION_NAMES = {
        'trk': ['VertexBarrelCollection', 'VertexEndcapCollection', 'InnerTrackerBarrelCollection', 'InnerTrackerEndcapCollection', 'OuterTrackerBarrelCollection', 'OuterTrackerEndcapCollection'],
        'cal': ['ECalBarrelCollection', 'ECalEndcapCollection', 'HCalBarrelCollection', 'HCalEndcapCollection'],
        'muo': ['YokeBarrelCollection', 'YokeEndcapCollection']
    }
    
    TIME_CUTS = [100, 10, 5, 2]
    
    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.histos = {}
        self.output_path = output_path
    
    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""

        for hit_name in ['trk', 'cal', 'muo']:
            histos = {}
            name = 'hit_time'
            histos[name] = R.TH1I('_'.join([hit_name, name]), ';Hit time [ns];Hits', 1100, -50, 500)
            name = 'hit_time_mt0'
            histos[name] = R.TH1I('_'.join([hit_name, name]), ';Hit time - T0 [ns];Hits', 1100, -50, 500)
            name = 'hit_time_mt0_e'
            histos[name] = R.TH2I('_'.join([hit_name, name]), ';Hit time - T0 [ns];Hit energy [MeV]', 1100,-50,500, 1000, 0, 2)
            for time_cut in self.TIME_CUTS:
                name = 'hit_e_tlt{0:d}'.format(time_cut)
                histos[name] = R.TH1I('_'.join([hit_name, name]), ';Hit energy [MeV];Hits', 10000, 0, 20)
                name = 'hit_zy_tlt{0:d}'.format(time_cut)
                histos[name] = R.TH2I('_'.join([hit_name, name]), ';Hit Z [mm];Hit Y [mm]', 2400,-6000,6000, 2400,-6000,6000)
            self.histos[hit_name] = histos
    
    def processEvent( self, event ):
        """Called by the event loop for each event"""

        # Loop over hits
        print('Event: {0:d}'.format(event.getEventNumber()))
        for col_type in ['trk', 'cal', 'muo']:
            histos = self.histos[col_type]
            for iCol, col_name in enumerate(self.HIT_COLLECTION_NAMES[col_type]):
                # print('Event: {0:d} Col: {1:s}'.format(event.getEventNumber(), col_name))
                col = event.getCollection(col_name)
                # Filling the Tracker hit properties
                if col_type == 'trk':
                    for iHit in range(col.getNumberOfElements()):
                        hit = col.getElementAt(iHit)
                        hit_time = hit.getTime()
                        histos['hit_time'].Fill(hit_time)
                        # Calculating the T0 based on the hit position in ns
                        hit_pos_vec = hit.getPositionVec()
                        hit_time_mt0 = hit_time - (hit_pos_vec.Mag() / (CONST_C / 1e6))
                        histos['hit_time_mt0'].Fill(hit_time_mt0)
                        histos['hit_time_mt0_e'].Fill(hit_time_mt0, hit.getEDep()*1e3)

                        for time_cut in self.TIME_CUTS:
                            if hit_time_mt0 > time_cut:
                                break
                            name = 'hit_e_tlt{0:d}'.format(time_cut)
                            histos[name].Fill(hit.getEDep()*1e3)
                            name = 'hit_zy_tlt{0:d}'.format(time_cut)
                            histos[name].Fill(hit_pos_vec.Z(), hit_pos_vec.Y())

                # Filling the Calorimeter-type hit properties
                elif col_type == 'cal' or col_type == 'muo':
                    for iHit in range(col.getNumberOfElements()):
                        hit = col.getElementAt(iHit)
                        hit_pos_vec = hit.getPositionVec()
                        hit_t0 = hit_pos_vec.Mag() / (CONST_C / 1e6)
                        nSubhits = hit.getNMCContributions()
                        # Hit from multiple particles (summing all particles in the time window)
                        if nSubhits > 1:
                            hit_time = 1e9
                            hit_e = {time_cut: 0.0 for time_cut in self.TIME_CUTS}
                            for iSubhit in range(nSubhits):
                                subhit_time = hit.getTimeCont(iSubhit)
                                if abs(subhit_time - hit_t0) < abs(hit_time - hit_t0):
                                    hit_time = subhit_time
                                for time_cut in self.TIME_CUTS:
                                    if subhit_time - hit_t0 > time_cut:
                                        break
                                    hit_e[time_cut] += hit.getEnergyCont(iSubhit)
                            histos['hit_time'].Fill(hit_time)
                            histos['hit_time_mt0'].Fill(hit_time - hit_t0)
                            histos['hit_time_mt0_e'].Fill(hit_time - hit_t0, hit_e[10]*1e3)
                            for time_cut in self.TIME_CUTS:
                                e = hit_e[time_cut]
                                if e > 0.0:
                                    name = 'hit_e_tlt{0:d}'.format(time_cut)
                                    histos[name].Fill(e*1e3)
                                    name = 'hit_zy_tlt{0:d}'.format(time_cut)
                                    histos[name].Fill(hit_pos_vec.Z(), hit_pos_vec.Y())
                        # All hit from a single particle (using the whole hit directly)
                        else:
                            hit_time = hit.getTimeCont(0)
                            hit_time_mt0 = hit_time - hit_t0
                            histos['hit_time'].Fill(hit_time)
                            histos['hit_time_mt0'].Fill(hit_time_mt0)
                            hit_e = hit.getEnergy()
                            histos['hit_time_mt0_e'].Fill(hit_time - hit_t0, hit_e*1e3)
                            for time_cut in self.TIME_CUTS:
                                if hit_time - hit_t0 > time_cut:
                                    break
                                name = 'hit_e_tlt{0:d}'.format(time_cut)
                                histos[name].Fill(hit_e*1e3)
                                name = 'hit_zy_tlt{0:d}'.format(time_cut)
                                histos[name].Fill(hit_pos_vec.Z(), hit_pos_vec.Y())

    def endOfData( self ):
        """Called by the event loop at the end of the loop"""
        
        # Storing histograms to the output ROOT file
        if self.output_path is not None:
            out_file = R.TFile(self.output_path, 'RECREATE')
            for hit_name, histos in self.histos.iteritems():
                for hname, histo in histos.iteritems():
                    histo.Write()
            out_file.Close()
