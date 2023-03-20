import ROOT as R
import numpy as np
from pyLCIO.drivers.Driver import Driver
from pyLCIO import EVENT, UTIL
from utils import get_oldest_mcp_parent, pdg_to_type

from pdb import set_trace as br

CONST_C = R.TMath.C()


class HitPropsDriver( Driver ):
    """Driver creating histograms of detector hits and their corresponding MCParticles"""

    HIT_COLLECTION_NAMES = {
        'SimTrackerHit': ['VertexBarrelCollection', 'VertexEndcapCollection', 'InnerTrackerBarrelCollection', 'InnerTrackerEndcapCollection', 'OuterTrackerBarrelCollection', 'OuterTrackerEndcapCollection'],
        'SimCalorimeterHit': ['ECalBarrelCollection', 'ECalEndcapCollection', 'HCalBarrelCollection', 'HCalEndcapCollection', 'HCalRingCollection']
    }
    N_LAYERS = {
        'SimTrackerHit': [3*2,3*2, 3,7, 3,4],
        'SimCalorimeterHit': [40,40, 60,60,68]
    }
    MCP_TYPES = {
        'SimTrackerHit': (9, 10, 11),
        'SimCalorimeterHit': (1, 2, 3, 4, 7, 8, 9, 10, 11)
    }
    
    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.histos = {}
        self.output_path = output_path
    
    
    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""
        # for iCol, col in enumerate(self.HIT_COLLECTION_NAMES['SimTrackerHit']):
        #     histos = {}
        #     for layer in range(1, self.N_LAYERS['SimTrackerHit'][iCol]+1):
        #         # Book histograms for hit multiplicity vs particle type
        #         name = 'h_hit_mult_vs_mcp_type'
        #         histos[(name, layer)] = R.TH1I('{0}_{1}_{2:d}'.format(name, col, layer), ';Oldest MCParticle type;Hit multiplicity', 40,0,40)
        #         name = 'p_hit_time_vs_mcp_type'
        #         histos[(name, layer)] = R.TProfile('{0}_{1}_{2:d}'.format(name, col, layer), ';Oldest MCParticle type;Hit time [ns]', 40,0,40, 'S')
        #         name = 'p_hit_e_vs_mcp_type'
        #         histos[(name, layer)] = R.TProfile('{0}_{1}_{2:d}'.format(name, col, layer), ';Oldest MCParticle type;Hit energy [KeV]', 40,0,40, 'S')
        #         for mcp_type in self.MCP_TYPES['SimTrackerHit']:
        #             name = 'h_hit_time_mcp_{0:d}'.format(mcp_type)
        #             histos[(name, layer)] = R.TH1I('{0}_{1}_{2:d}'.format(name, col, layer), ';Hit time [ns];Hits', 1100,-1,10)
        #     self.histos[col] = histos
        for iCol, col in enumerate(self.HIT_COLLECTION_NAMES['SimCalorimeterHit']):
            histos = {}
            # Book histograms for hit multiplicity vs particle type
            name = 'h_hit_mult_vs_mcp_type'
            histos[name] = R.TH1I('{0}_{1}'.format(name, col), ';Oldest MCParticle type;Hit multiplicity', 40,0,40)
            name = 'p_hit_time_vs_mcp_type'
            histos[name] = R.TProfile('{0}_{1}'.format(name, col), ';Oldest MCParticle type;Hit time [ns]', 40,0,40, 'S')
            name = 'p_hit_e_vs_mcp_type'
            histos[name] = R.TProfile('{0}_{1}'.format(name, col), ';Oldest MCParticle type;Hit energy [KeV]', 40,0,40, 'S')
            nlayers = self.N_LAYERS['SimCalorimeterHit'][iCol]
            name = 'h2_hit_layer_vs_mcp_type'
            histos[name] = R.TH2I('{0}_{1}'.format(name, col), ';Oldest MCParticle type;Hit layer', 40,0,40, nlayers+1, 0, nlayers+1)
            name = 'h2_hit_layer_vs_hit_time'
            histos[name] = R.TH2I('{0}_{1}'.format(name, col), ';Hit time [ns] type;Hit layer', 1100,-50,500, nlayers+1, 0, nlayers+1)
            name = 'h2_hit_e_vs_hit_time'
            histos[name] = R.TH2I('{0}_{1}'.format(name, col), ';Hit time [ns] type;Hit energy [KeV]', 1100,-50,500, 500,0,500)
            name = 'h2_hit_layer_vs_hit_time_zoom'
            histos[name] = R.TH2I('{0}_{1}'.format(name, col), ';Hit time [ns] type;Hit layer', 5000,-10,40, nlayers+1, 0, nlayers+1)
            name = 'h_hit_e'
            histos[name] = R.TH1I('{0}_{1}'.format(name, col), ';Hit energy [KeV];Hits', 10000,0,1000)
            name = 'h_hit_time'
            histos[name] = R.TH1I('{0}_{1}'.format(name, col), ';Hit time [ns];Hits', 2200,-10,40)
            for mcp_type in self.MCP_TYPES['SimCalorimeterHit']:
                    name = 'h_hit_time_mcp_{0:d}'.format(mcp_type)
                    histos[name] = R.TH1I('{0}_{1}'.format(name, col), ';Hit time [ns];Hits', 2200,-10,40)
                    name = 'h_hit_e_mcp_{0:d}'.format(mcp_type)
                    histos[name] = R.TH1I('{0}_{1}'.format(name, col), ';Hit energy [KeV];Hits', 10000,0,1000)
            self.histos[col] = histos
    
    
    def processEvent( self, event ):
        """Called by the event loop for each event"""
        
        # Get the MCParticle collection from the event
        mcParticles = event.getMcParticles()
        histos = self.histos
        hitMCParticles = set()

        # Loop over hits
        print('Event: {0:d}'.format(event.getEventNumber()))
        # for col_type in ['SimTrackerHit', 'SimCalorimeterHit']:
        for col_type in ['SimCalorimeterHit']:
            for iCol, col_name in enumerate(self.HIT_COLLECTION_NAMES[col_type]):
                # print('Event: {0:d} Col: {1:s}'.format(event.getEventNumber(), col_name))
                col = event.getCollection(col_name)
                # Creating the CellID decocder
                cellIdEncoding = col.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                cellIdDecoder = UTIL.BitField64( cellIdEncoding )
                # print('  N elements: {0:d}'.format(col.getNumberOfElements()))
                histos = self.histos[col_name]
                # Filling the Tracker hit properties
                if col_type == 'SimTrackerHit':
                    for iHit in range(col.getNumberOfElements()):
                        hit = col.getElementAt(iHit)
                        cellId = long(hit.getCellID0() & 0xffffffff) | (long( hit.getCellID1() ) << 32)
                        cellIdDecoder.setValue(cellId)
                        side = int(cellIdDecoder['side'].value())
                        layer = int(cellIdDecoder['layer'].value()) + 1
                        # Skipping negative side disks
                        if side < 0:
                            continue
                        # Calculating the T0 based on the hit position in ns
                        t0 = hit.getPositionVec().Mag() / (CONST_C / 1e6)
                        hit_time = hit.getTime()
                        # Getting the MCParticle of the hit
                        mcp = hit.getMCParticle()
                        mcp_type = pdg_to_type(mcp.getPDG())
                        # Getting the oldest MCParticle of the hit
                        mcp_o = get_oldest_mcp_parent(mcp)
                        mcp_o_type = pdg_to_type(mcp_o.getPDG())
                        histos[('h_hit_mult_vs_mcp_type', layer)].Fill(mcp_o_type)
                        histos[('p_hit_e_vs_mcp_type', layer)].Fill(mcp_o_type, hit.getEDep()*1e6)
                        histos[('p_hit_time_vs_mcp_type', layer)].Fill(mcp_o_type, hit_time - t0)
                        if mcp_o_type in self.MCP_TYPES[col_type]:
                            histos[('h_hit_time_mcp_{0:d}'.format(mcp_o_type), layer)].Fill(hit_time - t0)
                # Filling the Calorimeter hit properties
                if col_type == 'SimCalorimeterHit':
                    for iHit in range(min(col.getNumberOfElements(), 10000)):
                        hit = col.getElementAt(iHit)
                        cellId = long(hit.getCellID0() & 0xffffffff) | (long( hit.getCellID1() ) << 32)
                        cellIdDecoder.setValue(cellId)
                        side = int(cellIdDecoder['side'].value())
                        layer = int(cellIdDecoder['layer'].value()) + 1
                        # Looping over the MCParticles of the hit
                        nMcp = hit.getNMCContributions()
                        t0 = hit.getPositionVec().Mag() / (CONST_C / 1e6)
                        br()
                        for iM in range(nMcp):
                            # Getting the MCParticle of the hit
                            # mcp = hit.getParticleCont(iM)
                            hit_time = hit.getTimeCont(iM)
                            hit_e = hit.getEnergyCont(iM)*1e6
                            # mcp_type = pdg_to_type(mcp.getPDG())
                            # Getting the oldest MCParticle of the hit
                            # mcp_o = get_oldest_mcp_parent(mcp)
                            # mcp_o_type = pdg_to_type(mcp_o.getPDG())
                            # histos['h_hit_mult_vs_mcp_type'].Fill(mcp_o_type)
                            # histos['p_hit_e_vs_mcp_type'].Fill(mcp_o_type, hit.getEnergyCont(iM)*1e6)
                            # histos['p_hit_time_vs_mcp_type'].Fill(mcp_o_type, hit_time - t0)
                            # histos['h2_hit_layer_vs_mcp_type'].Fill(mcp_o_type, layer)
                            histos['h2_hit_layer_vs_hit_time'].Fill(hit_time - t0, layer)
                            histos['h2_hit_layer_vs_hit_time_zoom'].Fill(hit_time - t0, layer)
                            histos['h2_hit_e_vs_hit_time'].Fill(hit_time - t0, hit_e)
                            histos['h_hit_time'].Fill(hit_time - t0)
                            histos['h_hit_e'].Fill(hit_e)
                            # if mcp_o_type in self.MCP_TYPES[col_type]:
                            #     histos['h_hit_time_mcp_{0:d}'.format(mcp_o_type)].Fill(hit_time - t0)
                            #     histos['h_hit_e_mcp_{0:d}'.format(mcp_o_type)].Fill(hit_e)


    def endOfData( self ):
        """Called by the event loop at the end of the loop"""
        
        # Storing histograms to the output ROOT file
        if self.output_path is not None:
            out_file = R.TFile(self.output_path, 'RECREATE')
            for col, histos in self.histos.iteritems():
                for hname, histo in histos.iteritems():
                    histo.Write()
            out_file.Close()
