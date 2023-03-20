import ROOT as R
import numpy as np
import math
from pyLCIO.drivers.Driver import Driver
from pyLCIO import EVENT, UTIL

from pdb import set_trace as br

CONST_C = R.TMath.C()
T_MAX = 0.3 # ns
# T_MAX = 10.0 # ns
T_MIN = -1.0 # ns

def get_oldest_mcp_parent(mcp, nIters=0):
    """Recursively looks for the oldest parent of the MCParticle"""
    pars = mcp.getParents()
    if (len(pars) < 1):
        return mcp, nIters
    for par in pars:
        # Skipping if the particle is its own parent
        if par is mcp:
            continue
        return get_oldest_mcp_parent(par, nIters+1)

class CalHitsMCPDriver( Driver ):
    """Driver creating histograms of detector hits and their corresponding MCParticles"""

    HIT_COLLECTION_NAMES = ['ECalBarrelCollection', 'ECalEndcapCollection']
    # HIT_COLLECTION_NAMES = ['ECalBarrelCollection', 'ECalEndcapCollection',
    #                         'HCalBarrelCollection', 'HCalEndcapCollection']

    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.output_path = output_path


    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""

        names_F = ['edep', 'time', 'time0', 'path_len',
                   'pos_r', 'pos_z', 'pos_x', 'pos_y',
                   'mcp_vtx_r', 'mcp_vtx_z', 'mcp_vtx_x', 'mcp_vtx_y',
                   'mcp_bib_vtx_r', 'mcp_bib_vtx_z', 'mcp_bib_vtx_x', 'mcp_bib_vtx_y',
                   'mcp_theta', 'mcp_phi', 'mcp_bib_theta', 'mcp_bib_phi',
                   'mcp_time', 'mcp_bib_time',
                   'mcp_beta', 'mcp_gamma', 'mcp_e', 'mcp_p', 'mcp_pt', 'mcp_pz',
                   'mcp_bib_beta', 'mcp_bib_gamma', 'mcp_bib_e', 'mcp_bib_p', 'mcp_bib_pt', 'mcp_bib_pz'
                   ]
        names_I = ['layer', 'side', 'col_id',
                   'mcp_pdg', 'mcp_bib_pdg', 'mcp_bib_niters', 'mcp_gen', 'mcp_bib_gen']
        # Creating the TTree with branches
        self.data = {}
        self.tree = R.TTree('tree', 'SimTrackerHit properties')
        for name in names_F:
            self.data[name] = np.zeros(1, dtype=np.float32)
            self.tree.Branch(name, self.data[name], '{0:s}/F'.format(name))
        for name in names_I:
            self.data[name] = np.zeros(1, dtype=np.int32)
            self.tree.Branch(name, self.data[name], '{0:s}/I'.format(name))

    def processEvent( self, event ):
        """Called by the event loop for each event"""

        # Get the MCParticle collection from the event
        mcParticles = event.getMcParticles()

        # Loop over hits
        print('Event: {0:d}'.format(event.getEventNumber()))
        for iCol, col_name in enumerate(self.HIT_COLLECTION_NAMES):
            # print('Event: {0:d} Col: {1:s}'.format(event.getEventNumber(), col_name))
            col = event.getCollection(col_name)
            # print('  N elements: {0:d}'.format(col.getNumberOfElements()))
            # Creating the CellID decocder
            cellIdEncoding = col.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
            cellIdDecoder = UTIL.BitField64(cellIdEncoding)
            # Filling the Tracker hit properties
            data = self.data
            nHits = col.getNumberOfElements()
            # print('Checking {1:d} hits from: {0:s}'.format(col_name, nHits))
            for iHit in range(nHits):
                # if iHit % int(nHits/10) == 0:
                #     print('  hit {0:d} / {1:d}'.format(iHit, nHits))
                # Hit time information
                hit = col.getElementAt(iHit)
                # Decoding the CellID
                cellId = int(hit.getCellID0() & 0xffffffff) | (int( hit.getCellID1() ) << 32)
                cellIdDecoder.setValue(cellId)
                data['col_id'][0] = iCol
                data['side'][0] = int(cellIdDecoder['side'].value())
                data['layer'][0] = int(cellIdDecoder['layer'].value())
                # Hit general properties
                pos = hit.getPositionVec()
                data['pos_x'][0] = pos.X()
                data['pos_y'][0] = pos.Y()
                data['pos_z'][0] = pos.Z()
                data['pos_r'][0] = pos.Perp()
                t0 = pos.Mag() / (CONST_C / 1e6)
                data['time0'][0] = t0
                # Looping over hit contributions
                nC = hit.getNMCContributions()
                for iC in range(nC):
                    data['time'][0] = hit.getTimeCont(iC)
                    # Skipping hits outside of the time window
                    if (data['time'][0] - t0) > T_MAX:
                        continue
                    if (data['time'][0] - t0) < T_MIN:
                        continue
                    data['edep'][0] = hit.getEnergyCont(iC)
                    # MCParticle properties
                    mcp = hit.getParticleCont(iC)
                    mcp_bib, mcp_bib_niters = get_oldest_mcp_parent(mcp)
                    data['mcp_bib_niters'][0] = mcp_bib_niters
                    for prefix, part in {'mcp': mcp, 'mcp_bib': mcp_bib}.items():
                        pos = part.getVertex()
                        lv = part.getLorentzVec()
                        data[prefix+'_vtx_x'][0] = pos[0]
                        data[prefix+'_vtx_y'][0] = pos[1]
                        data[prefix+'_vtx_z'][0] = pos[2]
                        data[prefix+'_vtx_r'][0] = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1])
                        data[prefix+'_pdg'][0] = part.getPDG()
                        data[prefix+'_time'][0] = part.getTime()
                        data[prefix+'_gen'][0] = part.getGeneratorStatus()
                        data[prefix+'_theta'][0] = lv.Theta()
                        data[prefix+'_phi'][0] = lv.Phi()
                        data[prefix+'_p'][0] = lv.P()
                        data[prefix+'_pt'][0] = lv.Pt()
                        data[prefix+'_pz'][0] = lv.Pz()
                        data[prefix+'_beta'][0] = lv.Beta()
                        data[prefix+'_gamma'][0] = lv.Gamma()
                    self.tree.Fill()

        print('  Tree has {0:d} hits'.format(self.tree.GetEntries()))

    def endOfData( self ):
        """Called by the event loop at the end of the loop"""

        # Storing histograms to the output ROOT file
        if self.output_path is not None:
            out_file = R.TFile(self.output_path, 'RECREATE')
            self.tree.Write()
            out_file.Close()
