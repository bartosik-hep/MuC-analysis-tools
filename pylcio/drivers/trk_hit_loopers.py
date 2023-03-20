import ROOT as R
import numpy as np
import math
from copy import copy

from pyLCIO.drivers.Driver import Driver
from pyLCIO import EVENT, UTIL, IMPL, IO, IOIMPL

from pdb import set_trace as br

CONST_C = R.TMath.C()
# T_MAX = 0.18 # ns
T_MAX = 10e3 # ns
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

class TrkHitLoopersDriver( Driver ):
    """Driver creating histograms of detector hits associated to looper MCParticles"""

    # HIT_COLLECTION_NAMES = ['VertexBarrelCollection', 'VertexEndcapCollection']
    HIT_COLLECTION_NAMES = ['VertexBarrelCollection', 'VertexEndcapCollection',
                            'InnerTrackerBarrelCollection', 'InnerTrackerEndcapCollection',
                            'OuterTrackerBarrelCollection', 'OuterTrackerEndcapCollection']
    NHITS_MAX = 300
    HIT_F = ['hit_t',
             'hit_pos_z', 'hit_pos_x', 'hit_pos_y', 'hit_pos_r',
             ]
    HIT_I = ['hit_det']
    MCP_F = ['mcp_vtx_z', 'mcp_vtx_x', 'mcp_vtx_y', 'mcp_vtx_r',
             'mcp_theta', 'mcp_phi', 'mcp_t',
             'mcp_beta', 'mcp_gamma', 'mcp_e', 'mcp_p', 'mcp_pt', 'mcp_pz',
             ]
    MCP_I = ['mcp_pdg', 'mcp_nhits']


    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.output_path = output_path
        self.out_root = None
        self.out_lcio = None
        self.event = 0


    def clear_data(self):
        """Resets the data structures used for filling the TTree"""
        for name in self.MCP_F + self.MCP_I + self.HIT_I + self.HIT_F:
            self.data[name].fill(0)


    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""

        # Creating the TTree with branches
        self.data = {}
        self.tree = R.TTree('tree', 'SimTrackerHit properties')
        for name in self.MCP_F:
            self.data[name] = np.zeros(1, dtype=np.float32)
            self.tree.Branch(name, self.data[name], '{0:s}/F'.format(name))
        for name in self.MCP_I:
            self.data[name] = np.zeros(1, dtype=np.int32)
            self.tree.Branch(name, self.data[name], '{0:s}/I'.format(name))
        for name in self.HIT_I:
            self.data[name] = np.zeros(self.NHITS_MAX, dtype=np.int32)
            self.tree.Branch(name, self.data[name], '{0:s}[mcp_nhits]/I'.format(name))
        for name in self.HIT_F:
            self.data[name] = np.zeros(self.NHITS_MAX, dtype=np.float32)
            self.tree.Branch(name, self.data[name], '{0:s}[mcp_nhits]/F'.format(name))

        # Opening the output ROOT file
        if self.output_path is not None:
            self.out_root = R.TFile(self.output_path, 'RECREATE')

        # Opening the output LCIO file
        self.out_lcio = IOIMPL.LCFactory.getInstance().createLCWriter()
        self.out_lcio.open(self.output_path.replace('.root', '.slcio'), EVENT.LCIO.WRITE_NEW)


    def processEvent( self, event ):
        """Called by the event loop for each event"""

        eventNr = event.getEventNumber()
        print('Event: {0:d}'.format(eventNr))
        # Get the MCParticle collection from the event
        mcParticles = event.getMcParticles()
        # Loop over MCParticles
        nParticles = 0
        # Storing even as run in LCIO
        run = IMPL.LCRunHeaderImpl()
        run.setRunNumber(eventNr)
        self.out_lcio.writeRunHeader(run)
        for mcp in mcParticles:
            # Skipping non-electrons
            if abs(mcp.getPDG()) != 11:
                continue
            # Clearing the data structures
            self.clear_data()
            nHits_mcp = 0
            lcio_hits = {i: [] for i, _ in enumerate(self.HIT_COLLECTION_NAMES)}
            # Looking for matching hits in each collection
            for iCol, col_name in enumerate(self.HIT_COLLECTION_NAMES):
                col = event.getCollection(col_name)
                # Filling the Tracker hit properties
                data = self.data
                nHits = col.getNumberOfElements()

                # print('Checking {1:d} hits from: {0:s}'.format(col_name, nHits))
                for iHit in range(nHits):
                    # Checking if hit belong to the particle
                    hit = col.getElementAt(iHit)
                    mcp_h = hit.getMCParticle()
                    if mcp_h != mcp:
                        continue
                    # Storing hit for writing to LCIO
                    lcio_hits[iCol].append(hit)
                    # Storing hit properties
                    iH = nHits_mcp
                    pos = hit.getPositionVec()
                    data['hit_t'][iH] = hit.getTime()
                    # data['hit_t0'][iH] = pos.Mag() / (CONST_C / 1e6)
                    data['hit_det'][iH] = iCol
                    # data['hit_edep'][iH] = hit.getEDep()
                    data['hit_pos_x'][iH] = pos.X()
                    data['hit_pos_y'][iH] = pos.Y()
                    data['hit_pos_z'][iH] = pos.Z()
                    data['hit_pos_r'][iH] = pos.Perp()
                    nHits_mcp += 1

            # Storing info about MCParticle if it has 1+ hits
            if nHits_mcp == 0:
                continue

            pos = mcp.getVertex()
            lv = mcp.getLorentzVec()
            data['mcp_nhits'][0] = nHits_mcp
            data['mcp_vtx_x'][0] = pos[0]
            data['mcp_vtx_y'][0] = pos[1]
            data['mcp_vtx_z'][0] = pos[2]
            data['mcp_vtx_r'][0] = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1])
            data['mcp_pdg'][0] = mcp.getPDG()
            data['mcp_t'][0] = mcp.getTime()
            data['mcp_theta'][0] = lv.Theta()
            data['mcp_phi'][0] = lv.Phi()
            data['mcp_p'][0] = lv.P()
            data['mcp_pt'][0] = lv.Pt()
            data['mcp_pz'][0] = lv.Pz()
            data['mcp_beta'][0] = lv.Beta()
            data['mcp_gamma'][0] = lv.Gamma()
            nParticles += 1

            self.tree.Fill()

            # Adding MCParticle to the LCIO output
            print(f'Particle: {nParticles}')
            evt = IMPL.LCEventImpl()
            col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE)
            evt.setEventNumber(self.event)
            evt.setRunNumber(run.getRunNumber())
            evt.addCollection(col, "MCParticle")
            col.addElement(copy(mcp))
            col_hit = [IMPL.LCCollectionVec( EVENT.LCIO.SIMTRACKERHIT ) for col in self.HIT_COLLECTION_NAMES]
            # Adding hits to the LCIO output
            for idx, hits in lcio_hits.items():
                col = IMPL.LCCollectionVec(EVENT.LCIO.SIMTRACKERHIT)
                evt.addCollection(col, self.HIT_COLLECTION_NAMES[idx])
                for hit in hits:
                    col.addElement(copy(hit))
            self.out_lcio.writeEvent(evt)

            self.event += 1


        print(f'Saved  {nParticles}  particles')


    def endOfData( self ):
        """Called by the event loop at the end of the loop"""

        # Storing histograms to the output ROOT file
        if self.out_root:
            self.out_root.cd()
            self.tree.Write()
            self.out_root.Close()

        # Closing the LCIO file
        self.out_lcio.close()
