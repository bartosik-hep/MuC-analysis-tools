import ROOT as R
import numpy as np
import math
from pyLCIO.drivers.Driver import Driver
from pyLCIO import EVENT, UTIL

from pdb import set_trace as br

CONST_C = R.TMath.C()

class TrkHitPropsDriver( Driver ):
    """Driver creating histograms of tracker hits timing and energy"""

    HIT_COLLECTIONS = ['VXDTrackerHits', 'VXDEndcapTrackerHits', 'ITrackerHits', 'ITrackerEndcapHits', 'OTrackerHits', 'OTrackerEndcapHits']
    # SIMHIT_COLLECTIONS = ['VertexBarrelCollection', 'VertexEndcapCollection', 'InnerTrackerBarrelCollection', 'InnerTrackerEndcapCollection', 'OuterTrackerBarrelCollection', 'OuterTrackerEndcapCollection']
    SIMHIT_COLLECTIONS = ['VertexBarrelCollection']
    HIT_REL_COLLECTIONS = ['VXDTrackerHitRelations', 'VXDEndcapTrackerHitRelations', 'InnerTrackerBarrelHitsRelations', 'InnerTrackerEndcapHitsRelations', 'OuterTrackerBarrelHitsRelations', 'OuterTrackerEndcapHitsRelations']

    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.histos = {}
        self.output_path = output_path

    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""

        for trk_type in self.HIT_COLLECTIONS + self.SIMHIT_COLLECTIONS:
            histos = {}
            name = 'e'
            histos[name] = R.TH1I('_'.join([trk_type, name]), ';Hit energy [KeV];Hits', 300, 0, 150)
            name = 't_mt0_e'
            histos[name] = R.TH2I('_'.join([trk_type, name]), ';Hit time - T0 [ns];Hit energy [KeV]', 1100,-1,10, 300,0,150)
            name = 't'
            histos[name] = R.TH1I('_'.join([trk_type, name]), ';Hit time [ns];Hits', 1100, -1, 10)
            name = 't_mt0'
            histos[name] = R.TH1I('_'.join([trk_type, name]), ';Hit time - T0 [ns];Hits', 1100, -1, 10)
            name = 'pos_zy'
            histos[name] = R.TH2I('_'.join([trk_type, name]), ';Hit Z [mm];Hit Y [mm]', 4000, -2000, 2000, 4000, -2000, 2000)
            name = 'pos_xy'
            histos[name] = R.TH2I('_'.join([trk_type, name]), ';Hit X [mm];Hit Y [mm]', 4000, -2000, 2000, 4000, -2000, 2000)

            self.histos[trk_type] = histos
        self.eTot = np.zeros(8, dtype=np.float32)

    def processEvent( self, event ):
        """Called by the event loop for each event"""

        # Loop over hits
        print('Event: {0:d}'.format(event.getEventNumber()))

        # for iT, trk_type in enumerate(self.HIT_COLLECTIONS + self.SIMHIT_COLLECTIONS):
        for iT, trk_type in enumerate(self.SIMHIT_COLLECTIONS):
            histos = self.histos[trk_type]
            hits = event.getCollection(trk_type)
            # hits_rel = event.getCollection(self.HIT_REL_COLLECTIONS[iT])
            # Creating the CellID decocder
            cellIdEncoding = hits.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
            cellIdDecoder = UTIL.BitField64(cellIdEncoding)
            for iHit in range(hits.getNumberOfElements()):
                hit = hits.getElementAt(iHit)
                # Checking layer
                cellId = int(hit.getCellID0() & 0xffffffff) | (int( hit.getCellID1() ) << 32)
                cellIdDecoder.setValue(cellId)
                layer = int(cellIdDecoder['layer'].value())
                # Checking properties
                self.eTot[layer] += hit.getEDep()
                continue
                histos['e'].Fill(hit.getEDep()*1e6)
                histos['t'].Fill(hit.getTime())
                # Calculating T0
                pos_vec = hit.getPositionVec()
                t0 = pos_vec.Mag() / (CONST_C / 1e6)
                histos['t_mt0'].Fill(hit.getTime() - t0)
                histos['t_mt0_e'].Fill(hit.getTime() - t0, hit.getEDep()*1e6)
                histos['pos_zy'].Fill(pos_vec.Z(), pos_vec.Y())
                histos['pos_xy'].Fill(pos_vec.X(), pos_vec.Y())

    def endOfData( self ):
        """Called by the event loop at the end of the loop"""

        print('Total deposited energy in GeV')
        print(self.eTot)
        return

        # Storing histograms to the output ROOT file
        if self.output_path is not None:
            out_file = R.TFile(self.output_path, 'RECREATE')
            for trk_type, histos in self.histos.iteritems():
                for hname, histo in histos.iteritems():
                    histo.Write()
            out_file.Close()
