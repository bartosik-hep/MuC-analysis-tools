import ROOT as R
import numpy as np
import math
from pyLCIO import EVENT, UTIL
from pyLCIO.drivers.Driver import Driver

from pdb import set_trace as br

CONST_C = R.TMath.C()

class VtxHitPropsDriver( Driver ):
    """Driver creating histograms of Vertex hits timing and double-layer properties"""

    # HIT_COLLECTIONS = ['VXDTrackerHits', 'VXDEndcapTrackerHits']
    # HIT_COLLECTIONS = ['VertexBarrelCollection', 'VertexEndcapCollection']
    HIT_COLLECTIONS = ['VertexBarrelCollection']

    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.histos = {}
        self.output_path = output_path

    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""

        for trk_type in self.HIT_COLLECTIONS:
            histos = {}
            name = 'e'
            histos[name] = R.TH1F('_'.join([trk_type, name]), ';Hit energy [KeV];Hits', 300, 0, 150)
            name = 't'
            histos[name] = R.TH1F('_'.join([trk_type, name]), ';Hit time [ps];Hits', 700, -200, 500)
            name = 'layer'
            histos[name] = R.TH1F('_'.join([trk_type, name]), ';Layer ID;Hits', 10, 0, 10)
            name = 'side'
            histos[name] = R.TH1F('_'.join([trk_type, name]), ';Side ID;Hits', 10, -5, 5)
            name = 'module'
            histos[name] = R.TH1F('_'.join([trk_type, name]), ';Module ID;Hits', 50, 0, 50)
            name = 'sensor'
            histos[name] = R.TH1F('_'.join([trk_type, name]), ';Sensor ID;Hits', 10, 0, 10)

            self.histos[trk_type] = histos

    def processEvent( self, event ):
        """Called by the event loop for each event"""

        # Loop over hits
        print('Event: {0:d}'.format(event.getEventNumber()))
        for iT, trk_type in enumerate(self.HIT_COLLECTIONS):
            histos = self.histos[trk_type]
            hits = event.getCollection(trk_type)
            # Creating the CellID decocder
            cellIdEncoding = hits.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
            cellIdDecoder = UTIL.BitField64( cellIdEncoding )
            for iHit in range(hits.getNumberOfElements()):
                hit = hits.getElementAt(iHit)
                histos['e'].Fill(hit.getEDep()*1e6)
                histos['t'].Fill(hit.getTime()*1e3)
                # Filling hit position
                cellId = int(hit.getCellID0() & 0xffffffff) | (int( hit.getCellID1() ) << 32)
                cellIdDecoder.setValue(cellId)
                side = int(cellIdDecoder['side'].value())
                layer = int(cellIdDecoder['layer'].value())
                module = int(cellIdDecoder['module'].value())
                sensor = int(cellIdDecoder['sensor'].value())
                histos['layer'].Fill(layer)
                histos['side'].Fill(side)
                histos['module'].Fill(module)
                histos['sensor'].Fill(sensor)
                br()


    def endOfData( self ):
        """Called by the event loop at the end of the loop"""

        # Storing histograms to the output ROOT file
        if self.output_path is not None:
            out_file = R.TFile(self.output_path, 'RECREATE')
            for trk_type, histos in self.histos.items():
                for hname, histo in histos.items():
                    histo.Write()
            out_file.Close()
