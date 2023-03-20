import ROOT as R
import numpy as np
import math
from pyLCIO.drivers.Driver import Driver
from pyLCIO import EVENT, UTIL

from pdb import set_trace as br


class HitDensityDriver( Driver ):
    """Driver calculating per-layer maximum hit density"""

    # HIT_COLLECTION_NAMES = [ 'VertexBarrelCollection', 'VertexEndcapCollection',
    #                          'InnerTrackerBarrelCollection', 'InnerTrackerEndcapCollection',
    #                          'OuterTrackerBarrelCollection', 'OuterTrackerEndcapCollection']
    # HIT_COLLECTION_NAMES = ['VXDBTrackerHits', 'VXDETrackerHits',
    #                         'IBTrackerHits', 'IETrackerHits',
    #                         'OBTrackerHits', 'OETrackerHits']
    HIT_COLLECTION_NAMES = ['VXDBTrackerHits_DL2', 'VXDETrackerHits_DL2']
    # N_LAYERS = [8, 8, 3, 7, 3, 4]
    N_LAYERS = [8, 8]
    MODULE_MAX = 10000


    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.histos = {}
        self.output_path = output_path

    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""
        for name in ['min', 'max', 'mean', 'median', 'sum']:
            hname = 'h_nhits_{0:s}'.format(name)
            self.histos['h_nhits_'+name] = R.TH1F(hname, ';Layer;N hits', sum(self.N_LAYERS), 0, sum(self.N_LAYERS))


    def processEvent( self, event ):
        """Called by the event loop for each event"""

        layer_id_offset = 0
        # Loop over each hit collection
        for iCol, colName in enumerate(self.HIT_COLLECTION_NAMES):
            # Shifting the histogram layer ID for this collection
            if iCol > 0:
                layer_id_offset += self.N_LAYERS[iCol-1]
            hits = event.getCollection(colName)
            nhits = hits.getNumberOfElements()
            print('### Processing collection `{0:s}` with {1:d} hits'.format(colName, nhits))
            H = self.histos
            # Creating the arrays of sensor hit counts for each layer
            hit_counts = [{} for i in range(self.N_LAYERS[iCol])]
            # Creating the CellID decocder
            cellIdEncoding = hits.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
            cellIdDecoder = UTIL.BitField64(cellIdEncoding)
            # Loop over hits
            for iH in range(nhits):
                hit = hits.getElementAt(iH)
                cellId = int(hit.getCellID0() & 0xffffffff) | (int( hit.getCellID1() ) << 32)
                cellIdDecoder.setValue(cellId)
                layer = int(cellIdDecoder['layer'].value())
                if cellId not in hit_counts[layer]:
                    hit_counts[layer][cellId] = 1
                else:
                    hit_counts[layer][cellId] += 1
            # Converting dicts to arrays
            for i in range(len(hit_counts)):
                hit_counts[i] = np.array(list(hit_counts[i].values()), dtype=np.uint32)
            # Filling histograms
            for layer in range(self.N_LAYERS[iCol]):
                print('  {0:d} sensors hit in layer {1:d}'.format(len(hit_counts[layer]), layer))
                self.histos['h_nhits_min'].SetBinContent(1+layer_id_offset + layer, np.min(hit_counts[layer]))
                self.histos['h_nhits_max'].SetBinContent(1+layer_id_offset + layer, np.max(hit_counts[layer]))
                self.histos['h_nhits_mean'].SetBinContent(1+layer_id_offset + layer, np.mean(hit_counts[layer]))
                self.histos['h_nhits_median'].SetBinContent(1+layer_id_offset + layer, np.median(hit_counts[layer]))
                self.histos['h_nhits_sum'].SetBinContent(1+layer_id_offset + layer, np.sum(hit_counts[layer]))




    def endOfData( self ):
        """Called by the event loop at the end of the loop"""

        # Storing histograms to the output ROOT file
        if self.output_path is not None:
            out_file = R.TFile(self.output_path, 'RECREATE')
            for hname, histo in self.histos.items():
                histo.Write()
            out_file.Close()
