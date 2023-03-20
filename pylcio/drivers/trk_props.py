import ROOT as R
import numpy as np
from pyLCIO.drivers.Driver import Driver

from pdb import set_trace as br

CONST_C = R.TMath.C()

class TrkPropsDriver( Driver ):
    """Driver creating histograms of detector hits timing and energy"""

    TRK_COLLECTIONS = ['SiTracks', 'SiTracksCT', 'SiTracks_Refitted']
    HIT_COLLECTIONS = ['VXDTrackerHits', 'VXDEndcapTrackerHits', 'ITrackerHits', 'ITrackerEndcapHits', 'OTrackerHits', 'OTrackerEndcapHits']
    HIT_REL_COLLECTIONS = ['VXDTrackerHitRelations', 'VXDEndcapTrackerHitRelations', 'InnerTrackerBarrelHitsRelationsLCRelation', 'InnerTrackerEndcapHitsRelationsLCRelation', 'OuterTrackerBarrelHitsRelationsLCRelation', 'OuterTrackerEndcapHitsRelationsLCRelation']
    
    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.histos = {}
        self.output_path = output_path
    
    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""

        for trk_type in self.TRK_COLLECTIONS:
            histos = {}
            name = 'trk_pt'
            histos[name] = R.TH1I('_'.join([trk_type, name]), ';Track p_{T} [GeV];Tracks', 500, 0, 100)
            name = 'trk_d0'
            histos[name] = R.TH1I('_'.join([trk_type, name]), ';Track D0 [mm];Tracks', 500, 0, 10)
            name = 'trk_z0'
            histos[name] = R.TH1I('_'.join([trk_type, name]), ';Track Z0 [mm];Tracks', 500, 0, 10)
            name = 'trk_ndf'
            histos[name] = R.TH1I('_'.join([trk_type, name]), ';Track Ndf;Tracks', 50, 0, 50)
            name = 'trk_chi2'
            histos[name] = R.TH1I('_'.join([trk_type, name]), ';Track Chi2;Tracks', 1000, 0, 100)
            name = 'trk_chi2_norm'
            histos[name] = R.TH1I('_'.join([trk_type, name]), ';Track Chi2/Ndf;Tracks', 500, 0, 5)
        
            self.histos[trk_type] = histos
    
    def processEvent( self, event ):
        """Called by the event loop for each event"""

        # Loop over hits
        print('Event: {0:d}'.format(event.getEventNumber()))

        for trk_type in self.TRK_COLLECTIONS:
            histos = self.histos[trk_type]
            trks = event.getCollection(trk_type)
            for iTrk in range(trks.getNumberOfElements()):
                trk = trks.getElementAt(iTrk)
                histos['trk_d0'].Fill(trk.getD0())
                histos['trk_z0'].Fill(trk.getZ0())
                histos['trk_ndf'].Fill(trk.getNdf())
                histos['trk_chi2'].Fill(trk.getChi2())
                histos['trk_chi2_norm'].Fill(trk.getChi2()/trk.getNdf())
                # lv = trk.get
                # br()

    def endOfData( self ):
        """Called by the event loop at the end of the loop"""
        
        # Storing histograms to the output ROOT file
        if self.output_path is not None:
            out_file = R.TFile(self.output_path, 'RECREATE')
            for trk_type, histos in self.histos.iteritems():
                for hname, histo in histos.iteritems():
                    histo.Write()
            out_file.Close()
