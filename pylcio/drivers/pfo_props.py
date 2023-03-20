import ROOT as R
import numpy as np
from pyLCIO.drivers.Driver import Driver

from pdb import set_trace as br

CONST_C = R.TMath.C()

class PfoPropsDriver( Driver ):
    """Driver creating histograms of detector hits and their corresponding MCParticles"""
    
    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.histos = {}
        self.output_path = output_path
    
    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""
        print('Start of data')

        # Book histograms
        name = 'h_pfo_energy'
        self.histos[name] = R.TH1I(name, ';Energy [GeV];PFOs', 3000, 0, 150)
    
    def processEvent( self, event ):
        """Called by the event loop for each event"""
        
        # Get the PFOs
        pfos = event.getCollection('PandoraPFOs')

        # Loop over hits
        evnum = event.getEventNumber()
        if evnum % 100 == 0:
            print('Event: {0:d}'.format(evnum))
        for iP in range(pfos.getNumberOfElements()):
            pfo = pfos.getElementAt(iP)
            self.histos['h_pfo_energy'].Fill(pfo.getEnergy())


    def endOfData( self ):
        """Called by the event loop at the end of the loop"""
        
        # Storing histograms to the output ROOT file
        if self.output_path is not None:
            out_file = R.TFile(self.output_path, 'RECREATE')
            for hname, histo in self.histos.iteritems():
                histo.Write(hname)
            out_file.Close()
