import ROOT as R
import numpy as np
import math
from pyLCIO import EVENT, UTIL
from pyLCIO.EVENT import TrackState
from pyLCIO.drivers.Driver import Driver

from pdb import set_trace as br

CONST_C = R.TMath.C()

class TrkEfficiencyDriver( Driver ):
    """Driver calculating track-reconstruction efficiencies"""

    # TRK_COLLECTION_NAME = 'SiTracksCT'
    TRK_COLLECTION_NAME = 'SiTracks_Refitted'
    SIMHIT_COLLECTION_NAMES = ['VertexBarrelCollection', 'VertexEndcapCollection', 'InnerTrackerBarrelCollection', 'InnerTrackerEndcapCollection', 'OuterTrackerBarrelCollection', 'OuterTrackerEndcapCollection']
    HIT_COLLECTION_NAMES = ['VXDBTrackerHits', 'VXDETrackerHits', 'IBTrackerHits', 'IETrackerHits', 'OBTrackerHits', 'OETrackerHits']
    HIT_RELATION_NAMES = ['VXDBTrackerHitsRelations', 'VXDETrackerHitsRelations', 'IBTrackerHitsRelations', 'IETrackerHitsRelations', 'OBTrackerHitsRelations', 'OETrackerHitsRelations']
    N_LAYERS = [8, 8, 3, 7, 3, 4]
    N_LAYERS_VTX = 16
    MCP_PDGS = [-13, 13]
    MAG_FIELD = 4.0  # Tesla
    # Cuts
    MCP_PT_MIN = 0.1
    N_HIT_LAYERS_MIN = 4
    N_VTX_LAYERS_MIN = 4



    def __init__( self, output_path=None):
        """Constructor"""
        Driver.__init__(self)
        self.histos = {}
        self.output_path = output_path
        self.N_LAYERS_TOTAL = sum(self.N_LAYERS)

    def startOfData( self ):
        """Called by the event loop at the beginning of the loop"""

        # Create histograms for MCParticles
        name = 'h_mcp_pt'
        self.histos[name] = R.TH1F( name, ';p_{T} [GeV];MCParticles', 2200,0,110)
        name = 'h_mcp_theta'
        self.histos[name] = R.TH1F( name, ';#Theta [rad];MCParticles', 320,0,3.2)
        name = 'h_mcp_abstheta'
        self.histos[name] = R.TH1F( name, ';|#Theta| [deg];MCParticles', 90,0,90)
        name = 'h_mcp_abstheta_reco'
        self.histos[name] = R.TH1F( name, ';|#Theta| [deg];MCParticles', 90,0,90)
        name = 'h_trk_pt'
        self.histos[name] = R.TH1F( name, ';p_{T} [GeV];Tracks', 2200,0,110)
        name = 'h_trk_theta'
        self.histos[name] = R.TH1F( name, ';#Theta [rad];Tracks', 320,0,3.2)
        name = 'h_trk_dR'
        self.histos[name] = R.TH1F( name, ';dR(trk, MCParticle);Tracks', 500,0,0.1)
        name = 'h_trk_dpt'
        self.histos[name] = R.TH1F( name, ';dPt [%];Tracks', 200,-50,50)
        name = 'p_trk_theta_pt'
        self.histos[name] = R.TProfile( name, ';|#Theta| [deg];p_{T} [GeV]', 9,0,90, 's')
        name = 'h_ntrk'
        self.histos[name] = R.TH1F( name, ';# tracks;Events', 10, 0, 10)
        name = 'h_nhits'
        self.histos[name] = R.TH1F( name, ';# hits;Tracks', 30, 0, 30)


    def trk_vec(self, trk):
        """Returns a TVector3 for the track"""
        # Getting the track state at IP
        ts_ip = trk.getTrackState(TrackState.AtIP)
        # Calculating vector components
        trk_pt = 0.0003 * self.MAG_FIELD / abs(ts_ip.getOmega())
        trk_phi = ts_ip.getPhi()
        trk_theta = R.TMath.PiOver2() - math.atan(ts_ip.getTanLambda())
        v_t = R.TVector3()
        v_t.SetPtThetaPhi(trk_pt, trk_theta, trk_phi)
        return v_t

    def abstheta(self, theta):
        """Converts theta to the absolute value in degrees"""
        return abs(theta - R.TMath.PiOver2()) * R.TMath.RadToDeg()


    def processEvent( self, event ):
        """Called by the event loop for each event"""

        # if (event.getRunNumber(), event.getEventNumber()) != (1, 7678):
        #     return

        # Get the MCParticle collection from the event
        mcParticles = event.getMcParticles()
        trks = event.getCollection(self.TRK_COLLECTION_NAME)
        H = self.histos

        simhitcols = [event.getCollection(col) for col in self.SIMHIT_COLLECTION_NAMES]
        hitcols = [event.getCollection(col) for col in self.HIT_COLLECTION_NAMES]
        hitrels = [event.getCollection(col) for col in self.HIT_RELATION_NAMES]

        # Loop over all gen-level MCParticles
        nmcp = 0
        for mcp in mcParticles:
            if mcp.getGeneratorStatus() != 1:
                continue
            pdg = mcp.getPDG()
            if pdg not in self.MCP_PDGS:
                continue
            # Get the four vector of the MCParticle
            lv_m = mcp.getLorentzVec()
            v_m = lv_m.Vect()
            if lv_m.Pt() < self.MCP_PT_MIN:
                continue
            # # Checking the MCParticle Theta
            # if not math.radians(20) < abs(v_m.Theta()) < math.radians(25):
            #     continue
            # print(event.getRunNumber(), event.getEventNumber())
            # continue
            nmcp += 1
            # Getting the SimHits belonging to this particle
            simhits = []
            rechits = []
            nhits_sim = np.zeros(self.N_LAYERS_TOTAL, dtype=np.uint8)
            nhits_rec = np.zeros(self.N_LAYERS_TOTAL, dtype=np.uint8)
            for iCol, simhitcol in enumerate(simhitcols):
                # Creating the CellID decocder
                cellIdEncoding = simhitcol.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                cellIdDecoder = UTIL.BitField64(cellIdEncoding)
                layer_offset = sum(self.N_LAYERS[:iCol])
                nSimHits = simhitcol.getNumberOfElements()
                for iSimHit in range(nSimHits):
                    simHit = simhitcol.getElementAt(iSimHit)
                    if simHit.getMCParticle() != mcp:
                        continue
                    simhits.append(simHit)
                    # Getting the layer ID of the hit
                    cellId = int(simHit.getCellID0() & 0xffffffff) | (int( simHit.getCellID1() ) << 32)
                    cellIdDecoder.setValue(cellId)
                    layer = int(cellIdDecoder['layer'].value())
                    nhits_sim[layer_offset+layer] += 1
                    # Finding the corresponding RecHit
                    rels = hitrels[iCol]
                    nRels = rels.getNumberOfElements()
                    for iRel in range(nRels):
                        if rels.getElementAt(iRel).getTo() != simHit:
                            continue
                        recHit = rels.getElementAt(iRel).getFrom()
                        rechits.append(recHit)
                        nhits_rec[layer_offset+layer] += 1
            nlayers_sim = len(nhits_sim[nhits_sim>0])
            nlayers_rec = len(nhits_rec[nhits_rec>0])
            nhits_vtx = nhits_rec[:self.N_LAYERS_VTX]
            nlayers_vtx = len(nhits_vtx[nhits_vtx>0])
            if nlayers_sim < self.N_HIT_LAYERS_MIN:
                continue
            # print("{4:d} {5:d}: # sim hits (layers): {0:d} ({1:d})\treco hits (layers): {2:d} ({3:d})".format(np.sum(nhits_sim), nlayers_sim,
            #                                                                                                   np.sum(nhits_rec), nlayers_rec,
            #                                                                                                   event.getEventNumber(),
            #                                                                                                   event.getRunNumber()))
            # br()

            H['h_mcp_pt'].Fill(lv_m.Pt())
            H['h_mcp_theta'].Fill(lv_m.Theta())
            theta_m = 90 - abs(lv_m.Theta() * R.TMath.RadToDeg() - 90)
            H['h_mcp_abstheta'].Fill( theta_m )
            # print('MCP: pt: {0:2f}  theta: {1:.2f}  phi: {2:.2f}'.format(lv_m.Pt(), lv_m.Theta(), lv_m.Phi()))

            # Finding the closest track
            ntrk = trks.getNumberOfElements()
            H['h_ntrk'].Fill(ntrk)
            # print('  # tracks: {0:d}'.format(ntrk))
            dR_min = 999
            trk_id = -1
            for iTrk in range(ntrk):
                trk = trks.getElementAt(iTrk)
                v_t = self.trk_vec(trk)
                dR = v_t.DeltaR(v_m)
                if dR < dR_min:
                    dR_min = dR
                    trk_id = iTrk
                # print('trk {0:d}: pt: {1:.2f} dR: {2:.2f}'.format(iTrk, v_t.Pt(), dR))
            if trk_id == -1:
                # Checking the MCParticle Theta
                if math.radians(22) < abs(v_m.Theta()) < math.radians(27):
                    print(event.getEventNumber(), event.getRunNumber())
                continue
            # Filling histograms with track properties
            trk = trks.getElementAt(iTrk)
            v_t = self.trk_vec(trk)
            nhits = len(trk.getTrackerHits())
            # print('  trk {0:d}: pt: {1:2f}  theta: {2:.2f}  phi: {3:.2f}  dR: {4:.2f}'.format(iTrk, v_t.Pt(), v_t.Theta(), v_t.Phi(), dR_min))
            H['h_trk_dR'].Fill(dR_min)
            if dR_min > 0.01:
                continue
            dpt = (v_t.Pt() - v_m.Pt()) / v_m.Pt()
            if abs(dpt) > 1.0:
                continue

            H['h_mcp_abstheta_reco'].Fill(theta_m)
            H['p_trk_theta_pt'].Fill(theta_m, v_t.Pt())
            H['h_trk_pt'].Fill(v_t.Pt())
            H['h_trk_theta'].Fill(v_t.Theta())
            H['h_trk_dpt'].Fill(100.0*dpt)
            H['h_nhits'].Fill(nhits)

        if nmcp != 1:
            print('#######  MCParticles: {0:d}'.format(nmcp))



    def endOfData( self ):
        """Called by the event loop at the end of the loop"""

        # Storing histograms to the output ROOT file
        if self.output_path is not None:
            out_file = R.TFile(self.output_path, 'RECREATE')
            for hname, histo in self.histos.items():
                histo.Write(hname)
            out_file.Close()
