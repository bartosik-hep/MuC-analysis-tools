import argparse

parser = argparse.ArgumentParser(description='Process hits from a file')
parser.add_argument('input', metavar='input.root', type=str, help='List of input files', nargs="+")
parser.add_argument('-m', '--max_events', metavar='N', type=int, help='Maximum number of events to process', default=-1)
parser.add_argument('-o', dest='output', metavar='OUT.root', type=str, help='Path to the output ROOT file')
parser.add_argument('-s', '--skip_events', metavar='N', type=int, help='Number of events to skip', default=0)

opts = parser.parse_args()

from pyLCIO.io.EventLoop import EventLoop
# from drivers.hits_timing import HitsTimingDriver as TheDriver
# from drivers.trk_props import TrkPropsDriver as TheDriver
# from drivers.trk_hit_props import TrkHitPropsDriver as TheDriver
# from drivers.hit_props import HitPropsDriver as TheDriver
# from drivers.pfo_props import PfoPropsDriver as TheDriver
# from drivers.vtx_hit_props import VtxHitPropsDriver as TheDriver
# from drivers.trk_efficiency import TrkEfficiencyDriver as TheDriver
# from drivers.trk_hit_density import HitDensityDriver as TheDriver
from drivers.trk_hits_mcp import TrkHitsMCPDriver as TheDriver
# from drivers.trk_hit_loopers import TrkHitLoopersDriver as TheDriver
# from drivers.cal_hits_mcp import CalHitsMCPDriver as TheDriver



print('### Starting analysis with {0:d} input files:'.format(len(opts.input)))

evLoop = EventLoop()
for infile in opts.input:
    print('  {0:s}'.format(infile))
    evLoop.addFile(infile)
# evLoop.addFile('/home/bartosik/clic/test3_py/v2.slcio')
# evLoop.addFile('/home/bartosik/clic/out/digi_bkg_QGSP_BERT_HP/c0_25ns_nEkin150MeV/sim_mod1_mumi-1e3x500-26m-lowth-excl_j1.slcio')
# evLoop.addFile('/home/bartosik/clic/out/digi_bkg_QGSP_BERT_HP/c0_25ns_nEkin150MeV/sim_mod1_mumi-1e3x500-26m-lowth-excl_j2.slcio')
# evLoop.addFile('/home/bartosik/clic/out/digi_bkg_QGSP_BERT_HP/c0_25ns_nEkin150MeV/sim_mod1_mumi-1e3x500-26m-lowth-excl_j3.slcio')
# evLoop.addFile('/home/bartosik/clic/out/digi_bkg_QGSP_BERT_HP/c0_25ns_nEkin150MeV/sim_mod1_mumi-1e3x500-26m-lowth-excl_j4.slcio')
# evLoop.addFile('/home/bartosik/clic/out/digi_bkg_QGSP_BERT_HP/c0_25ns_nEkin150MeV/sim_mod1_mumi-1e3x500-26m-lowth-excl_j5.slcio')
# evLoop.addFile('/home/bartosik/clic/out/digi_bkg_QGSP_BERT_HP/c0_25ns_nEkin150MeV/sim_mod1_mumi-1e3x500-26m-lowth-excl_j6.slcio')
# evLoop.addFile('/home/bartosik/clic/out/digi_bkg_QGSP_BERT_HP/c0_25ns_nEkin150MeV/sim_mod1_mumi-1e3x500-26m-lowth-excl_j7.slcio')
# evLoop.addFile('/home/bartosik/clic/out/digi_bkg_QGSP_BERT_HP/c0_25ns_nEkin150MeV/sim_mod1_mumi-1e3x500-26m-lowth-excl_j8.slcio')
# driver = HitsTimingDriver('/home/bartosik/clic/test3_py/plots/hits_timing_c0_25ns_nEkin150MeV.root')
print('### Will store output in: {0:s}'.format(opts.output))
nEvents = evLoop.reader.getNumberOfEvents()
print('### Total number of events in the files: {0:d}'.format(nEvents))
driver = TheDriver(opts.output)
evLoop.add(driver)

if opts.max_events > 0:
	nEvents = opts.max_events

print('### Starting the loop over {0:d} events'.format(nEvents))
# event = evLoop.reader.next()
if opts.skip_events:
	print('### Skipping {0:d} events'.format(opts.skip_events))
	evLoop.skipEvents(opts.skip_events)
evLoop.loop(nEvents)
evLoop.printStatistics()

print('### Finished')
