import ROOT as R
import copy

def draw(histos, config, out_file=None):
    C = R.TCanvas("canvas", "", config['canvas'][0], config['canvas'][1])
    C.SetRightMargin(0.05)
    C.SetGridx()
    C.SetGridy()
    R.SetOwnership(C, False)
    n = len(histos)
    leg = None
    if 'leg' in config:
        leg = R.TLegend(0.75,0.9-n*0.06, 0.95,0.9)
        R.SetOwnership(leg, 0)
    h_axis = None
    if 'h_axis' in config:
        h_axis = config['h_axis']
        h_axis.Draw('AXIS')
        h_axis.Draw('AXIG')
    for iH, h_ori in enumerate(histos):
        drawopt = 'hist'
        if 'drawopt' in config:
            drawopt = config['drawopt']
        if iH > 0 or 'h_axis' in config:
            drawopt += "same"
        h = h_ori.Clone()
        if iH == 0:
            h_axis = h
        if 'rebin' in config:
            h.Rebin(config['rebin'])
#         h.Rebin(2)
        if 'norm' in config:
            h.Scale(config['norm'] / h.GetEntries())
        R.SetOwnership(h, 0)
        h.SetLineWidth(2)
        if 'style' in config:
            h.SetLineStyle(config['style'][iH])
        h.SetLineColor(config['color'][iH])
        h.SetMarkerColor(config['color'][iH])
        if leg:
            leg.AddEntry(h, config['leg'][iH], 'l')
        h.Draw(drawopt)
    if leg:
        leg.Draw()
    if 'y' in config:
        h_axis.GetYaxis().SetRangeUser(*config['y'])
    if 'x' in config:
        h_axis.GetXaxis().SetRangeUser(*config['x'])
    h_axis.GetYaxis().SetMaxDigits(4)
    if 'logY' in config:
        C.SetLogy(config['logY'])
    if 'logX' in config:
        C.SetLogx(config['logX'])
    C.RedrawAxis()
    C.Draw()
    if out_file:
        C.Print(out_file)


def read_root_obj(file_path, obj_path, path_delimiter='/', file_delimiter=':'):
    """Finds an object inside a file supporting ROOT file inputs"""

    if isinstance(file_path, str):
        # Stopping if input file doesn't exist
        if not os.path.isfile(file_path):
            return None
        # Opening the ROOT file
        file_in = R.TFile(file_path)
    else:
        # Using the already opened ROOT file
        file_in = file_path

    obj = None

    path_sequence = obj_path.split(path_delimiter)
    # Trying every element of the sequence
    for i_el, path_el in enumerate(path_sequence):
        # Getting the object from current directory when it's the last element
        if i_el == len(path_sequence) - 1:
            if path_el.startswith('*'):
                obj = R.gDirectory.FindObjectAny(path_el)
            else:
                obj = R.gDirectory.Get(path_el)
            break
        # Entering directly into the directory if full name is provided
        if '*' not in path_el:
            R.gDirectory.Get(path_el).cd()
        else:
            dir_list = R.gDirectory.GetListOfKeys()
            for key in dir_list:
                dir_name = key.GetName()
                if not re.match(path_el, dir_name):
                    continue
                R.gDirectory.Get(dir_name).cd()

    # Copying the object into memory and closing the file
    obj = copy.deepcopy(obj)
    file_in.Close()

    return obj
