def get_oldest_mcp_parent(mcp):
    """Recursively looks for the oldest parent of the MCParticle"""
    pars = mcp.getParents()
    if (len(pars) < 1):
        return mcp
    for par in pars:
        if par is mcp:
            continue
        return get_oldest_mcp_parent(par)

def pdg_to_type(pdgId):
    """Convert PDG ID of a particle to its type for histograms"""
    if pdgId == 2212:
        return 1
    if pdgId == 2112:
        return 2
    if pdgId == 211:
        return 3
    if pdgId == -211:
        return 4
    if pdgId == 321:
        return 5
    if pdgId == -321:
        return 6
    if pdgId == -13:
        return 7
    if pdgId == 13:
        return 8
    if pdgId == 22:
        return 9
    if pdgId == 11:
        return 10
    if pdgId == -11:
        return 11
    if pdgId == -2212:
        return 12
    if pdgId == 111:
        return 13
    if pdgId == 1000010020:
        return 14
    if pdgId == 1000010030:
        return 15
    if pdgId == 1000020030:
        return 16
    if pdgId == 1000020040:
        return 17
    if pdgId == 14:
        return 18
    if pdgId == -14:
        return 19
    if pdgId == 12:
        return 20
    if pdgId == -12:
        return 21
    if pdgId == 130:
        return 22
    if pdgId == 310:
        return 23
    if pdgId == 311:
        return 24
    if pdgId == -311:
        return 25
    if pdgId == 3122:
        return 26
    if pdgId == -3122:
        return 27
    if pdgId == 3222:
        return 28
    if pdgId == 3212:
        return 29
    if pdgId == 3112:
        return 30
    if pdgId == -2112:
        return 31
    if pdgId == 3322:
        return 32
    if pdgId == 3312:
        return 33
    if pdgId == 3334:
        return 34
    if pdgId == 5112:
        return 35
    if pdgId == 5212:
        return 36
    if pdgId == 5222:
        return 37
    if pdgId == -3322:
        return 38
    if pdgId == -5132:
        return 39
    if pdgId == -5332:
        return 40
    return 0
