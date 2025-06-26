MAX_SECTOR = 18
SECTOR_Z = 1
SECTOR_CUSTOM = 18
worms = {'WormMeasGiw': 18,
         'WormMeasGtau': 18,
         'WormMeasGSigmaiw': 3,
         'WormMeasG4iw': 18,
         'WormMeasG4tau': 18,
         'WormMeasH4iw': 5,
         'WormMeasP2iwPH': 6,
         'WormMeasP2tauPH': 6,
         'WormMeasP2iwPP': 7,
         'WormMeasP2tauPP': 7,
         'WormMeasP3iwPH': 18,
         'WormMeasP3tauPH': 18,
         'WormMeasP3iwPP': 18,
         'WormMeasP3tauPP': 18,
         'WormMeasQQ': 10,
         'WormMeasQQtau': 10,
         'WormMeasQQQQ': 11,
         'WormMeasNQQdag': 12,
         'WormMeasQQdd': 13,
         'WormMeasUcaca': 14,
         'WormMeasUcacatau': 14,
         'WormMeasUccaa': 15,
         'WormMeasUccaatau': 15,
         'WormMeasQUDdag': 16,
         'WormMeasRaman': 17,
         'WormMeasCustomMat': 18,
         'WormMeasCustomTau': 18}

def get_sector_index(cfg):
    """
    Get the worm sector number corresponding to the parameters WormMeas...
    """
    sectors = []
    for worm in worms.keys():
        if cfg[worm] != 0:
            sectors.append(worms[worm])
    sectors = list(set(sectors))
    if len(sectors) > 1:
        raise ValueError('Multiple worm sectors {} given.'.format(sectors))
    elif len(sectors) < 1:
        raise ValueError('No worm sector given')
    else:
        sector = sectors[0]

    return sector
