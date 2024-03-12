import numpy as np
from numba import jit
import awkward as ak

class loop_matcher():
     
    def __init__(self, ids, dr = None, dpt = None) -> None:
        self.matchID = ids
        self.maxDR = dr
        self.maxPt = dpt
        self.checkCharge = False   
    # @jit
    def process(self, src: ak.Array, matched: ak.Array, builder):
        
        nev = len(src)
        builder.begin_record()
        builder.field("DgenPartIdx")
        for iev in range(nev):
            builder.begin_list()
            isrc = src[iev]
            imatch = matched[iev]
            nsrc = len(isrc)
            nmatch = len(imatch)
            for i0 in range(nsrc):
                # tmpids = []            
                with builder.list():
                    for i1 in range(nmatch):
                        mat_id = imatch.pdgId[i1] in self.matchID
                        mat_dr = isrc[i0].delta_r(imatch[i1]) < self.maxDR
                        mat_dpt = abs(isrc.pt[i0] - imatch.pt[i1])/isrc.pt[i0] < self.maxPt
                        if mat_id and mat_dr and mat_dpt:
                            # print(i1, type(i1))
                            builder.integer(np.int16(i1))
                # src.begin_tuple(len(tmpids))
                # for i in range(len(tmpids)):
                #     src.index(i).integer(tmpids[i])
                # src.end_tuple()
      
            builder.end_list()
            # print(iev)
        builder.end_record()
        return builder       