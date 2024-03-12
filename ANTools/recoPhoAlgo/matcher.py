import awkward as ak



class matcher():
    
    def __init__(self, ids, dr = None, dpt = None) -> None:
        self.matchID = ids
        self.maxDR = dr
        self.maxPt = dpt
        self.checkCharge = False
        
    def process(self, src: ak.Array, matched: ak.Array):
        
            
        combination = ak.argcartesian([src, matched], axis=1)
        
        conID = ak.where(src[combination["0"]].pdgId == 0, False, False) # initialize
        for idnum in self.matchID:
            conID = conID | (abs(src[combination["0"]].pdgId) == idnum )
        # require the angulardistance
        if self.maxDR is not None:
            # print("dr is working")
            condr = src[combination["0"]].delta_r(matched[combination["1"]]) < self.maxDR
        # require the momentum difference
        if self.maxPt is not None:
            conpt = (abs(src[combination["0"]].pt - matched[combination["1"]].pt) / src[combination["0"]].pt) < self.maxPt
            

        combination = combination[condr & conpt & conID]
        # idxa, idxb = ak.unzip(combination)
        match_table = []
        # remove double match
        for i in range(ak.max(ak.num(src, axis=1), axis=0)):
            idx_mask = combination["0"] == i
            # find the minimum angular distance 
            sub_combination = combination[idx_mask]
            dr_list = src[sub_combination["0"]].delta_r(matched[sub_combination["1"]])
            if ak.max(dr_list, axis=None) == None:
                match_table.append([None]*len(src))
                continue
            # print(ak.max(dr_list, axis=None))
            ismin = ak.min(dr_list, axis=1) == dr_list    
            match_table.append(ak.min(sub_combination["1"][ismin], axis=1))
           
        builder = ak.ArrayBuilder()
        for iev in range(len(src)):
            builder.begin_list()
            nsrc = len(src[iev])
            for isrc in range(nsrc):
                matidx = match_table[isrc][iev]
                matidx = -999 if matidx == None else matidx
                builder.integer(matidx)
            builder.end_list()
        results = builder.snapshot()
        
        return results