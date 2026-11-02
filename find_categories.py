import json
from TVBW import SphericalFusionCategory

if __name__ == "__main__":
    import sys

    json_path = "Rank1-7.json"

    #cat = SphericalFusionCategory(sys.argv[1], sys.argv[2])

    cats = ["0" for _ in range(7)]
    cats[0] = '37fa4c58-97c6-11f1-a6a6-639d9c4bace6'
    with open(json_path) as fh:
        data = json.load(fh)
    order = data["order"]
    for c in order:
        #cat = SphericalFusionCategory(c, json_path)
        cat = data["data"][c]
        r = len(cat["quantum_dimensions"])
        if False not in [d in [[-1.0, 0.0], [1.0, 0.0]] for d in cat["numeric_quantum_dimensions"]] and cats[r - 1] == "0":
            print("found all ones", r)
            cats[r-1] = c
            if "0" not in cats:
                break
    print(cats)
    '''
    for cate in cats:
        cat = SphericalFusionCategory(cate, sys.argv[2])
        Z = turaev_viro(tri, cat)
        print(f'{sys.argv[1]} \t Rank:{cat.rank} \t Sum:{Z}')
    #print("Z (numeric)=", sp.N(Z))
    '''
