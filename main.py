
if __name__ == "__main__":
    import sys
    #cat = SphericalFusionCategory(sys.argv[1], sys.argv[2])
    ## tri = regina.Triangulation3.fromSig(sys.argv[3]) if len(sys.argv) > 3 \

    tri = regina.Example3.lens(5, 1)
    #tri = regina.Example3.poincare()
    cats = ['37fa4c58-97c6-11f1-a6a6-639d9c4bace6',
            '37fa4e2e-97c6-11f1-aa03-3b817be52f72',
            '37fa4ea8-97c6-11f1-8b25-ad9b74351a3e',
            '37fa4f78-97c6-11f1-90cf-c3d624fdb861',
            '37fa5126-97c6-11f1-8476-3f003633fef2',
            '37fa5642-97c6-11f1-9458-d142162199ed',
            '37fa5c3e-97c6-11f1-8693-91f67133af96']
    for cate in cats:
        cat = SphericalFusionCategory(cate, sys.argv[2])
        Z = turaev_viro(tri, cat)
        print(f'{sys.argv[1]} \t Rank:{cat.rank} \t Sum:{Z}')
    #print("Z (numeric)=", sp.N(Z))
