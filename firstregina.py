import regina

# A layered solid torus
T = regina.Example3.lst(3, 4)

print(T)
print("tetrahedra:", T.countTetrahedra())
print("vertices:", T.countVertices())
print("edges:", T.countEdges())
print("triangles:", T.countTriangles())
print("isoSig:", T.isoSig())

