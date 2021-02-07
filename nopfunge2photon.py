import sys
"""NopFunge to Photon(Python) Transpiler"""

def transpile(listing):
    lines = listing.splitlines()
    hsplit = [i for i,e in enumerate(lines) if e.startswith('=')][0]
    maxlen = len(lines[hsplit])
    vsplit = lines[0].index(';')
    lines = list(map(lambda s: s.ljust(maxlen),lines))
    hrepeat = maxlen-vsplit-1
    vrepeat = len(lines)-hsplit-1
    #v<^>
    selector = ["""(1,1,0) if (x{2}=={0} and (y{3}=={1} else
(-1,0,0) if (x{2}=={0}+1 and (y{3}=={1} else
""","""(-1,1,0) if (x{2}=={0} and (y{3}=={1} else
(0,-1,0) if (x{2}=={0} and (y{3}=={1}+1 else
""","""(1,-1,0) if (x{2}=={0} and (y{3}=={1} else
(-1,0,0) if (x{2}=={0}+1 and (y{3}=={1} else
""","""(1,1,0) if (x{2}=={0} and (y{3}=={1} else
(0,-1,0) if (x{2}=={0} and (y{3}=={1}+1 else
"""]
    modder="-{})%{}"
    xmod=modder.format(2*vsplit,2*hrepeat)
    ymod=modder.format(2*hsplit,2*vrepeat)
    output = """(0,1,0,1,0,0)
(0,0,0) if x*y==0 else
"""
    for i in range(maxlen):
        if i == vsplit:
            continue
        for j in range(len(lines)):
            if j==hsplit or lines[j][i] not in "v<^>":
                continue
            refl = selector["v<^>".index(lines[j][i])]
            if i<vsplit and j<hsplit:
                output+=refl.format(2*i+1,2*j+1,")",")")
            elif i>vsplit and j>hsplit:
                output+=refl.format(2*(i-vsplit),2*(j-hsplit),xmod,ymod)
            elif i>vsplit and j<hsplit:
                output+=refl.format(2*(i-vsplit),2*j+1,xmod,")")
            elif i<vsplit and j>hsplit:
                output+=refl.format(2*i+1,2*(j-hsplit),")",ymod)
    output+="(0,0,1)"
    return output



if __name__=="__main__":
    if len(sys.argv)<3:
        print("Usage: "+argv[0]+" inputprog.frac outputprog.pho")
    else:
        with open(sys.argv[1]) as f:
            listing = f.read()
        ptn = transpile(listing)
        with open(sys.argv[2],'w') as f:
            f.write(ptn)