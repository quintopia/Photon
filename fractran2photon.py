import sys
"""Fractran to Photon(Python) Transpiler"""

def transpile(listing):
    fracs,start = listing.splitlines()[:2]
    fracs = fracs.split()
    for i,frac in enumerate(fracs):
        fracs[i]=eval("("+frac.replace("/",",")+")")
    output = "(0,{},0,1,0,0)\n(".format(start)
    for i,frac in enumerate(fracs):
        output+="""
 (
  (-1,0,1) if y%{}==0 else
  (1,0,0)
 ) if x=={} else""".format(frac[1],i+1)
    output+="""
 (1,0,1) if x==0 else
 (0,0,0) if x=={} else
 (1,0,0)
) if z==0 else
(
 (1,0,-1) if x==0 else""".format(len(fracs)+1)
    for i,frac in enumerate(fracs):
        dir = (-1)**(frac[0]<frac[1])
        output += """
 (
  (0,{3},-1) if y%{0}==0 and y/{0}==z else
  (-1,{4},0) if y%{1}==0 and y/{1}==z else
  (0,{3},0)
 ) if x=={2} else""".format(frac[1],frac[0],i+1,dir,-dir)
    output+="""
 (1,0,0)
)"""
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