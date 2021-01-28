import sys
import numpy as np
import math
import random

#utility function to make programs easier to write
def sgn(x):
    return x and x/abs(x)

class Photon:

    class Halt(RuntimeError):
        pass
        
    @staticmethod
    def normalized(a):
        if any(a!=0):
            return a/np.linalg.norm(a)
        return a
        
    def doprint(self):
        if self.verbose:
            print(self.printbuffer)
        self.printbuffer = ""
            
    def vprint(self,text):
        self.printbuffer += text+"\n"
        
    def maybeclearbuffer(self):
        if not self.showmiss:
            self.printbuffer = ""
    
    def __init__(self,listing,filename,verbose=False,showmiss=False):
        self.verbose = verbose
        self.showmiss = showmiss
        self.printbuffer = ""
        parts = listing.splitlines()
        if "[verbose]" in parts[0]:
            if "[showmiss]" in parts[0]:
                self.showmiss = True
            self.verbose = True
            parts.pop(0)
        x,y,z,vx,vy,vz = eval(parts[0])
        self.start=np.array([x,y,z])
        self.grad=np.array([vx,vy,vz])
        self.center=np.array([round(x) for x in self.start])
        expr="lambda x,y,z:"+"".join(parts[1:])
        if "for" in expr:
            raise SyntaxError("Generator expressions are BANNED",(filename,2,expr.index("for"),expr+"\n"))
        self.fx = eval(expr)
        self.t = 0
    
    def curpoint(self):
        return self.start+self.t*self.grad
    
    def update_t(self):
        next_t=np.inf
#        if round(self.curpoint()[1]+4.5)==0:
#            import pdb;pdb.set_trace()
        for i,dim in enumerate(self.curpoint()):
            if self.grad[i]==0:
                continue
            if abs(dim-0.5-round(dim-0.5))<0.000000000001:
                nextboundary=dim+np.sign(self.grad)[i]
            elif self.grad[i]<0:
                nextboundary=math.floor(dim+0.5)-0.5
            else:
                nextboundary=math.ceil(dim+0.5)-0.5
            candidate_t = np.roots([self.grad[i],self.start[i]-nextboundary])[0]
            if next_t>candidate_t>self.t:
                next_t=candidate_t
        self.t=next_t
        
    def identify_cube(self):
        candidates=[[],[],[]]
        for i,dim in enumerate(self.curpoint()):
            if abs(dim-0.5-round(dim-0.5))<0.00000000001:
                candidates[i]+=[round(dim-0.5),round(dim+0.5)]
            else:
                candidates[i]+=[round(dim)]
        center=[]
        dist=np.inf
        for x in candidates[0]:
            for y in candidates[1]:
                for z in candidates[2]:
                    point=np.array([x,y,z])
                    candist=np.linalg.norm(self.curpoint()-point)
                    if candist<dist or (candist==dist and np.dot(self.grad,point-self.curpoint())>np.dot(self.grad,center-self.curpoint())):
                        dist=candist
                        center=point
        return center
        
    def intersect(self,center,n):
        t=np.dot(center-self.start,n)/np.dot(self.grad,n)
        return self.start+self.grad*t
        
    def clamp(self):
        for i,x in enumerate(self.start):
            if abs(x-round(x))<0.0000000000001:
                self.start[i]=round(x)
        for i,x in enumerate(self.grad):
            if abs(x-round(x))<0.0000000000001:
                self.grad[i]=round(x)
    
    def step(self):
        self.update_t()
        self.vprint("\npoint:"+str(self.curpoint().tolist()))
        self.center = self.identify_cube()
        self.vprint("cube:"+str(self.center.tolist()))
        n=self.normalized(np.array(self.fx(*self.center)))
        self.vprint("normal:"+str(n.tolist()))
        pause = True
        if all(n==[0,0,0]):
            if self.verbose:
                self.doprint()
            raise self.Halt()
        if np.dot(n,self.grad)<0:
            point = self.intersect(self.center,n)
            if all(self.center-0.5<=point) and all(self.center+0.5>=point):
                self.vprint("hit reflector")
                self.start = point
                self.t=0
                self.grad = self.grad - 2*np.dot(self.grad,n)*n
                #self.clamp()
                self.vprint("new point:"+str(self.curpoint().tolist()))
                self.vprint("new direction:"+str(self.grad.tolist()))
            else:
                self.vprint("missed reflector")
                self.maybeclearbuffer()
                pause = False
        else:
            self.vprint("missed reflector")
            self.maybeclearbuffer()
            pause = False
        if self.verbose:
            if self.showmiss or pause:
                self.doprint()
                input()
            else:
                if random.randrange(1000)==0:
                    sys.stdout.write(".")
                    sys.stdout.flush()
            
    def run(self):
        self.vprint("start point:"+str(self.curpoint().tolist()))
        self.vprint("start direction:"+str(self.grad.tolist()))
        while True:
            try:
                self.step()
            except self.Halt:
                print("Halted at {} in cube {} with direction {}".format(self.curpoint().tolist(),self.center.tolist(),self.grad.tolist()))
                return
            except KeyboardInterrupt:
                print("Interrupted by user at {} in cube {} with direction {}".format(self.curpoint().tolist(),self.center.tolist(),self.grad.tolist()))
                sys.exit(1)
 
if __name__=="__main__":
    if len(sys.argv)<2:
        print("No filename to execute.")
    else:
        with open(sys.argv[1]) as f:
            listing = f.read()
        ptn = Photon(listing,sys.argv[1])
        ptn.run()