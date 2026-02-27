from pyIMSRG import *
import numpy as np
class rdm_el:
    def __init__(self):
        self.n1b=0
        self.n2b=0
        self.n3b=0
        self.jtotal=0
        self.norb=0
        self.rd1b=[]
        self.rd2b=[]
        self.rd3b=[]
        self.ms=0
    def read_tdm(self,tdm_file,ms):
        rank_j, parity, rank_Tz, particle_rank, herm= 0,0,0,2,1
        self.ms=ms

        with open(tdm_file, 'r') as f:
            lines = f.readlines()
        ## line 0 is the total J of the state
        lidx=0
        line=lines[lidx]
        values = line.strip().split()
        self.jtotal = float(values[0])
        factor=np.sqrt(2*self.jtotal +1.)

        ## line 1 is the number of single particle orbits
        lidx+=1
        line=lines[lidx]
        values = line.strip().split()
        norb=int(values[0])
        self.norb=norb
        ob_idx=np.zeros([norb,3],dtype=np.int8)

        for obs in range(norb):
            lidx+=1
            line=lines[lidx]
            values = line.strip().split()
            nn=int(values[1])
            ll=int(values[2])
            jj=int(values[3])
            tt=int(values[4])
            ips = ms.GetOrbitIndex(nn,ll,jj,tt)
            ob_idx[obs,0]=ips
            ob_idx[obs,1]=ll
            ob_idx[obs,2]=tt
        ## read in OBTD

        lidx+=1
        line=lines[lidx]
        values = line.strip().split()
        n_obrd=int(values[0])
        self.n1b=n_obrd

        for obs in range(n_obrd):
            lidx+=1
            line=lines[lidx]
            values = line.strip().split()
            aa=ob_idx[int(values[1])-1,0]
            bb=ob_idx[int(values[2])-1,0]
            rd=float(values[-1])/factor
            self.rd1b.append([aa,bb,rd])

        lidx+=1
        line=lines[lidx]
        values = line.strip().split()
        n_tbrd=int(values[0])
        self.n2b=n_tbrd


        for obs in range(n_tbrd):
            lidx+=1
            line=lines[lidx]
            values = line.strip().split()
            aa=ob_idx[int(values[1])-1,0]
            bb=ob_idx[int(values[2])-1,0]
            cc=ob_idx[int(values[3])-1,0]
            dd=ob_idx[int(values[4])-1,0]
            jij=int(values[5])
            pij=ob_idx[int(values[1])-1,1]+ob_idx[int(values[2])-1,1]
            tij=ob_idx[int(values[1])-1,2]+ob_idx[int(values[2])-1,2]
            pij=int(pij%2)
            tij=int(tij/2)
            jkl=int(values[5])
            pkl=ob_idx[int(values[3])-1,1]+ob_idx[int(values[4])-1,1]
            tkl=ob_idx[int(values[3])-1,2]+ob_idx[int(values[4])-1,2]
            pkl=int(pkl%2)
            tkl=int(tkl/2)
            rd=float(values[-1])/factor
            self.rd2b.append([jij,pij,tij, jkl,pkl,tkl, aa,bb,cc,dd, rd])

        lidx+=1
        line=lines[lidx]
        values = line.strip().split()
        n_tbrd=int(values[0])
        self.n3b=n_tbrd

        for  obs in range(n_tbrd):
            lidx+=1 
            line = lines[lidx]
            values = line.strip().split()
            aa=ob_idx[int(values[1])-1,0]
            bb=ob_idx[int(values[2])-1,0]
            cc=ob_idx[int(values[3])-1,0]
            ee=ob_idx[int(values[4])-1,0]
            ff=ob_idx[int(values[5])-1,0]
            kk=ob_idx[int(values[6])-1,0]
            jab=int(int(values[7])/2)
            jef=int(int(values[8])/2)
            jtot=int(values[9])
            rd=float(values[-1])/factor
            self.rd3b.append([jab,jef,jtot,aa,bb,cc,ee,ff,kk,rd])

    def GetVSEOM_Overlap_rd(self,nop):
        ## one body
        rst=0
        rst1=0
        rst2=0
        rst3=0
        nme = self.n1b
        for i in range(nme):
            mel=self.rd1b[i]
            val=nop.GetOneBody(mel[0], mel[1])
            obr=self.ms.GetOrbit(mel[0])
            rst1+= val*mel[2]*np.sqrt((obr.j2+1))
        nme = self.n2b
        for i in range(nme):
            mel=self.rd2b[i]
            j_bra=mel[0]
            p_bra=mel[1]
            t_bra=mel[2]
            j_ket=mel[3]
            p_ket=mel[4]
            t_ket=mel[5]
            a=mel[6]
            b=mel[7]
            c=mel[8]
            d=mel[9]
            val=nop.TwoBody.GetTBME_J (j_bra, j_ket, a, b, c, d)
            norm_fact=1.
            if(a==b):
                norm_fact/=np.sqrt(2)
            if(c==d):
                norm_fact/=np.sqrt(2)

#            print(a,b,c,d,j_bra,val,mel[-1],'new')
            rst2+= val*mel[-1]*np.sqrt(2*j_bra+1)*norm_fact
        
        
        nme = self.n3b
        for i in range(nme):
            mel=self.rd3b[i]
            Jab_in =int(mel[0])
            Jde_in =int(mel[1])
            j0     =mel[2] 
            j1     =mel[2] 
            a      =mel[3] 
            b      =mel[4] 
            c      =mel[5] 
            d      =mel[6] 
            e      =mel[7] 
            f      =mel[8] 
            val=nop.ThreeBody.GetME_pn (Jab_in, Jde_in, j0, a, b, c, d, e, f)
            rst3+= val*mel[-1]**np.sqrt(j0+1)
        return(rst1+rst2+rst3)
