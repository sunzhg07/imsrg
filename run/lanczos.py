from pyIMSRG import *
import numpy as np
from rdm_fun import *
cm=Commutator
#cm.SetIMSRG3Noqqq(True)
#cm.SetIMSRG3Onlyvvv(True)
gm=Generator()



    ## initialize the T and D^dagger


def norm_single(eom, T1, T2):
    return(eom.GetVSEOM_Overlap_single(T1,T2))

def norm_multiref(eom, T1, T2):
    T1d=T1*0.
    T1d.SetHermitian()
    T1d=eom.GetVSEOM_ladder_multiref(T1,0)
    nop=T1*0.
    nop.SetHermitian()
    nop= cm.Commutator(T1d,T2)
    rst = eom.GetVSEOM_Overlap_multiref(nop)
    return(rst/2.)

def norm3_multiref(eom, t1,t2,haml,ms):

    # tested

    rank_j, parity, rank_Tz, particle_rank = 0,0,0,3
    t3 = Operator(ms,rank_j, parity, rank_Tz, particle_rank)
    t3.ThreeBody.SetMode("pn")
    
    t3*=0.
    rst=0.

    nop=t1*0.0
    nop.SetHermitian()
    t3.SetHermitian()

    cm.comm223ss(t2,haml,t3)
    cm.comm232ss(t1,t3,nop)
    cm.comm231ss(t1,t3,nop)
    cm.comm132ss(t1,t3,nop)
    rst = eom.GetVSEOM_Overlap_multiref(nop)
    return(rst/2)




def htc_single(eom,Haml, chi,proj=None):
    ht_plus= chi*0
    ht_plus.SetAntiHermitian()
    ht_plus = cm.Commutator(Haml, chi )
    heom1= eom.GetVSEOM_ladder_single(ht_plus, 0)
    hod = heom1
    return(hod)

def htc_multiref(eom, Haml, chi,prjop=None):

    ht_minus= chi*0
    ht_minus.SetHermitian()  
    ht_minus = cm.Commutator(Haml, chi )
    heom= eom.GetVSEOM_ladder_multiref(ht_minus, 1)
    if(prjop != None):
        prjop(heom)
    return(heom)



import numpy as np

def lanczos_proc( hv_func, norm_func,haml, vi,max_iter,state_want,ms,eom,norm_three=None,prjop=None):
    lanczos_vector = []
    hall = np.zeros([max_iter,max_iter])
    hall[0,0]=0.

    ## normalize it to 1
    nn=norm_func(eom,vi,vi)
    print('initial norm vector', nn)
    vi=vi/np.sqrt(nn)
    lanczos_vector.append(vi)
    norm_e_old=-1000
    norm_e_new=-1000


    for j in range(max_iter):

        w = hv_func(eom,haml,lanczos_vector[j],prjop)
        ai=norm_func(eom, w,lanczos_vector[j])
        hall[j,j]=ai

        if(j>0):
            w=w-ai*lanczos_vector[j]-bj*lanczos_vector[j-1]
        else:
            w=w-ai*lanczos_vector[j]

        nm=norm_func(eom, w,w)
        bj = np.sqrt(nm)
        w1=w/bj
        if(j<max_iter-1):
            hall[j,j+1]=bj
            hall[j+1,j]=bj

        if bj < 0.01 :
            print('bj is small')
            break
        lanczos_vector.append(w/bj)
        if(j > state_want and j%4 == 0):
            e,v = np.linalg.eig(hall[0:j,0:j])
            e=np.sort(e)
            print("Energy on ", j , ' th iteration: ', e[0:state_want] )
            norm_e_new=0.0
            for k in range(state_want):
                norm_e_new += e[k]*e[k]
            if(abs(norm_e_new-norm_e_old) < 0.01):
                print('energy converge')
                break
            norm_e_old=norm_e_new

    print(hall)
    print('Lanczos converged with ', j, ' step')
    for k in range(state_want):
        print('E(',k,')= ', e[k])
    print('Lanczos converged with ', j, ' step')
    return(e[0:state_want], v[0:state_want,:],lanczos_vector)








def read_tdm(tdm_file,ms):


    rank_j, parity, rank_Tz, particle_rank, herm= 0,0,0,2,1
    ops = Operator(ms,rank_j, parity, rank_Tz, particle_rank)
    ops*=0.

    data = []
    with open(tdm_file, 'r') as f:
        lines = f.readlines()
    ## line 0 is the total J of the state
    lidx=0
    line=lines[lidx]
    values = line.strip().split()
    jtotal = float(values[0])
    factor=np.sqrt(2*jtotal +1.)

    ## line 1 is the number of single particle orbits
    lidx+=1
    line=lines[lidx]
    values = line.strip().split()
    norb=int(values[0])

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

    for obs in range(n_obrd):
        lidx+=1
        line=lines[lidx]
        values = line.strip().split()
        aa=ob_idx[int(values[1])-1,0]
        bb=ob_idx[int(values[2])-1,0]
        rd=float(values[-1])/factor
        ops.SetOneBody(aa,bb,rd)

    lidx+=1
    line=lines[lidx]
    values = line.strip().split()
    n_tbrd=int(values[0])


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
        ops.SetTwoBody(jij,pij,tij, jkl,pkl,tkl, aa,bb,cc,dd, rd)

    lidx+=1
    line=lines[lidx]
    values = line.strip().split()
    n_tbrd=int(values[0])

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
        ops.ThreeBody.SetME_pn(jab,jef,jtot,aa,bb,cc,ee,ff,kk,rd)
    return(ops)






def print_op(chi, ms):
    chiv=[]
    pp=0
    tt=0
    for itz in [1]:
        for pp in [0]:
            for jj in [0,2]:
                ch=ms.GetTwoBodyChannelIndex(jj,pp,itz)
                kcf=ms.GetTwoBodyChannel(ch)
                bras=kcf.GetKetIndex_qq()+kcf.GetKetIndex_qv()
                kets=kcf.GetKetIndex_vv()
                for ibra in bras:
                    dbra=kcf.GetKet(ibra)
                    for iket in kets:
                        dket=kcf.GetKet(iket)
                        chiv.append(chi.GetTwoBody(ch,ch,ibra,iket))
    return(np.array(chiv))
    ## ppvh
 #   nch=ms.GetNumberTwoBodyChannels()

 #   for ich in range(nch):
 #       dch = ms.GetTwoBodyChannel(ich)
 #       jj=dch.J
 #       pp=dch.parity
 #       tt=dch.Tz
 #       ch=ich

 #       kcf=dch
 #       if(kcf.GetNumberKets() == 0):
 #           continue
 #       bras=kcf.GetKetIndex_qq()+kcf.GetKetIndex_qv()+kcf.GetKetIndex_vv()
 #       kets=kcf.GetKetIndex_vc()
 #       for ibra in bras:
 #           dbra=kcf.GetKet(ibra)
 #           for iket in kets:
 #               dket=kcf.GetKet(iket)
 #               chiv.append(chi.GetTwoBody(ch,ch,ibra,iket))



def arnoldi_proc(hv_func, norm_func, haml, vi, max_iter, state_want, ms, eom,
                 norm_three=None, rdmat=None, prjop=None):
    """
    Arnoldi eigensolver for H = H1 + H2 where H is hermitian but H1, H2 are not.

    Builds a Krylov subspace using H1 (via hv_func), then computes the full
    H = H1 + H2 matrix in that subspace and diagonalizes it.

    Parameters
    ----------
    hv_func    : callable(eom, haml, v, prjop) -> H1*v   (e.g. htc_multiref)
    norm_func  : callable(eom, v1, v2) -> <v1|v2>        (e.g. norm_multiref)
    haml       : Hamiltonian operator
    vi         : initial vector
    max_iter   : maximum Arnoldi steps
    state_want : number of lowest eigenvalues to target
    ms         : model space
    eom        : EOM object
    norm_three : callable(eom, v1, v2, haml, ms) -> <v1|H2|v2>  (e.g. norm3_multiref)
    rdmat      : unused (reserved for future transition density matrix usage)
    prjop      : projection operator passed to hv_func to remove null space
    """
    def _norm(a, b):
        return norm_func(eom, a, b)

    lanczos_vector = []
    h1v_cache = []          # cache H1*v_j so we don't recompute it
    hall = np.zeros([max_iter, max_iter])

    # Normalize initial vector
    nn = _norm(vi, vi)
    print('arnoldi: initial vector norm =', nn)
    vi = vi / np.sqrt(nn)
    lanczos_vector.append(vi)

    prev_e = None
    tol    = 1e-8
    min_iter = state_want + 1
    e  = np.zeros(state_want)
    vs = np.zeros([1, 1])

    for j in range(max_iter - 1):

        # --- compute H1*v_j ---
        # For matrix elements we need the RAW H1*v (no prjop), matching mr_eom.py.
        # For subspace expansion we will project separately after Gram-Schmidt.
        h1v_j = hv_func(eom, haml, lanczos_vector[j],prjop)   # no prjop
        h1v_cache.append(h1v_j)

        # --- H matrix elements for column j (and row j) ---
        # Both directions are computed independently.
        # H = H1+H2 is hermitian, so hall[i,j] == hall[j,i] must hold exactly.
        # If they differ, something upstream is wrong — print a warning.
        for i in range(j + 1):
            h1ij = _norm(lanczos_vector[i], h1v_j)          # <v_i | H1 | v_j>
            h1ji = _norm(lanczos_vector[j], h1v_cache[i])   # <v_j | H1 | v_i>
            if norm_three is not None:
                h2ij = norm_three(eom, lanczos_vector[i], lanczos_vector[j], haml, ms)
                h2ji = norm_three(eom, lanczos_vector[j], lanczos_vector[i], haml, ms)
            else:
                h2ij = h2ji = 0.0
            hall[i, j] = h1ij + h2ij
            hall[j, i] = h1ji + h2ji
            # Hermiticity check: report H1 and H2 separately to locate the source
            diff_h1 = h1ij - h1ji
            diff_h2 = h2ij - h2ji
            diff_H  = hall[i, j] - hall[j, i]
            scale   = max(1.0, abs(hall[i, j]), abs(hall[j, i]))
            if abs(diff_H) / scale > 1e-6:
                print(f'  HERM BREAK ({i},{j}):  '
                      f'H[i,j]={hall[i,j]:.8f}  H[j,i]={hall[j,i]:.8f}  '
                      f'|diff|={abs(diff_H):.3e}  '
                      f'H1_ij={h1ij:.6f} H1_ji={h1ji:.6f} dH1={diff_h1:.3e}  '
                      f'H2_ij={h2ij:.6f} H2_ji={h2ji:.6f} dH2={diff_h2:.3e}')

        # --- expand subspace: project a copy of h1v_j then Gram-Schmidt ---
        # Use a copy so h1v_cache[j] (used for matrix elements) stays unmodified.
        # Project BEFORE GS (removes null space from the start direction).
        # Project AGAIN AFTER GS (GS subtraction can reintroduce null-space drift
        # via floating-point; without this the euclidean norm blows up over steps).
        w = h1v_j * 1.0
        if prjop is not None:
            prjop(w)
        for i in range(j + 1):
            cij = _norm(lanczos_vector[i], w)
            w   = w - cij * lanczos_vector[i]
        if prjop is not None:
            prjop(w)

        bj = _norm(w, w)
        if abs(bj) < 1e-12:
            print('arnoldi: early breakdown at step', j + 1)
            break
        lanczos_vector.append(w / np.sqrt(bj))

        # --- eigenvalue check / convergence ---
        if j + 1 >= min_iter:
            e, v = np.linalg.eig(hall[:j+1, :j+1])
            idx  = np.argsort(e.real)
            e    = e.real[idx]
            vs   = v[:, idx]

            if (j + 1) % 5 == 0:
                print('arnoldi eigenvalues @ step', j + 1, ':', e[:state_want])

            if prev_e is not None:
                delta = np.max(np.abs(e[:state_want] - prev_e))
                scale = max(1.0, np.max(np.abs(e[:state_want])))
                if delta < tol or delta / scale < tol:
                    print('Arnoldi converged with', j + 1, 'steps')
                    break
            prev_e = e[:state_want].copy()

    print('Arnoldi finished with', j + 1, 'steps')
    for k in range(state_want):
        print('E(', k, ')=', e[k])

    # build Ritz vectors
    ritz_vecs = []
    nb = len(lanczos_vector)
    for k in range(state_want):
        vec = lanczos_vector[0] * 0.0
        for m in range(min(nb, vs.shape[0])):
            vec = vec + float(vs[m, k]) * lanczos_vector[m]
        ritz_vecs.append(vec)

    return e[:state_want], vs, ritz_vecs


def dcom222312(eom,Hs,chi):
    opa=chi*0
    opa.SetHermitian()
    cm.FactorizedDoubleCommutator.SetUse_1b_Intermediates(True);
    cm.FactorizedDoubleCommutator.SetUse_2b_Intermediates(True);
    cm.FactorizedDoubleCommutator.comm223_231(chi, Hs, opa);
    cm.FactorizedDoubleCommutator.comm223_232(chi, Hs, opa);
    cm.FactorizedDoubleCommutator.comm223_132(chi, Hs, opa);
    rst=eom.GetVSEOM_Overlap_multiref(opa)
    return(rst/2,opa)
