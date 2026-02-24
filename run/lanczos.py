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
    T2d=T1*0.
    T2d.SetAntiHermitian()
    T2d=eom.GetVSEOM_ladder_multiref(T2,1)
    nop=T1*0.
    nop.SetHermitian()
    nop= cm.Commutator(T1,T2d)
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
#    cm.comm232ss(t1,t3,nop)
#    cm.comm231ss(t1,t3,nop)
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

    chi_d=chi*0
    chi_d.SetAntiHermitian()
    chi_d = eom.GetVSEOM_ladder_multiref(chi, 1)
    ht_minus= chi*0
    ht_minus.SetHermitian()  
    ht_minus = cm.Commutator(Haml, chi_d )
    heom= eom.GetVSEOM_ladder_multiref(ht_minus, 0)
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



def arnoldi_proc_new( hv_func, norm_func,norm_three, haml, vi, max_iter,state_want,ms,eom, prjop):
    lanczos_vector = []
    hall = np.zeros([max_iter,max_iter])
    hall[0,0]=0.

    converged = 0
    ## normalize it to 1
    nn=norm_func(eom,vi,vi)
    print('initial norm vector', nn)
    vi=vi/np.sqrt(nn)
    lanczos_vector.append(vi)
    norm_e_old=-10000
    norm_e_new=-1000
    e=np.zeros(state_want)

    for j in range(max_iter-1):
        w = hv_func(eom,haml,lanczos_vector[j],prjop)
        ai=norm_func(eom,w,lanczos_vector[j])
        hall[j,j]=ai

        w=w-ai*lanczos_vector[j]

        if(j>0):
            for i in range(j):
                bj=norm_func(eom,w,lanczos_vector[i])
                w=w-bj*lanczos_vector[i]

        bj = norm_func(eom,w,w)

        if abs(bj) < 0.00001:
            print('coupling bj is small: ', bj, j)
            break
        w=w/np.sqrt(bj)
        lanczos_vector.append(w)
        print('lanczos vector euclidean norm: ',w.Norm())
        nmm=0.
        nm3=0.
        vs=np.zeros([1,1])
        for m in range(j):
            vl=lanczos_vector[m]
            vm=lanczos_vector[j]
            vec=vl*0.0
            vec.SetHermitian()
            vec=hv_func(eom,haml, vm,prjop)
            nmm=norm_func(eom,vl,vec)
            nm3=norm_three(eom,vl,vm,haml, ms)
            hall[m,j] =nmm+nm3
            hall[j,m] =nmm+nm3
        if(j+1 >= state_want ):
            hsub=hall[0:j+1,0:j+1]
            e,v = np.linalg.eig(hsub[0:j+1,0:j+1])
            idx = np.argsort(e)   
            e = e[idx]
            vs = v[:,idx]
            vec_a=lanczos_vector[j]*0.
            for kk in range(len(e)):
                vec_a = vec_a + lanczos_vector[kk]*vs[kk,0]

            vec_b=hv_func(eom,haml, vec_a,prjop)
            nmm=norm_func(eom,vec_a,vec_b)
            nm3=norm_three(eom,vec_a,vec_a,haml, ms)
            nm=norm_func(eom,vec_a,vec_a)

            print('energy check: ', nmm+nm3,nm)
            if(abs(nm-1.)>0.0001):
                print('Lost Orthogonality, converged')
                break
            print("Energy on ", j+1 , ' th iteration: ', e[:] )

            norm_e_new=0.0
            for k in range(state_want):
                norm_e_new += e[k]*e[k]
            if(abs(norm_e_new-norm_e_old) < 0.0000001):
                print('energy converge')
                break
            norm_e_old=norm_e_new


    print('Arnoldi converged with ', j+1, ' step')
    lanczos_vector.pop()
    print(hall)

    return(e, vs,lanczos_vector)



import numpy as np

def arnoldi_solver(hv_func, norm_func, haml, vi, max_iter, state_want, ms, eom,norm_three=None, rdmat=None, prjop=None):
    """
    Final Professional Version:
    - Renamed basis to lanczos_vector.
    - Double Orthogonalization (Twice is Enough) to kill spurious zero-states.
    - Periodic energy reporting.
    """
    tol=1e-3
    res_tol=1e-5
    lanczos_vector = []
    hall = np.zeros([max_iter, max_iter], dtype=complex)
    
    # Normalize initial guess

    nn = norm_func(eom, vi, vi)
    if nn <= 0: 
        raise ValueError("Initial vector has zero or negative norm.")
    
    vi = vi / np.sqrt(nn)
    lanczos_vector.append(vi)
    
    e_vals = np.zeros(state_want, dtype=complex)
    subspace_vs = None
    norm_e_old = -1000.0

    # --- 2. Main Iteration Loop ---
    for j in range(max_iter - 1):
        # A. Apply Operator (Once per step)
        w = hv_func(eom, haml, lanczos_vector[j], prjop)
        
        # B. Modified Gram-Schmidt Orthogonalization
        # Fills the j-th column of the Hessenberg matrix
        for i in range(j + 1):
            h_ij = norm_func(eom, w, lanczos_vector[i])
            hall[i, j] = h_ij
            w = w - h_ij * lanczos_vector[i]
            
        # C. Sub-diagonal element and normalization
        bj_sq = norm_func(eom, w, w)
        bj = np.sqrt(max(0, bj_sq.real))

        if bj < 1e-1:
            print(f"Breakdown at step {j}: Subspace converged.")
            break
            
        if (j + 1) < max_iter:
            hall[j + 1, j] = bj
            
        w_next = w / bj
        lanczos_vector.append(w_next)
        
        # --- 3. Convergence & Reporting ---
        if (j + 1) >= state_want:
            h_sub = hall[0:j+1, 0:j+1]
            h_sub=(h_sub+h_sub.T)/2

            vals, vecs = np.linalg.eig(h_sub)
            
            # Sort energies (real part)
            idx = np.argsort(np.real(vals))
            current_evals = vals[idx]
            current_vs = vecs[:, idx]
            
            # D. Periodic Energy Report (Every 5 steps)
            if (j + 1) % state_want == 0:
                print(f"Iteration {j+1}: Target Energies = {current_evals[:state_want].real}")

            # E. Convergence Check
            norm_e_new = np.sum(np.abs(current_evals[:state_want])**2)
            if abs(norm_e_new - norm_e_old) < tol:
                print(f"Energy converged at iteration {j+1}")
                e_vals, subspace_vs = current_evals, current_vs
                break
            
            norm_e_old = norm_e_new
            e_vals, subspace_vs = current_evals, current_vs

    # --- 4. Final Ritz Vector Construction ---
    physical_states = []
    if subspace_vs is not None:
        for k in range(state_want):
            # Transform coefficients back to full basis
            ritz_vec = lanczos_vector[0] * 0.0
            for i in range(len(lanczos_vector)-1):
                ritz_vec += subspace_vs[i, k].real * lanczos_vector[i]
            physical_states.append(ritz_vec)

    print(f"Arnoldi converged with {len(lanczos_vector)} steps")
    return e_vals[:state_want], physical_states, lanczos_vector



import numpy as np

def subspace_solver(h1_func, norm_func, h2_expect_func, haml, vi, max_basis_size, state_want, ms, eom, rdmat=None ,  prjop=None):
    tol=1e-8
    """
    Subspace expansion for H = H1 + H2.
    Only <v1|H2|v2> is used for H2.
    Elements are cached to prevent redundant computation.
    """
    basis = []
    # Cache for H1 and H2 matrix elements
    h1_cache = np.zeros((max_basis_size, max_basis_size), dtype=complex)
    h2_cache = np.zeros((max_basis_size, max_basis_size), dtype=complex)

    # Normalize initial guess
    nn = norm_func(eom, vi, vi)
    if nn <= 0:
        raise ValueError("Initial vector norm is 0 or negative.")

    v = vi / np.sqrt(nn)

    basis.append(v)

    prev_energy = 0.0

    for j in range(max_basis_size):
        m = len(basis)
        new_idx = m - 1 # The index of the newest vector added

        # 1. Update the Cache for the new row/column
        # We only need to compute the interaction of the NEW vector with all vectors
        h1_v_new = h1_func(eom, haml, basis[new_idx], prjop)

        for i in range(m):
            # Compute H1 part: <basis_i | H1 | basis_new>
            val_h1 =  norm_func(eom, basis[i], h1_v_new)
            h1_cache[i, new_idx] = val_h1
            h1_cache[new_idx, i] = np.conj(val_h1)

            # Compute H2 part: <basis_i | H2 | basis_new>
            val_h2 = h2_expect_func(eom,basis[i], basis[new_idx],haml,ms)
            h2_cache[i, new_idx] = val_h2
            h2_cache[new_idx, i] = np.conj(val_h2)

        # 2. Extract the current active subspace matrix
        H_sub = h1_cache[:m, :m] + h2_cache[:m, :m]

        # 3. Solve the Hermitian eigenvalue problem
        evals, evecs = np.linalg.eigh(H_sub)

        current_energy = evals[0].real
        print(f"Iteration {j+1}: Basis Size = {m}, Energy = {current_energy:.10f}")

        # Convergence check
        if j > 0 and abs(current_energy - prev_energy) < tol:
            print("Converged!")
            break
        prev_energy = current_energy

        # 4. Generate next basis vector using H1 (Krylov expansion)
        # We apply H1 to the most recent vector to explore new directions
        w_next =  h1_func(eom, haml, basis[new_idx], prjop)

        # Modified Gram-Schmidt Orthonormalization
        for b in basis:
            nm=norm_func(eom, b, w_next)
            w_next -= nm * b

        norm=norm_func(eom, w_next, w_next)

        # Handling linear dependency or H1-invariant subspaces
        if norm < 1e-10:
            if m >= max_basis_size: break
            # Inject random noise to find new directions H2 might care about
            w_next = np.random.normal(0, 1, vi.shape).astype(vi.dtype)
            for b in basis:
                nm=norm_func(eom, b, w_next)
                w_next -= nm * b
            norm = norm_func(eom, w_next, w_next)

        basis.append(w_next* ( 1./ norm))

    # Construct the final Ritz vectors (physical states)
    physical_states = []
    for k in range(min(state_want, m)):
        state = basis[0] * 0.0
        for i in range(m):
            state += evecs[i, k] * basis[i]
        physical_states.append(state)

    return evals[:state_want], physical_states


def arnoldi_proc( hv_func, norm_func,haml, vi,max_iter,state_want,ms,eom,norm_three=None,rdmat=None,prjop=None):
    lanczos_vector = []
    hall = np.zeros([max_iter,max_iter])
    hall[0,0]=0.

    converged = 0
    ## normalize it to 1
    nn=norm_func(eom,vi,vi,rdmat)
    print('initial norm vector', nn)
    vi=vi/np.sqrt(nn)
    lanczos_vector.append(vi)
    norm_e_old=-10000
    norm_e_new=-1000
    e=np.zeros(state_want)

    for j in range(max_iter-1):
        w = hv_func(eom,haml,lanczos_vector[j],prjop)
        ai=norm_func(eom,w,lanczos_vector[j])
        hall[j,j]=ai
        w=w-ai*lanczos_vector[j]

        if(j>0):
            for i in range(j):
                bj=norm_func(eom,w,lanczos_vector[i])
                w=w-bj*lanczos_vector[i]


        ## generate the new vector
        bj = norm_func(eom,w,w)
    

        if abs(bj) < 0.00001:
            print('coupling bj is small: ', bj, j)
            break
        w=w/np.sqrt(bj)
        lanczos_vector.append(w)
        print('lanczos vector euclidean norm: ',w.Norm())
        nmm=0.
        vs=np.zeros([1,1])
        for m in range(j):
            vl=lanczos_vector[m]
            vm=lanczos_vector[j]
            vec=vl*0.0
            vec.SetHermitian()
            vec=hv_func(eom,haml, vm,prjop)
            nmm=norm_func(eom,vl,vec)
            hall[m,j] =nmm
            hall[j,m] =nmm
        if(j+1 >= state_want ):
            hsub=hall[0:j+1,0:j+1]
            e,v = np.linalg.eig(hsub[0:j+1,0:j+1])
            idx = np.argsort(e)   
            e = e[idx]
            vs = v[:,idx]
            vec_a=lanczos_vector[j]*0.
            for kk in range(len(e)):
                vec_a = vec_a + lanczos_vector[kk]*vs[kk,0]


            print('energy check: ', nmm)
            if(abs(nmm-1.)>0.0001):
                print('Lost Orthogonality, converged')
                break
            #print("Energy on ", j+1 , ' th iteration: ', e[:] )

            norm_e_new=0.0
            for k in range(state_want):
                norm_e_new += e[k]*e[k]
            if(abs(norm_e_new-norm_e_old) < 0.0000001):
                print('energy converge')
                break
            norm_e_old=norm_e_new


    print('Arnoldi converged with ', j+1, ' step')
    lanczos_vector.pop()
    print(hall)

    return(e, vs,lanczos_vector)





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



def dcom222312(eom,Hs,chi):
    opa=chi*0
    opa.SetHermitian()
#    cm.FactorizedDoubleCommutator.SetUse_1b_Intermediates(True);
#    cm.FactorizedDoubleCommutator.SetUse_2b_Intermediates(True);
#    cm.FactorizedDoubleCommutator.comm223_231(chi, Hs, opa);
#    cm.FactorizedDoubleCommutator.comm223_232(chi, Hs, opa);
    cm.FactorizedDoubleCommutator.comm223_132(chi, Hs, opa);
    rst=eom.GetVSEOM_Overlap_multiref(opa)
    return(rst/2,opa)
