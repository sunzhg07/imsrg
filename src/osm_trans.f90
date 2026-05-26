subroutine trans_rdme(lv,rv,op_jexp)
    use mpi_mapping
    use parallel
    use lanczos
    use m_sp
    use constants
    use tiny_functions
    use storage
    use m_configurations
    use ang_mom_functions
    use file_list
    use operators
    use m_osm_basis_new
    use truncations

    implicit none
    integer:: lv,rv
    integer:: bra,ket, number_channels,ch,bra_dim
    real*8:: coef, phase_ab,phase_cd
    integer:: iph
    real*8:: dij
    integer:: itz,p_parity,ang_mom,a,b,c,d,nconfs
    integer:: aa,bb,cc,dd, ja,jb,jc,jd,ma,mb,mc,md
    integer:: jmin,jmax, jt, channel, bra1,ket1
    integer:: jab
    integer:: e,f,k,je,jf,jk,me,mf,mk,ee,ff,kk
    integer:: mef, jef,mab,i,j,ich
    real*8::angmon_fact,val
    real*8,allocatable,dimension(:,:):: amat
    !!  integer:: idx(4)=(/2,4,3,5/)
    type(m_configuration_descriptor)::bra_configs
    type(oper_jm), intent(inout) :: op_jexp

    ! New loop variables
    integer :: t, p_i, p_k, n_j, n_l
    integer :: dist_p, dist_n
    integer :: mz_idx, mz_idx_ket
    integer :: p_bra_start, p_bra_end, n_bra_start, n_bra_end
    integer :: p_ket_start, p_ket_end, n_ket_start, n_ket_end
    integer :: p_bra_idx, n_bra_idx, p_ket_idx, n_ket_idx
    integer :: bra_global, ket_global
    integer :: mz_p_val, mz_n_val, mz_p_ket_val, mz_n_ket_val
    
    if (iam == 0) then
       write(16,*)'!reduced transition density matrix'
       write(16,*)'!one body, Twobody follow the same convention as kshell'
       write(16,*)'!threebody is normalized, for energy : sqrt(ang_mom+1.)*rdm_j%v3b*v_int_j%v3b '
       write(16,*)"Computing densitymatrix of state: |", lv,'>'

       write(6,*)"Computing transition bettwen state: ", lv,rv
       call flush(6)
    endif
    
    call init_mscheme_operator(rdm_m)
    call init_jscheme_operator(rdm_j)


    ! Partition-based loop for Transition Density Matrix with Task Distribution
    ! Loop over assigned tasks (distributed proton partition pairs)
    do t = 1, size(my_tasks)
       p_i = my_tasks(t)%pi
       p_k = my_tasks(t)%pk
       dist_p = partition_map_proton(p_i, p_k)
       
       ! Loop over neutron partition pairs
       do n_j = 1, size(partition_blocks_n)
          if (pb_mbasis_meta(p_i, n_j, 3) == -1) cycle 

          do n_l = 1, size(partition_blocks_n)
             if (pb_mbasis_meta(p_k, n_l, 3) == -1) cycle
             
             dist_n = partition_map_neutron(n_j, n_l)
             if (dist_p + dist_n > 3) cycle
             
             ! Iterate over bra states in (p_i, n_j) block
             bra_global = pb_mbasis_meta(p_i, n_j, 1) - 1
          
             do mz_idx = 1, partition_blocks_p(p_i)%n_mz_blocks
                mz_p_val = partition_blocks_p(p_i)%mz_values(mz_idx)
                mz_n_val = jtot - mz_p_val
                
                p_bra_start = lookup_pid_p(p_i, mz_p_val, 1)
                p_bra_end   = lookup_pid_p(p_i, mz_p_val, 2)
                n_bra_start = lookup_pid_n(n_j, mz_n_val, 1)
                n_bra_end   = lookup_pid_n(n_j, mz_n_val, 2)
                
                if (p_bra_start == 0 .or. n_bra_start == 0) cycle

                do p_bra_idx = p_bra_start, p_bra_end
                   if (mbsp_valid(p_bra_idx)==0) cycle
                   do n_bra_idx = n_bra_start, n_bra_end
                      if (mbsn_valid(n_bra_idx)==0) cycle
                      if (trunc .and. .not.check_occ(p_bra_idx, n_bra_idx)) cycle

                      bra_global = bra_global + 1
                      
                      ! For each bra state, iterate over ket states
                      ! Start index check for triangle condition optimization
                      if (pb_mbasis_meta(p_k, n_l, 2) < bra_global) cycle

                      ket_global = pb_mbasis_meta(p_k, n_l, 1) - 1
                      
                      do mz_idx_ket = 1, partition_blocks_p(p_k)%n_mz_blocks
                         mz_p_ket_val = partition_blocks_p(p_k)%mz_values(mz_idx_ket)
                         mz_n_ket_val = jtot - mz_p_ket_val
                         
                         p_ket_start = lookup_pid_p(p_k, mz_p_ket_val, 1)
                         p_ket_end   = lookup_pid_p(p_k, mz_p_ket_val, 2)
                         n_ket_start = lookup_pid_n(n_l, mz_n_ket_val, 1)
                         n_ket_end   = lookup_pid_n(n_l, mz_n_ket_val, 2)
                         
                         if (p_ket_start == 0 .or. n_ket_start == 0) cycle
                         
                         do p_ket_idx = p_ket_start, p_ket_end
                            if (mbsp_valid(p_ket_idx)==0) cycle
                            do n_ket_idx = n_ket_start, n_ket_end
                               if (mbsn_valid(n_ket_idx)==0) cycle
                               if (trunc .and. .not.check_occ(p_ket_idx, n_ket_idx)) cycle
                               
                               ket_global = ket_global + 1
                               
                               if (ket_global < bra_global) cycle
                               
                               coef = wvfs(bra_global, lv) * wvfs(ket_global, rv)
                               if (abs(coef) > 1.0d-20) then
                                  call dmoverlap(mbsp(p_bra_idx), mbsn(n_bra_idx), mbsp(p_ket_idx), mbsn(n_ket_idx), coef)
                               endif
                               
                            enddo ! n_ket_idx
                         enddo ! p_ket_idx
                      enddo ! mz_idx_ket
                   enddo ! n_bra_idx
                enddo ! p_bra_idx
             enddo 
          enddo
       enddo
    enddo ! End task loop

    ! MPI Reduction: Sum rdm_m contributions from all processes
    call mpi_allreduce(mpi_in_place, rdm_m%v1b, size(rdm_m%v1b), &
                       mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
    
    ! Reduce 2-body density matrix
    do ich = 1, num_m2_channels
       if (allocated(rdm_m%v2b(ich)%val)) then
          call mpi_allreduce(mpi_in_place, rdm_m%v2b(ich)%val, size(rdm_m%v2b(ich)%val), &
                             mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
       endif
    enddo
    
    ! Reduce 3-body density matrix (if rank_rdm == 3)
    if (rank_rdm == 3) then
       do ich = 1, num_m3_channels
          if (allocated(rdm_m%v3b(ich)%val3)) then
             call mpi_allreduce(mpi_in_place, rdm_m%v3b(ich)%val3, size(rdm_m%v3b(ich)%val3), &
                                mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
          endif
       enddo
    endif 

    ! Only master process performs post-processing and output
    if (iam == 0) then




    do bra= 1,jjsp%total_orbits
        do ket= 1,jjsp%total_orbits
            if(jjsp%itzp(bra) /= jjsp%itzp(ket))cycle
            if(jjsp%jj(bra) /= jjsp%jj(ket))cycle
            if(jjsp%ll(bra) /= jjsp%ll(ket))cycle
            rdm_m%v1b(bra,ket)=rdm_m%v1b(bra,ket)/sqrt(jjsp%jj(bra)+1.)
            if(abs(rdm_m%v1b(bra,ket)) > 1e-6)&
                write(16,'(A6,2(I8,2x),f13.6)')"OBTD ", bra,ket, rdm_m%v1b(bra,ket)*(sqrt(jtot+1.))
        enddo
    enddo


    !!! antisymmetrize rdm_m%v2b
    do ich=1,num_m2_channels
        itz=m2_channel_list(ich,1)
        p_parity=m2_channel_list(ich,2)
        ang_mom=m2_channel_list(ich,3)

        nconfs=m2_channel_list(ich,4)

        allocate(amat(nconfs,nconfs))

        amat=0.

        do bra=1,mab_configs(ich)%number_confs
            a=mab_configs(ich)%config_ab(2*bra-1)
            b=mab_configs(ich)%config_ab(2*bra)
            if(a>b)cycle
            bra1=locate_m2_channel(b,a,2)
            amat(bra1,:)=-rdm_m%v2b(ich)%val(bra,:)
        enddo
        rdm_m%v2b(ich)%val = amat+rdm_m%v2b(ich)%val

        amat=0.
        do bra=1,mab_configs(ich)%number_confs
            a=mab_configs(ich)%config_ab(2*bra-1)
            b=mab_configs(ich)%config_ab(2*bra)
            if(a>b)cycle
            bra1=locate_m2_channel(b,a,2)
            amat(:,bra1)=-rdm_m%v2b(ich)%val(:,bra)
        enddo

        rdm_m%v2b(ich)%val = amat+rdm_m%v2b(ich)%val

        deallocate(amat)
    enddo



    do ich=1,num_m2_channels
        itz=m2_channel_list(ich,1)
        p_parity=m2_channel_list(ich,2)
        ang_mom=m2_channel_list(ich,3)

        nconfs=m2_channel_list(ich,4)  
        do bra=1,mab_configs(ich)%number_confs
            a=mab_configs(ich)%config_ab(2*bra-1)
            b=mab_configs(ich)%config_ab(2*bra)
            aa=all_orbits%jord(a)
            bb=all_orbits%jord(b)
            ma=all_orbits%mm(a)
            mb=all_orbits%mm(b)
            ja=jjsp%jj(aa)
            jb=jjsp%jj(bb)

            do ket =1, mab_configs(ich)%number_confs
                c=mab_configs(ich)%config_ab(2*ket-1)
                d=mab_configs(ich)%config_ab(2*ket)
                cc=all_orbits%jord(c)
                dd=all_orbits%jord(d)
                mc=all_orbits%mm(c)
                md=all_orbits%mm(d)
                jc=jjsp%jj(cc)
                jd=jjsp%jj(dd)

                if(abs(rdm_m%v2b(ich)%val(bra,ket))<1e-6)cycle

                jmin = max(abs(ja-jb)/2, abs(jc-jd)/2)
                jmax = max(abs(ja+jb)/2, abs(jc+jd)/2)

                do jt= jmin,jmax

                    if((aa==bb .or. cc==dd) .and. mod(jt,2)==1)cycle
                    channel = locate_j2_channel(itz, p_parity, jt)
                    if(channel ==0)cycle
                    bra1 = lookup_ab_configs(channel)%ival(aa,bb)
                    ket1 = lookup_ab_configs(channel)%ival(cc,dd)
                    if(bra1*ket1 ==0 )cycle
                    rdm_j%v2b(channel)%val(bra1,ket1) = rdm_j%v2b(channel)%val(bra1,ket1) +&
                        rdm_m%v2b(ich)%val(bra,ket)*cgc(ja,jb,jt*2,ma,mb,ang_mom*2)*&
                        cgc(jc,jd,jt*2,mc,md,ang_mom*2)/dij(aa,bb)/dij(cc,dd)

                enddo
            enddo
        enddo
    enddo





    do ich=1,num_j2_channels
        itz=j2_channel_list(ich,1)
        p_parity=j2_channel_list(ich,2)
        ang_mom=j2_channel_list(ich,3)
        nconfs=j2_channel_list(ich,4)


        rdm_j%v2b(ich)%val=rdm_j%v2b(ich)%val/sqrt(ang_mom*2.+1.)

        do bra=1,nconfs
            a=ab_configs(ich)%config_ab(2*bra-1)
            b=ab_configs(ich)%config_ab(2*bra)
            ja=jjsp%jj(a)
            jb=jjsp%jj(b)
            phase_ab = iph((ja+jb)/2+ang_mom+1)
            do ket =1, nconfs
                c=ab_configs(ich)%config_ab(2*ket-1)
                d=ab_configs(ich)%config_ab(2*ket)
                jc=jjsp%jj(c)
                jd=jjsp%jj(d)
                phase_cd = iph((jc+jd)/2+ang_mom+1)
                if(a>b)cycle
                if(c>d)cycle
                if(abs(rdm_j%v2b(ich)%val(bra,ket)) < 1e-8)cycle
                if(a<=b .and. c<=d)then
                    write(16,'(A6, 5(I4, 1x),f13.6)')'TBTD',a,b,c,d,ang_mom, &
                        rdm_j%v2b(ich)%val(bra,ket)*sqrt(jtot+1.)
                endif
            enddo
        enddo
    enddo

    !!! fully antisymitrize 3b in mscheme
    do ich=1,num_m3_channels
        itz=m3_channel_list(ich,1)
        p_parity=m3_channel_list(ich,2)
        ang_mom=m3_channel_list(ich,3)
        nconfs=m3_channel_list(ich,4)


        allocate(amat(nconfs,nconfs))

        amat=0.d0


        do  bra=1,nconfs
            i=mabc_configs(ich)%config_abc(3*bra-2)
            j=mabc_configs(ich)%config_abc(3*bra-1)
            k=mabc_configs(ich)%config_abc(3*bra)
            if(i>=j)cycle
            if(j>=k)cycle
            bra1=locate_m3_channel(i,k,j,2)
            amat(bra1,:) = -rdm_m%v3b(ich)%val(bra,:)
            bra1=locate_m3_channel(k,i,j,2)
            amat(bra1,:) = rdm_m%v3b(ich)%val(bra,:)
            bra1=locate_m3_channel(j,i,k,2)
            amat(bra1,:) = -rdm_m%v3b(ich)%val(bra,:)
            bra1=locate_m3_channel(j,k,i,2)
            amat(bra1,:) = rdm_m%v3b(ich)%val(bra,:)
            bra1=locate_m3_channel(k,j,i,2)
            amat(bra1,:) = -rdm_m%v3b(ich)%val(bra,:)
        enddo

        rdm_m%v3b(ich)%val = rdm_m%v3b(ich)%val + amat
        amat=0.d0
        do  bra=1,nconfs
            i=mabc_configs(ich)%config_abc(3*bra-2)
            j=mabc_configs(ich)%config_abc(3*bra-1)
            k=mabc_configs(ich)%config_abc(3*bra)

            if(i>=j)cycle
            if(j>=k)cycle

            bra1=locate_m3_channel(i,k,j,2)
            amat(:,bra1) = -rdm_m%v3b(ich)%val(:,bra)
            bra1=locate_m3_channel(k,i,j,2)
            amat(:,bra1) = rdm_m%v3b(ich)%val(:,bra)
            bra1=locate_m3_channel(j,i,k,2)
            amat(:,bra1) = -rdm_m%v3b(ich)%val(:,bra)
            bra1=locate_m3_channel(j,k,i,2)
            amat(:,bra1) = rdm_m%v3b(ich)%val(:,bra)
            bra1=locate_m3_channel(k,j,i,2)
            amat(:,bra1) = -rdm_m%v3b(ich)%val(:,bra)
        enddo

        rdm_m%v3b(ich)%val = rdm_m%v3b(ich)%val + amat
        deallocate(amat)
    enddo



do ich=1,num_m3_channels
    itz=m3_channel_list(ich,1)
    p_parity=m3_channel_list(ich,2)
    ang_mom=m3_channel_list(ich,3)
    nconfs=m3_channel_list(ich,4)
    do  bra=1,nconfs
        a=mabc_configs(ich)%config_abc(3*bra-2)
        b=mabc_configs(ich)%config_abc(3*bra-1)
        c=mabc_configs(ich)%config_abc(3*bra)
        aa=all_orbits%jord(a)
        bb=all_orbits%jord(b)
        cc=all_orbits%jord(c)
        ja=all_orbits%jj(a)
        jb=all_orbits%jj(b)
        jc=all_orbits%jj(c)
        ma=all_orbits%mm(a)
        mb=all_orbits%mm(b)
        mc=all_orbits%mm(c)
        mab=ma+mb


        do ket = 1,nconfs
            e=mabc_configs(ich)%config_abc(3*ket-2)
            f=mabc_configs(ich)%config_abc(3*ket-1)
            k=mabc_configs(ich)%config_abc(3*ket)
            je=all_orbits%jj(e)
            jf=all_orbits%jj(f)
            jk=all_orbits%jj(k)
            me=all_orbits%mm(e)
            mf=all_orbits%mm(f)
            mk=all_orbits%mm(k)
            ee=all_orbits%jord(e)
            ff=all_orbits%jord(f)
            kk=all_orbits%jord(k)
            mef=me+mf


            do jab = abs(ja-jb), ja+jb, 2
                do jef= abs(je-jf), je+jf,2

                    jmin = max(abs(jab-jc),abs(jef-jk))
                    jmax = min(jab+jc, jef+jk)


                    do jt = jmin, jmax ,2

                        channel = locate_j3_channel(itz,p_parity,jt)
                        if(channel==0)cycle
                        bra1=lookup_abc_configs(channel)%ival4(aa,bb,cc,jab)
                        ket1=lookup_abc_configs(channel)%ival4(ee,ff,kk,jef)
                        if(bra1*ket1==0)cycle
                        angmon_fact=cgc(ja,jb,jab,ma,mb,mab)*cgc(jab,jc,jt, mab,mc, ang_mom)*&
                            cgc(je,jf,jef,me,mf,mef)*cgc(jef,jk,jt,mef,mk,ang_mom)/sqrt(jt+1.)
                        rdm_j%v3b(channel)%val(bra1,ket1)= rdm_j%v3b(channel)%val(bra1,ket1) + angmon_fact*&
                            rdm_m%v3b(ich)%val(bra,ket)
                    enddo
                enddo
            enddo
        enddo
    enddo
enddo






    !! check the density matrix

    !! onebody
    val=0.d0
    do bra =1,jjsp%total_orbits
        do ket = 1,jjsp%total_orbits
            val = val + rdm_m%v1b(bra,ket)*op_jexp%v1b(bra,ket)*sqrt(jjsp%jj(bra)+1.)
        enddo
    enddo

!    write(6,*)'1b contribute to energy: ', val
!    call flush(6)


    do ich=1,num_j2_channels
        itz=j2_channel_list(ich,1)
        p_parity=j2_channel_list(ich,2)
        ang_mom=j2_channel_list(ich,3)
        nconfs=j2_channel_list(ich,4)

        do bra=1,nconfs
            a=ab_configs(ich)%config_ab(2*bra-1)
            b=ab_configs(ich)%config_ab(2*bra)
            if(a>b)cycle
            do ket=1,nconfs
                c=ab_configs(ich)%config_ab(2*ket-1)
                d=ab_configs(ich)%config_ab(2*ket)
                if(c>d)cycle
                if(abs(rdm_j%v2b(ich)%val(bra,ket)*op_jexp%v2b(ich)%val(bra,ket)) <1e-6)cycle

                val = val + rdm_j%v2b(ich)%val(bra,ket)*op_jexp%v2b(ich)%val(bra,ket)*sqrt(2.*ang_mom+1.)
                !              write(*,'(5(I4,1x),2(f13.6,1x))')a,b,c,d,ang_mom, rdm_j%v2b(ich)%val(bra,ket),v_int_j%v2b(ich)%val(bra,ket)

            enddo
        enddo
    enddo

!    write(6,*)'1b+2b contribute to energy: ', val
!    call flush(6)



    do ich=1,num_j3_channels
        itz=j3_channel_list(ich,1)
        p_parity=j3_channel_list(ich,2)
        ang_mom=j3_channel_list(ich,3)
        nconfs=j3_channel_list(ich,4)
        do  bra=1,nconfs
            a=abc_configs(ich)%config_abc(4*bra-3)
            b=abc_configs(ich)%config_abc(4*bra-2)
            c=abc_configs(ich)%config_abc(4*bra-1)
            jab=abc_configs(ich)%config_abc(4*bra)
            if(a>b)cycle
            if(b>c)cycle
            do  ket=1,nconfs
                e=abc_configs(ich)%config_abc(4*ket-3)
                f=abc_configs(ich)%config_abc(4*ket-2)
                k=abc_configs(ich)%config_abc(4*ket-1)
                jef=abc_configs(ich)%config_abc(4*ket)
                if(e>f)cycle
                if(f>k)cycle
                if(abs(rdm_j%v3b(ich)%val(bra,ket))< 1.e-8)cycle
                val = val + sqrt(ang_mom+1.)*rdm_j%v3b(ich)%val(bra,ket)*op_jexp%v3b(ich)%val(bra,ket)/dabc(a,b,c)/dabc(e,f,k)

                write(16,'(A8,10(I2,1x),f13.6,1x)')'TRBTD',a,b,c,e,f,k,jab,jef,ang_mom,itz,&
                    rdm_j%v3b(ich)%val(bra,ket)/dabc(a,b,c)/dabc(e,f,k)
            enddo
        enddo
    enddo

    write(16,*)"END densitymatrix of state: |", lv,'>'

    write(6,*)'1b+2b+3b contribute to energy: ', val
    call flush(6)

    endif ! End iam == 0
    
    ! All processes clean up their operators
    call finalize_mscheme_operator(rdm_m)
    call finalize_jscheme_operator(rdm_j)

contains
    real*8 function dabc(a,b,c)
        integer:: a,b,c
        dabc=1.
        if(a==b .and. b==c)then
            dabc=6.
            return
        elseif(a==b .or. b==c)then
            dabc=2.
            return
        else
            dabc=1.
            return
        endif
    end function dabc

end subroutine trans_rdme
