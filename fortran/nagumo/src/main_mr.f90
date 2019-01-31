      Program MAIN_MR

      use mod_common
      use mod_structure
      use mod_initial
      use mod_arbre
      use mod_recherche
      use mod_integration
      use mod_reaction
      use mod_diffusion

      IMPLICIT NONE

      CHARACTER(len=80) :: nom_fich

      INTEGER :: i, j, k
      INTEGER :: l, m, comp_f

      INTEGER, DIMENSION(3) :: ind_pos
      INTEGER, DIMENSION(ndim) :: position

      DOUBLE PRECISION, DIMENSION(ndim) :: x_deb, x_fin
      DOUBLE PRECISION :: t1, t2, t_cpu, t3, t4, t_reac, t5, t6,&
                          t_reac1, t_reac2, t_diff, t7, t8

      TYPE(structure_maille), Pointer :: root, courant


! ====================== Initialisation ============================

! ..... Initialisation du temps de calcul

       call cpu_time(time=t1)

! ======================= lecture des donnees =========================

      call lect_don
     
! .....Initialisation d'epsilon

     epsilon = eps_max

!.... Pour resolution splittee sans MR

      if (dabs(epsilon-0D0).le. zero) then
       allocate(sol_split(nvar*nbre_maille,nint(tps_fin/dtinc)+1))
       allocate(err_max(nvar),error_l2(nvar))
      endif

! ========================== Initialisations ==========================

! .....Pas de temps de splitting	

          dt = tps_fin/dble(npdt)


      Print*, ' '
      Print*, ' Splitting timestep : ', dt

! .....etendu du support pour la graduation
!!      sgrad = s + (nordre+1)/2
!!      sgrad = s + (nordre)/2 + 1
       !!sgrad = s + 1
       sgrad = s

      !Print*, ' '
      !Print*, ' Graduation order : ', sgrad

! .....pas d'espace sur le niveau le plus fin

      do l = 1,ndim
        dxmin(l) = ( lim_fin(l) - lim_deb(l) ) / &
                   ( dble(2**niv_max * nbre_racine(l)) ) 
      enddo

      !Print*, ' '
      !Print*, ' Dx_min : ', dxmin


! +++++Initialisation des coefficients d'interpolation selon l'ordre

! .....Allocation tableaux des coefficients d'interpolation d'ordre 2*s
      ALLOCATE(coef(s))

! .....ordre s = 1
      if(s == 1) then

        coef(1)=-0.125d0

! .....ordre s = 2
      elseif(s == 2) then

        coef(1)=-22.d0/128.d0
        coef(2)=3.d0/128.d0

! .....ordre s = 3
      elseif(s == 3) then

        coef(1)=-201.d0/1024.d0
        coef(2)=11.d0/256.d0
        coef(3)=-5.d0/1024.d0

      endif


! +++++Initialisation des racines de l'arbre

! .....Allocation des racines
      ALLOCATE( racine(nbre_racine(1),nbre_racine(2), nbre_racine(3)) )

! ============== Definition de l'arbre et de la solution ==============



! +++++Creation et chainage de l'arbre

      do k = 1, nbre_racine(3)
        ind_pos(3) = k
        do j = 1, nbre_racine(2)
          ind_pos(2) = j
          do i = 1, nbre_racine(1)
            ind_pos(1) = i

            do l = 1, ndim
              position(l) = ind_pos(l)
            enddo
            call chainage(racine(i,j,k)%ptr, niv_max, 0, position, courant)

          enddo
        enddo
      enddo

! +++++Creation et chainage du support

! .....creation du support pour la racine
 

      call support_racine(racine)



! .....support pour tous l'arbre


      do k = 1, nbre_racine(3)
        do j = 1, nbre_racine(2)
          do i = 1, nbre_racine(1)

            root => racine(i, j, k)%ptr
            call chainage_support(root)

           enddo
        enddo
      enddo


! ++++++Definition du maillage

      do k = 1, nbre_racine(3)
        ind_pos(3) = k
        do j = 1, nbre_racine(2)
          ind_pos(2) = j
          do i = 1, nbre_racine(1)
            ind_pos(1) = i

            do l = 1, ndim
              x_deb(l) = lim_deb(l) + (lim_fin(l) - lim_deb(l)) * &
                         (dble(ind_pos(l) - 1) / dble(nbre_racine(l)))

              x_fin(l) = lim_deb(l) + (lim_fin(l) - lim_deb(l)) * &
                         (dble(ind_pos(l)) / dble(nbre_racine(l)))
            enddo
            root => racine(i, j, k)%ptr
            call multimesh(root, x_deb, x_fin)

          enddo
        enddo
      enddo
       
!..............................................................
!..............BOUCLE EN EPSILON...............................
!..............................................................
      Do while (epsilon .ge. eps_min .or. dabs(epsilon-0D0).le. zero)

        if (dabs(epsilon-0D0).gt. zero .and. allocated(sol_split)) then
          write(nom_fich,1002) niv_max,s,epsilon
          open(33,file=nom_fich,status='unknown')
          call cpu_time(time=t1)
        endif


! +++++Initialisation de la solution sur le niveau le plus fin
        
! .....Initialisation du temps de simulation
        temps=0.d0
        i_err = 0

 
        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
              call init_sol(root)

            enddo
          enddo
        enddo



!+++++Initialisation du solveur RADAU5

      LENWKREAC = 4*nvar*nvar + 12*nvar + 20
      LENIWKREAC = 3*nvar + 20
      
      allocate( ELWRKREAC(LENWKREAC),  IELWRKREAC(LENIWKREAC))

      ELWRKREAC(1:20) = 0.d0
      IELWRKREAC(1:20) = 0

!+++++Initialisation du solveur ROCK4

      IELWRKDIFF(1) = 0 ! 0 -> ROCK calculates max eigenvalue
      IELWRKDIFF(2) = 0 ! 1 -> the Jacobian is constant
      IELWRKDIFF(3) = 0 
      IELWRKDIFF(4) = 0

 
! .............................................................
! .......................Boucle en temps.......................
! .............................................................
        
      Do nt = 1, npdt


       if (dabs(epsilon - 0d0) .gt. zero) then 
! .....intialisation de arbre = .false.

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr

              call init_flags(root)

            enddo
          enddo
        enddo


! +++++Codage par valeurs moyennes 
!      u connue par valeurs moyennes u(i,niv_max) sur grille la plus fine

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr

              call encoding_moy(root)

            enddo
          enddo
        enddo


! +++++Calcul des valeurs de details sur tout l'arbre

        det_max = 0.d0

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
  
              call encoding_det(root)

            enddo
          enddo
        enddo



! ============== Generation de la grille hybride initiale ============

! +++++Construction de l'arbre hybride

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
  
              call seuillage(root)

            enddo
          enddo
        enddo


! +++++Construction de l'arbre hybride

        do
          compteur_create = 0

! .....verification de la graduation et creation si necessaire
          do k = 1, nbre_racine(3)
            do j = 1, nbre_racine(2)
              do i = 1, nbre_racine(1)

                root => racine(i, j, k)%ptr
  
                call graduation(root)

              enddo
            enddo
          enddo

! .....valorisation ==> decodage par valeur moyenne

          do k = 1, nbre_racine(3)
            do j = 1, nbre_racine(2)
              do i = 1, nbre_racine(1)

                root => racine(i, j, k)%ptr
  
                call decoding_moy(root)

              enddo
            enddo
          enddo


          if( compteur_create == 0 ) exit
        enddo


! +++++Elagage de l'arbre hybride

          do k = 1, nbre_racine(3)
            do j = 1, nbre_racine(2)
              do i = 1, nbre_racine(1)

                racine(i, j, k)%ptr%arbre = .true.

                root => racine(i, j, k)%ptr

                call elagage(root)

              enddo
            enddo
          enddo

! +++++Etablissament du support sur le nouvel arbre

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
              call chainage_support(root)

            enddo
          enddo
        enddo

      endif
! +++++Etablissement de la liste des feuilles

        compteur_feuille = 0

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
              call liste_feuille(root, compteur_feuille)

            enddo
          enddo
        enddo
 

      if (dabs(epsilon - 0d0) .gt. zero) then 
    
! +++++Etablissament des feuilles fictives : pour calcul de la diffusion

        compteur_fictive = 0

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
  
              call feuille_fictive(root)

            enddo
          enddo
        enddo

! +++++Etablissament du support sur le nouvel arbre

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
              call chainage_support(root)

            enddo
          enddo
        enddo

       if (save_res == 1 .or. dabs(temps+dt-tps_fin).le.zero) then

   
!  Compute new details to save them

! +++++Codage par valeurs moyennes 
!      u connue par valeurs moyennes u(i,niv_max) sur grille la plus fine

            do k = 1, nbre_racine(3)
               do j = 1, nbre_racine(2)
                  do i = 1, nbre_racine(1)

                   root => racine(i, j, k)%ptr

                   call encoding_moy(root)

                  enddo
               enddo
            enddo


! +++++Calcul des valeurs de details sur tout l'arbre

             det_max = 0.d0

             do k = 1, nbre_racine(3)
                do j = 1, nbre_racine(2)
                   do i = 1, nbre_racine(1)

                    root => racine(i, j, k)%ptr
  
                    call encoding_det(root)

                   enddo
                enddo
             enddo



         endif 


       endif
! +++++Ecriture de la solution initiale 

        if( nt == 1 ) then

          if (save_res == 1) then
             !!write(nom_fich,999) int(temps),niv_max,s,epsilon
             write(nom_fich,999) int(temps)
             open(31,file=nom_fich,status='unknown')
          endif

!--------x,y,...,dx,dy,..,a,b,c

          do k = 1, nbre_racine(3)
            do j = 1, nbre_racine(2)
              do i = 1, nbre_racine(1)

                root => racine(i, j, k)%ptr
                call print_leaves(root)

              enddo
            enddo
          enddo

        if (save_res == 1) close(31)
            print*
            !!Print*,' Time : ',temps
            write(*,'(a, f10.6)') '  Time : ', temps
            Print*, ' Number of leaves : ', compteur_feuille, &
                    ' Compression rate: ', &
                  100-100.*dble(compteur_feuille)/dble(nbre_maille),' %'


        endif

! +++++integration sur les feuilles uniquement non fictives

! .....Allocation du vector de feuilles phi

      ALLOCATE( tab_phi(compteur_feuille) )

! .....Construction du vector phi

              comp_f = 0

              do k = 1, nbre_racine(3)
                do j = 1, nbre_racine(2)
                  do i = 1, nbre_racine(1)
  
                    root => racine(i, j, k)%ptr
                    call construction_phi(root, comp_f)

                  enddo
                enddo
              enddo


! .....Premier pas de reaction

       call cpu_time(time=t3)

! .....Integration : RADAU5

!..... Integration de la rection a partir du vecteur de feuilles

       call integration_reaction

       call cpu_time(time=t4)
       
       t_reac1 = t4 - t3

!======================================================================

! .....Pas de temps de splitting (DIFFUSION)	

! .....Integration temporelle -> ROCK4

      call cpu_time(time=t7)

      call integration_diffusion

      call cpu_time(time=t8)
       
       t_diff = t8 - t7

! .....Fin du traitement de la diffusion

!======================================================================

! .....Deuxieme pas de reaction

! .....Integration : RADAU5

       call cpu_time(time=t5)

!..... Integration de la rection a partir du vecteur de feuilles

       call integration_reaction

       call cpu_time(time=t6)
       
       t_reac2 = t6 - t5


! .....Desallocation du vector de variables phi

       DEALLOCATE(tab_phi)

! .....Accumulation du temps de simulation
     
          temps = temps + dt

! +++++Ecriture et sauvegarde de la solution aux temps intermediaires

          if( (temps >= tps_deb ) .or. &
              (abs(temps-tps_deb) <= zero) .or. temps == tps_fin) then

            print*
            !!Print*,' Time = ',temps
            write(*,"(a, f10.6)") '  Time : ', temps
            Print*, ' Number of leaves : ', compteur_feuille, &
                    ' Compression rate: ', &
                 100.-100.*dble(compteur_feuille)/dble(nbre_maille),' %'

            t_reac = t_reac1 + t_reac2 

            Print*, ' Integration time REACTION   = ',t_reac, ' sec.'
            Print*, ' Integration time DIFFUSION  = ',t_diff, ' sec.'

! .....ecriture des feuilles reelles
            i_err = i_err + 1


      if (save_res == 1 .or. dabs(temps-tps_fin).le.zero .or.& 
           temps.gt.tps_fin .or. dabs(epsilon-0D0).le. zero) then

         if (save_res == 1 .or. dabs(temps-tps_fin).le.zero .or. temps.gt.tps_fin) then

            write(nom_fich,999) int(100.*temps)
            !!write(nom_fich,999) int(100.*temps),niv_max,s,epsilon

            open(31,file=nom_fich,status='unknown')
    
         endif 

!--------x,y,...,dx,dy,..,a,b,c
            i_sav = 0 ! used to save split solution for eps=0
          
            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)
  
                  root => racine(i, j, k)%ptr
                  call print_leaves(root)
     
                enddo
              enddo
            enddo

          if (save_res == 1 .or. dabs(temps-tps_fin).le.zero .or. temps.gt.tps_fin)  close(31)

        endif

!... reconstruction au niveau le plus fin 

           if ( (eps_max .ne. eps_min .and.  temps == tps_fin) .or. &
                (dabs(epsilon-0D0).gt. zero) .and. (allocated(sol_split))  ) then


! +++++Codage par valeurs moyennes 
!      u connue par valeurs moyennes u(i,niv_max) sur grille la plus fine

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr
                  call encoding_moy(root)

                enddo
              enddo
            enddo


        do
          compteur_create = 0


            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)
  
                  root => racine(i, j, k)%ptr
                  call built_fine_sol(root)
     
                enddo
              enddo
            enddo
!          print*,'comp create = ',compteur_create,&
!               'fine-leaves = ',nbre_maille-compteur_feuille

          if( compteur_create == 0 ) exit
        enddo

        endif



!... calcul de l'erreur MR

      if ( dabs(epsilon-0D0).gt. zero .and. allocated(sol_split)  ) then

! .....ecriture des erreurs L2

         if (save_res == 1 .or. dabs(temps-tps_fin).le.zero) then 

            write(nom_fich,1001) int(100.*temps),niv_max,s,epsilon

            open(32,file=nom_fich,status='unknown')

         endif

            err_max(:) = 0d0
            error_l2(:) = 0d0 
            i_sav = 0

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)
  
                  root => racine(i, j, k)%ptr
                  call compute_error(root)
     
                enddo
              enddo
            enddo

!.........write t, L2 error, L infinity error, data compression
           write(33,900) temps,(dsqrt(error_l2(m)*dxmin(1)),m=1,nvar),&
              (err_max(m),m=1,nvar), 100.*dble(compteur_feuille)/dble(nbre_maille)

          if (save_res == 1 .or. dabs(temps-tps_fin).le.zero)  close(32)

           endif


            tps_deb=tps_deb+dtinc
    
          endif

! .....Fin de boucle en temps

      Enddo

! ..... Finalisation du temps de calcul

       call cpu_time(time=t2)
       
       t_cpu = t2 - t1

       Print*, ' '
       Print*, '================================================= '
       Print*, ' Time CPU   = ',t_cpu, ' sec.'
       !!Print*, ' Threshold value   = ',epsilon

       deallocate(ELWRKREAC,IELWRKREAC)

! .....Mis a jour d'epsilon
       if (epsilon .eq. 0d0) then
          epsilon = 0.1d0
       else
          epsilon = epsilon/10d0
       endif
       !!print *, "LSLSLSLS: epsilon", epsilon

       if ( dabs(epsilon-0D0).gt. zero .and. allocated(sol_split)  ) close(33)
       tps_deb=0d0

      Enddo

!====================================================================

!......FORMATS

1000  format('Sol_init_',i8.8,'.dat')

!!999   format('Champs_t',i4.4,'_niv',i2.2,'_s',i1,'_eps',1pe6.0,'.dat')
999   format('Champs_t',i4.4,'.dat')

1001  format('Error_t',i4.4,'_niv',i2.2,'_s',i1,'_eps',1pe6.0,'.dat')

1002  format('L2_error_niv',i2.2,'_s',i1,'_eps',1pe6.0,'.dat')

900   format(4(1x,1pe15.8))


      end program MAIN_MR
