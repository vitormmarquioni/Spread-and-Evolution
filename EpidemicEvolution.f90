PROGRAM epidemic

IMPLICIT NONE

INTEGER :: i,j,k,L,indiv,m,k1,k2
INTEGER :: N,B,initial,total,N2
INTEGER :: numModules
INTEGER :: seed(12)
INTEGER :: susceptibleTot,infectedTot,recoveredTot
INTEGER, ALLOCATABLE :: susceptible(:),infected(:),recovered(:),modul(:)
INTEGER, ALLOCATABLE :: sus(:),inf(:),recov(:),dea(:),sizeModules(:)
INTEGER, ALLOCATABLE :: virus(:)!,dist(:,:)
INTEGER :: dist
INTEGER, ALLOCATABLE :: adj(:,:),genome(:,:),how(:),how2(:),indtime(:)
REAL, ALLOCATABLE :: incub(:),prob(:),rprob(:)
REAL :: p,r,d
REAL*8 :: mu
REAL :: aux
character(100) :: network
character(4) :: string

OPEN(UNIT=7,FILE='auxi.in',STATUS='OLD',POSITION='REWIND')
	!model parameters
	READ(7,*)N       !population size
	READ(7,*)B       !genome size (times 2)
	READ(7,*)total   !time horizon
	READ(7,*)p       !transmission probability	
	READ(7,*)r       !recovery probability
	READ(7,*)mu      !mutation rate
	READ(7,*)initial !initial infected individual
	READ(7,*)network !network -> adjacency matrix file
	
	CLOSE(7)
OPEN(UNIT=14,FILE='OUT4.txt',STATUS='UNKNOWN')
OPEN(UNIT=17,FILE='ancestor.txt',STATUS='UNKNOWN')
OPEN(UNIT=18,FILE='extinction.txt',STATUS='UNKNOWN')
OPEN(UNIT=19,FILE='phylo.out',STATUS='UNKNOWN')
OPEN(UNIT=9,FILE='auxi2.in',STATUS='OLD',POSITION='REWIND')

	READ(9,*)numModules
	ALLOCATE(sizeModules(numModules))
	ALLOCATE(susceptible(numModules))
	ALLOCATE(infected(numModules))
	ALLOCATE(recovered(numModules))
	Do i=1,numModules
		READ(9,*)sizeModules(i)
	End do
	
CLOSE(9)

!-----------------------------------------------------------
OPEN(UNIT=7,FILE='seed.in',STATUS='OLD',POSITION='REWIND')
	READ(7,*)seed    !random seed "put"
CLOSE(7)

write(*,*)seed

!--------------------incubation times-----------------------
ALLOCATE(incub(N))
OPEN(UNIT=7,FILE='incub.in',STATUS='OLD',POSITION='REWIND')
Do i=1,N
	READ(7,*)incub(i)
end do
CLOSE(7)
!-----------------------------------------------------------

write(*,*) '----------------------------------------- '
write(*,*)'r = ',r
write(*,*)'p = ',p
write(14,*) '----------------------------------------- '
write(14,*)'r = ',r
write(14,*)'p = ',p
!-----------------------------------------------------------

CALL RANDOM_SEED(put=seed)

Allocate(adj(N,N),genome(N,B))
Allocate(how(N),how2(N),indtime(N),prob(N),rprob(N))
Allocate(sus(N),inf(N),recov(N),dea(N),modul(N))

!----------defining the module of each individual-----------
k1=1
k2=sizeModules(1)
Do i=1,numModules
	do j=k1,k2
		modul(j)=i
	end do
	k1=k2+1
	k2=k2+sizeModules(i+1)
end do
!-----------------------------------------------------------

how = 1!individual state vector
indtime = 0!individual clock vector - the clock starts after it gets exposed

! how(individual) = 1 : Susceptible
! how(individual) = 0 : Infected
! how(individual) = 2 : Recovered

!---------------------Defining the network------------------
Write(*,*)'Defining the network'

OPEN(UNIT=8,File=TRIM(adjustl(network)),status='old',position='rewind')

Do i=1,N-1
	adj(i,i)=0
	do j=i+1,N
		READ(8,*)adj(i,j)
		adj(j,i)=adj(i,j)
	end do
end do
adj(n,n)=0
!-----------------------------------------------------

!Defining the initial virus genome
genome = 1

!-----------------------------------------------------
!Defining the initial conditions

susceptibleTot=N-1
infectedTot=1
recoveredTot=0

susceptible = sizeModules 
infected = 0
recovered = 0

susceptible(modul(initial))=susceptible(modul(initial))-1
infected(modul(initial))=+1

Write(17,*)initial,0,0

Do i=1,numModules
	write(string,*)i
	Open(unit=20+i,File='susMod'//Trim(adjustl(string))//'.out',status='UNKNOWN',position='rewind')
		write(20+i,*)float(susceptible(i))*100./float(sizeModules(i))
		write(*,*)float(susceptible(i))*100./float(sizeModules(i))
	Open(unit=30+i,File='infMod'//Trim(adjustl(string))//'.out',status='UNKNOWN',position='rewind')
		write(30+i,*)float(infected(i))*100./float(sizeModules(i))
		write(*,*)float(infected(i))*100./float(sizeModules(i))
	Open(unit=40+i,File='recMod'//Trim(adjustl(string))//'.out',status='UNKNOWN',position='rewind')
		write(40+i,*)float(recovered(i))*100./float(sizeModules(i))
		write(*,*)float(recovered(i))*100./float(sizeModules(i))
end do

Write(11,*)float(susceptibleTot)*100./float(N)
Write(12,*)float(infectedTot)*100./float(N)
Write(13,*)float(recoveredTot)*100./float(N)

!-----------------------------------------------------
prob=0
rprob=0
how(initial)=0
how2=how!updated state vector

Do m=1,N
	Write(100,*)how(m)
End Do

Do L=0,total-1
	Write(*,*)'infecting time = ',L
	Do indiv=1,N
		if (how(indiv)==0) then
			If (indtime(indiv).GE.incub(indiv))then
				prob(indiv)=p
				rprob(indiv)=r
			end if
			Do j=1,N
				if (how(j)*how2(j)*adj(indiv,j).NE.0) then
					CALL RANDOM_NUMBER(aux)
					
					if (aux.LE.prob(indiv)) then
						how2(j)=0
						infected(modul(j))=infected(modul(j))+1
						susceptible(modul(j))=susceptible(modul(j))-1
						infectedTot=infectedTot+1
						susceptibleTot=susceptibleTot-1
						Write(19,*)indiv,'->',j
						Write(17,*)j,indiv,L+1
						!copying the virus
						Do k=1,B
							genome(j,k)=genome(indiv,k)
						End Do
					End if
				End if
			End Do
			
			CALL RANDOM_NUMBER(aux)
			!attempt to recovery
			if (aux.LE.rprob(indiv)) then
				how2(indiv)=2
				recovered(modul(indiv))=recovered(modul(indiv))+1
				infected(modul(indiv))=infected(modul(indiv))-1
				recoveredTot=recoveredTot+1
				infectedTot=infectedTot-1
				Write(18,*)L+1,indiv
				Do i=1,N
					adj(indiv,i)=0
					adj(i,indiv)=0
				end do
				else
					Do k=1,B/2
						CALL RANDOM_NUMBER(aux)
						if (aux.LE.mu) then
							CALL RANDOM_NUMBER(aux)
							if (aux.LE.(1./3.)) then
								genome(indiv,2*k)=1-genome(indiv,2*k)
								genome(indiv,2*k-1)=1-genome(indiv,2*k-1)
								else
									if (aux.LE.(2./3.)) then
										genome(indiv,2*k)=1-genome(indiv,2*k)
										else
											genome(indiv,2*k-1)=1-genome(indiv,2*k-1)
									end if
							end if
						end if
					End Do
					indtime(indiv)=indtime(indiv)+1
			end if
		End if
	End Do
	how=how2
	
	Do m=1,N
		Write(100+L+1,*)how(m)
	End Do
	
	Do m=1,numModules
		write(string,*)m
		write(20+m,*)float(susceptible(m))*100./float(sizeModules(m))
		write(30+m,*)float(infected(m))*100./float(sizeModules(m))
		write(40+m,*)float(recovered(m))*100./float(sizeModules(m))
	end do
	
	Write(11,*)float(susceptibleTot)*100./float(N)
	Write(12,*)float(infectedTot)*100./float(N)
	Write(13,*)float(recoveredTot)*100./float(N)
	!-----------------------------------------------------
	!!Total number of virus within the population
	N2=infectedTot+recoveredTot
	!-----------------------------------------------------
	allocate(virus(N2))

	j=1
	Do i=1,N
		If (how(i).NE.1) then
			virus(j)=i
			j=j+1
		end if
	end do


	
	!
	dist = 0
	!---------Calculating the pairwise distances----------
	!---------among the viruses only at instants----------
	!--------------that are multiples of 15---------------
	if (mod(L,15).EQ.0) then
	write(*,*)'Calculating the Distances'
		Do i=1,N2-1
			write(*,*)i,' / ',N2
			k1=virus(i)
			Do j=i+1,N2
				k2=virus(j)
				dist=B/2
				DO k=1,B/2
					dist=dist-(abs(genome(k1,2*k)-genome(k2,2*k))-1)*(abs(genome(k1,2*k-1)-genome(k2,2*k-1))-1)
				End DO
				write(1000+L,*)dist
			End DO
		End Do
	end if
	!-----------------------------------------------------
	deallocate(virus)
End Do

Close(19)

write(*,*) '----------------------------------------- '
Write(*,*)float(susceptibleTot)*100./float(N),'% susceptible individuals'
Write(*,*)float(infectedTot)*100./float(N),'% infected individuals'
Write(*,*)float(recoveredTot)*100./float(N),'% recovered individuals'
write(*,*) '----------------------------------------- '

write(*,*)'non-susceptible = ',N2
write(*,*) '----------------------------------------- '

write(14,*) '----------------------------------------- '
Write(14,*)float(susceptibleTot)*100./float(N),'% susceptible individuals'
Write(14,*)float(infectedTot)*100./float(N),'% infected individuals'
Write(14,*)float(recoveredTot)*100./float(N),'% recovered individuals'
write(14,*) '----------------------------------------- '

write(14,*)'non-susceptible = ',N2
write(14,*) '----------------------------------------- '


!---------------Updating the random seed------------------
CALL RANDOM_SEED(get=seed)
OPEN(UNIT=7,FILE='seed.in',STATUS='OLD',POSITION='REWIND')
	WRITE(7,*)seed    !random seed "put"
CLOSE(7)

CLOSE(14)

END PROGRAM epidemic