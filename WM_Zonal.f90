program MorphZonal

    implicit none
    
    integer,parameter :: sp=selected_real_kind(p=6)                 ! determine compiler KIND value for 4-byte (single precision) floating point real numbers
    integer,parameter :: dp=selected_real_kind(p=15)                ! determine compiler KIND value for 8-byte (double precision) floating point real numbers
    
    character*1000 :: rundir
    character*1000 :: dump_txt
    character*100 :: rasLW_bin_pth
    character*1000 :: rasZone_bin_pth
    character*1000 :: output_csv_pth
    character*1000 :: overlapping_zone_pth
    character*20 :: n_overlap_str
    character*20 :: noData_str
    character*20 :: nras_str
    character*3 :: S_str
    character*4 :: G_str
    character*2 :: sy_str
    character*2 :: ey_str
    character*2 :: y_str
    character*20 :: fnc_tag
    
    integer :: dump_int
    integer :: nras
    integer :: dem_res
    integer :: noData
    integer :: y,sy,ey
    integer :: i,j,k,iz
    integer :: lt
    integer :: n_overlap
    integer :: nzones
    integer :: nzones_max
    integer :: zone
    integer :: zoneID
    integer :: zoneID_index
    
    integer,dimension(:),allocatable :: rasLW
    integer,dimension(:),allocatable :: rasZone
    integer, dimension(:),allocatable :: zoneIDs
    integer, dimension(:),allocatable :: overlap_zones
    integer, dimension(:,:),allocatable :: overlap
    integer, dimension(:,:),allocatable :: zone_counts
    character*3,dimension(:),allocatable :: lnd_codes               ! array of land-type codes used in output summary file
    

    nzones_max = 1000
    dem_res = 30
    

    call GET_COMMAND_ARGUMENT( 1,rasZone_bin_pth)
    call GET_COMMAND_ARGUMENT( 2,overlapping_zone_pth)
    call GET_COMMAND_ARGUMENT( 3,n_overlap_str)
    call GET_COMMAND_ARGUMENT( 4,nras_str)
    call GET_COMMAND_ARGUMENT( 5,noData_str)
    call GET_COMMAND_ARGUMENT( 6,S_str)
    call GET_COMMAND_ARGUMENT( 7,G_str)
    call GET_COMMAND_ARGUMENT( 8,sy_str)
    call GET_COMMAND_ARGUMENT( 9,ey_str)
    
    rundir = '/ocean/projects/bcs200002p/ewhite12/MP2023/ICM/'//trim(adjustL(S_str))//'/'//trim(adjustL(G_str))//'/geomorph/output'
    fnc_tag = 'MP2023_'//trim(adjustL(S_str))//'_'//trim(adjustL(G_str))//'_C000_U00_V00_SLA'
    
    output_csv_pth = trim(adjustL(rundir))//trim(adjustL(fnc_tag))//'_O_'//trim(adjustL(sy_str))//'_'//trim(adjustL(ey_str))//'_W_MCdir.csv'
        
    read(n_overlap_str,*) n_overlap
    read(nras_str,*) nras
    read(noData_str,*) noData

    allocate(rasLW(nras))
    allocate(rasZone(nras))
    allocate(zoneIDs(nzones_max))
    allocate(overlap_zones(n_overlap))
    allocate(overlap(n_overlap,2))
    
    rasZone         = 0
    zoneIDs         = 0
    overlap_zones   = 0
    overlap         = 0

    ! open output file and write header    
    write(*,'(A,A)') 'Opening output file for writing: ',trim(adjustL(output_csv_pth))
    open(unit=909, file=trim(adjustL(output_csv_pth)), position='append')
    write(909,'(A)') 'prj_no,S,year,code,ElementID,value'
    
    write(*,'(A)') 'Identifying zone IDs for use in zonal summaries.'
    ! open raster file containing map of summary zones
    !open(unit=101, file = trim(adjustL(rasZone_bin_pth)),form='unformatted')
    !read(101) rasZone
    !close(101)
    open(unit=101, file = trim(adjustL(rasZone_bin_pth)) )
    do i=1,nras
        read(101,*) dump_int, dump_int, rasZone(i)
    end do
    close(101)
    
    nzones = 0
    do i=1,nras
        zone = rasZone(i)
        if ( ANY(zoneIDs == zone) ) then
            ! zoneIDs already contains value for current zone - skip ahead
        else
            nzones = nzones + 1
            zoneIDs(nzones) = zone
            !write(*,'(A,I,A,I)') 'ZoneID ',nzones,': ',zone
        end if
        
        if (nzones > nzones_max) then
            write(*,'(A)') '******  FOUND MORE ZONES THAN ALLOWED BY CODE ******'
            write(*,'(A,I,A)') ' Code is limited to ',nzones_max,' summary zones. Exiting now.'
            exit
        end if
    end do

    write(*,'(A)') 'Reading in table of overlapping zoneIDs.'
    open(unit=1111, file=trim(adjustL(overlapping_zone_pth)))
    read(1111,*) dump_txt                                           ! dump header
    do i = 1,n_overlap
        read(1111,*) overlap_zones(i),                      &       ! zoneID equal to the sum of the two overlapping ElementID values
   &                 overlap(i,1),                          &       ! ElementID for first of two overlapping zones
   &                 overlap(i,2)                                   ! ElementID for second of two overlapping zones
    end do
    close(1111)
    
    ! allocate array that will store landtype pixel counts for each zone
    allocate(zone_counts(nzones,5))
    zone_counts = 0
    
    allocate(lnd_codes(5))
    lnd_codes(1) = 'LND'        !  lt=1: vegetated wetland
    lnd_codes(2) = 'WAT'        !  lt=2: water                                           
    lnd_codes(3) = 'BRG'        !  lt=3: unvegetated wetland/new subaerial unvegetated mudflat (e.g., bare ground)
    lnd_codes(4) = 'UPL'        !  lt=4: developed land/upland/etc. that are not modeled in ICM-LAVegMod
    lnd_codes(5) = 'FLT'        !  lt=5: flotant marsh 
    
    read(sy_str,*) sy
    read(ey_str,*) ey
    
    do y = sy,ey
        if (y < 10) then
            y_str = '0'//y
        else
            y_str = ''//y
        end if
                
        rasLW_bin_pth = trim(adjustL(rundir))//trim(adjustL(fnc_tag))//'_N_'//y_str//'_'//y_str//'_W_lndtyp30.xyz.b'
        rasLW = 0
        
        write(*,'(a,a)') 'tabulating zonal counts for:  ', trim(adjustL(rasLW_bin_pth))
        open(unit=100, file = trim(adjustL(rasLW_bin_pth)),form='unformatted')
        read(100) rasLW
        close(100)
        
        ! loop over DEM and analyze landtype for each pixel that is located within a summary zone
        do i=1,nras
            if (rasZone(i) > 0) then
                lt = rasLW(i)
                if (lt /= noData) then
                    zoneID = rasZone(i)
                    zoneID_index = findloc(zoneIDs,zoneID,1)
                    ! check to see if zone is one that is comprised of two overlapping zones
                    ! if so, add the area to both zone counts for the two overlapping zones
                    if ( ANY(overlap_zones == zoneID) ) then
                        iz = findloc(overlap_zones,zoneID,1)
                        do k = 1,2
                            zoneID = overlap(iz,k)
                            zoneID_index = findloc(zoneIDs,zoneID,1)
                            zone_counts(zoneID_index,lt) = zone_counts(zoneID_index,lt) + 1
                        end do
                    ! if zone is not one that is overlapping, simply add the area to the zone count
                    else
                        zone_counts(zoneID_index,lt) = zone_counts(zoneID_index,lt) + 1 
                    end if                                                  
                end if
            end if
        end do
        
        ! write zonal counts to output file
        do i = 1,nzones
            do j = 1,5
                zoneID = zoneIDs(i)
                write(909,1909) trim(adjustL(G_str)),trim(adjustL(S_str)),y,lnd_codes(j),zoneID,zone_counts(i,j)*dem_res*dem_res  
            end do
        end do
    
    end do
    
    close(909)        
    
1909    format(A,',',A,',',I,',',A,',',I,',',I)    
    
end program
