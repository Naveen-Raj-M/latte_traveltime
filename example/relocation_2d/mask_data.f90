
program test

    use libflit

    implicit none

    real, allocatable, dimension(:) :: tp, ts
    integer :: i
    integer :: ns, nr
    integer, allocatable, dimension(:) :: rp, rs
    real :: delta = 3.0

    ! Mask data
    call make_directory('./data_masked')

    ns = 40
    nr = 101
    do i = 1, ns

        tp = load('./data/shot_'//num2str(i)//'_traveltime_p.bin', nr)
        ts = load('./data/shot_'//num2str(i)//'_traveltime_s.bin', nr)

        rp = random_permute(regspace(1, 1, nr), seed=i)
        rp = rp(1:irand(range=[2, nint(nr/2.0)], seed=100 + i))

        rs = random_permute(regspace(1, 1, nr), seed=200 + i)
        rs = rs(1:irand(range=[2, nr], seed=300 + i))

        tp(rp) = -1
        ts(rs) = -1

        call output_array(tp, './data_masked/shot_'//num2str(i)//'_traveltime_p.bin')
        call output_array(ts, './data_masked/shot_'//num2str(i)//'_traveltime_s.bin')

        print *, i

    end do


    ! Mask more data
    call make_directory('./data_masked_more')

    ns = 40
    nr = 101
    do i = 1, ns

        tp = load('./data/shot_'//num2str(i)//'_traveltime_p.bin', nr)
        ts = load('./data/shot_'//num2str(i)//'_traveltime_s.bin', nr)

        rp = random_permute(regspace(1, 1, nr), seed=500 + i)
        rp = rp(1:irand(range=[nint(0.5*nr), nint(0.75*nr)], seed=1000 + i))

        rs = random_permute(regspace(1, 1, nr), seed=2000 + i)
        rs = rs(1:irand(range=[nint(0.75*nr), nint(0.9*nr)], seed=3000 + i))

        tp(rp) = -1
        ts(rs) = -1

        call output_array(tp, './data_masked_more/shot_'//num2str(i)//'_traveltime_p.bin')
        call output_array(ts, './data_masked_more/shot_'//num2str(i)//'_traveltime_s.bin')

        print *, i

    end do


    ! Mask data
    call make_directory('./data_masked_no_t0')

    ns = 40
    nr = 101
    do i = 1, ns

        tp = load('./data/shot_'//num2str(i)//'_traveltime_p.bin', nr)
        ts = load('./data/shot_'//num2str(i)//'_traveltime_s.bin', nr)

        tp = tp - minval(tp) + delta
        ts = ts - minval(ts) + delta

        rp = random_permute(regspace(1, 1, nr), seed=i)
        rp = rp(1:irand(range=[2, nint(nr/2.0)], seed=100 + i))

        rs = random_permute(regspace(1, 1, nr), seed=200 + i)
        rs = rs(1:irand(range=[2, nr], seed=300 + i))

        tp(rp) = -1
        ts(rs) = -1

        call output_array(tp, './data_masked_no_t0/shot_'//num2str(i)//'_traveltime_p.bin')
        call output_array(ts, './data_masked_no_t0/shot_'//num2str(i)//'_traveltime_s.bin')

        print *, i

    end do


    ! Mask more data
    call make_directory('./data_masked_more_no_t0')

    ns = 40
    nr = 101
    do i = 1, ns

        tp = load('./data/shot_'//num2str(i)//'_traveltime_p.bin', nr)
        ts = load('./data/shot_'//num2str(i)//'_traveltime_s.bin', nr)

        tp = tp - minval(tp) + delta
        ts = ts - minval(ts) + delta

        rp = random_permute(regspace(1, 1, nr), seed=500 + i)
        rp = rp(1:irand(range=[nint(0.5*nr), nint(0.75*nr)], seed=1000 + i))

        rs = random_permute(regspace(1, 1, nr), seed=2000 + i)
        rs = rs(1:irand(range=[nint(0.75*nr), nint(0.9*nr)], seed=3000 + i))

        tp(rp) = -1
        ts(rs) = -1

        call output_array(tp, './data_masked_more_no_t0/shot_'//num2str(i)//'_traveltime_p.bin')
        call output_array(ts, './data_masked_more_no_t0/shot_'//num2str(i)//'_traveltime_s.bin')

        print *, i

    end do

end program test
