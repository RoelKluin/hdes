
# not needed, already done
#source ~/.gdbinit

set breakpoint pending on
break -source fa.cpp -function handle_range -label uq_1st_occurance if *sk > k
break -source fa.cpp -function extd_uq_by_p -label reeval_keys_in_scope \
  if kc->kct[*ndxkct] & DUP_BIT

# back to default behavior. When gdb cannot find the breakpoint location, it queries whether a pending is ok
set breakpoint pending auto

break -source fa.cpp -function uq_bounded       if prevk && prev >= p
break -source fa.cpp -function handle_range     if prev >= p
break -source fa.cpp -function get_nextnt       if p >= (kc->s_l << 2)
break -source fa.cpp -function extd_uq_by_p     if pend >= (kc->s_l << 2) || p >= pend
#break -source fa.cpp -function 

define hook-stop
shell make cleartest
end
