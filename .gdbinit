
# not needed, already done
#source ~/.gdbinit

# back to default behavior. When gdb cannot find the breakpoint location, it queries whether a pending is ok


#break -source fa.cpp -function get_next_nt_seq  if p >= (kc->s_l << 2)

#ASSERT(seq.p < pend, return -EFAULT, "%lu >= %lu?", seq.p, pend);\
#ASSERT(pend < (kc->s_l << 2), return -EFAULT, "%lu/%lu?", pend, kc->s_l);\
#break -source fa.cpp -function 

define pdna
    if $argc != 1
        help pdna
    else
        call print_seq(&$arg0, KEY_WIDTH)
    end
end

define pseq
    if $argc == 0 || $argc > 2
        help pdna
    else
        if $argc == 1
            printf "[\t+%u]:\t", seq.p
        else
            printf "[\t%u]:\t", seq.p + $arg1
        end
        call print_seq(&$arg0, KEY_WIDTH)
    end
end


document pdna
	Prints dna.
	Syntax: pvector <seq>
end

define hook-run
shell make cleartest
end

#tb main
#r

## test key_init.c
#wa seq.dna

b key_init.cpp:190

commands 1
silent
pseq seq h->s_s
c
end
r





#b fa.cpp:201
#
## if the assigned break was 1:
#
#command 1
#pdna seq 3
#p seq.p + (s - kc->s)
#end
