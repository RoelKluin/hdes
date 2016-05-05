
source gdb_py/frame_filter.py

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

document pdna
	Prints dna.
	Syntax: pdna <dna (or number)>
end


define pseq
    if $argc != 0
        help pseq
    else
        printf "[\t+%u]:\t", seq.p
        call print_seq(&seq, KEY_WIDTH)
    end
end

document pseq
	Prints seq.
	Syntax: pseq
end

#define hook-quit
#    set confirm off
#end

define hook-run
    shell make cleartest
end

#set prompt \033[01;31mgdb$ \033[0m

define break_re
    if $argc != 3 && $argc != 4
        help break_re
    else
        if $argc == 3
            python break_re_py($arg0, $arg1, $arg2)
        else
            python break_re_py($arg0, $arg1, $arg2, $arg3)
        end
    end
end

document break_re
	Insert specified break in file where match occurs [must be only one].
	Syntax: break_re <match> <file> <break|tbreak> [command]
end



# leave this: is for buffers.
#tb key_init.cpp:190
#break_re '_addtoseq(kc->s, seq.t); // kc->s_l grows here' 'key_init.cpp' 'tbreak'
#commands
#    silent
#    call print_seq(&seq, KEY_WIDTH)
#    call print_seq(&seq, KEY_WIDTH)
#    call print_seq(&seq, KEY_WIDTH)
#    call print_seq(&seq, KEY_WIDTH)
#    call print_seq(&seq, KEY_WIDTH)
#    call print_seq(&seq, KEY_WIDTH)
#    printf "----------[ start debugging ]------------\n"
#    continue
#end
#
#################################################################
# key_init.cpp

#b key_init.cpp:190
break_re '_addtoseq(kc->s, seq.t); // kc->s_l grows here' 'key_init.cpp' 'break'
commands
    silent
    pseq
    c
end

#TODO use grep -n to get this:
#cscope -L -0 main main.cpp

#b key_init.cpp:223
break_re 'if (kc->s_l != h->s_s) {' 'key_init.cpp' 'break'
commands
    silent
    if corr == 0
        printf "Header "
        # because printf states: Cannot access memory at address ...
        output  kc->id[h->part[0]]
        printf "\n"
    else
        if kc->s_l != h->s_s
            printf "N stretch of %u\n", corr
        else
            printf "N stretch at start of contig (not inserted)\n"
        end
    end
    c
end

break_re 'kc->uqct += kc->kct_l;' 'key_init.cpp' 'tbreak'
commands
    silent
    print show_mantras(kc)
    c
end


#################################################################
# fa.cpp


define handle_non_uniques
    #bt 2
    #printf "prev:%u\tseq.p:%u\tpend:%u\text-1:%u\t", prev, seq.p, pend, kc->ext + 1
    pseq
end

#b fa.cpp:157
break_re 'for (;;) {' 'fa.cpp' 'tbreak'
commands
    silent
    handle_non_uniques
    wa seq.dna
    commands
        silent
        pseq
        #c
    end
    break_re 'for(;;) {' 'fa.cpp' 'break'
    commands
        silent
        handle_non_uniques
        #c
    end
end

#b fa.cpp:251
break_re 'update_header(kc, k, hk, b);' 'fa.cpp' 'break'
commands
    silent
    printf "uniq at\t%u\n", *k >> 1
    #c
end

#b fa.cpp:201
#
## if the assigned break was 1:
#
#command 1
#pdna seq 3
#p seq.p + (s - kc->s)
#end

#b fa.cpp:200
break_re '// next contig' 'fa.cpp' 'break'
commands
    silent
    printf "header update (can be late?):\t"
    #c
end

#b fa.cpp:216
break_re '// next boundary' 'fa.cpp' 'break'
commands
    silent
    printf "boundary update:\t"
    #c
end

break_re '// stored first occurance matches current position and contig' 'fa.cpp' 'break'
commands
    silent
    printf "excised:\t"
    pseq
    #c
end


define bt_p_locals
    bt 1
    info locals
end
break_re '// GDB$' 'fa.cpp' 'break' 'bt_p_locals'








#leave:
r

