
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
    if $argc == 0 || $argc > 2
        help pdna
    else
        #echo \033[01;31
        #echo hoi
        #resets the color
        #echo \033[0m
        if $argc == 1
            printf "[\t+%u]:\t", seq.p >> 1
        else
            printf "[\t%u]:\t", (seq.p >> 1) + $arg1
        end
        call print_seq(&$arg0, KEY_WIDTH)
    end
end


document pseq
	Prints seq.
	Syntax: pseq <seq>
end

define hook-quit
    set confirm off
end

define hook-run
    shell make cleartest
end

#set prompt \033[01;31mgdb$ \033[0m
# Preventing GDB from Pausing during Long Output
set height 0
set width 0

set disassembly-flavor intel

define break_re
    if $argc != 3
        help ngrep
    else
        python break_re_py($arg0, $arg1, $arg2)
    end
end

document break_re
	Prints line in file where match occurs [must be only one].
	Syntax: ngrep <match> <file>
end



# leave this: is for buffers.
#tb key_init.cpp:190
break_re 'seq_t ndx;' 'key_init.cpp' 'tbreak'
commands
    silent
    call print_seq(&seq, KEY_WIDTH)
    call print_seq(&seq, KEY_WIDTH)
    call print_seq(&seq, KEY_WIDTH)
    call print_seq(&seq, KEY_WIDTH)
    call print_seq(&seq, KEY_WIDTH)
    call print_seq(&seq, KEY_WIDTH)
    printf "----------[ start debugging ]------------\n"
    continue
end

#################################################################
# key_init.cpp

#b key_init.cpp:190
break_re 'seq_t ndx;' 'key_init.cpp' 'break'
commands
    silent
    pseq seq h->s_s
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

#################################################################
# fa.cpp


#b fa.cpp:155
break_re 'for(;;) {' 'fa.cpp' 'break'
commands
    silent
    # display caller
    bt 2
    pseq seq b.h->s_s
    printf "prev:%u\tpend:%u\text-1:%u\n", prev, pend, kc->ext + 1
    wa seq.dna
    commands
        silent
        pseq seq
        c
    end
    c
end

#b fa.cpp:245
break_re 'update_header(kc, k, hk, b);' 'fa.cpp' 'break'
commands
    printf "uniq at\t%u", *k >> 1
    c
end

#b fa.cpp:201
#
## if the assigned break was 1:
#
#command 1
#pdna seq 3
#p seq.p + (s - kc->s)
#end



r

