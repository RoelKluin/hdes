
source gdb_py/frame_filter.py

# not needed, already done
#source ~/.gdbinit

# back to default behavior. When gdb cannot find the breakpoint location, it queries whether a pending is ok


#break -source fa.cpp -function get_next_nt_seq  if p >= (kc->s_l << 2)

#ASSERT(seq.p < pend, return -EFAULT, "%lu >= %lu?", seq.p, pend);\
#ASSERT(pend < (kc->s_l << 2), return -EFAULT, "%lu/%lu?", pend, kc->s_l);\
#break -source fa.cpp -function


define pbuf
    if $argc != 2
        help pbuf
    else
        set $i = $arg0[$arg1]
        printf "%c%c%c%c\n", b6(($i&3)<<1), b6(($i&0xc)>>1), b6(($i&0x30)>>3), b6(($i&0xc0)>>5)
    end
end

document pbuf
	Prints buf as dna quad.
	Syntax: pbuf <buffer> <offset>
end

define pkct
    if $argc != 0
        call print_kct(kc, b, $arg0)
    else
        call print_kct(kc, b, 0)
    end
end

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

#define hook-run
#    shell make cleartest
#end


#define hook-quit
#    set confirm off
#end

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
b fa_read

#b key_init.cpp:190
break_re 'uint32_t ndx;' 'key_init.cpp' 'break'
commands
    silent
    pseq
    c
end

#TODO use grep -n to get this:
#cscope -L -0 main main.cpp

#b key_init.cpp:223
break_re '// key after header/stretch to be rebuilt' 'key_init.cpp' 'break'
commands
    silent
    if corr == 0
        printf "Header "
        # because printf states: Cannot access memory at address ...
        output  kc->id[h->ido]
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
    pkct kc->kct + *contxt_idx
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
break_re '// update contig' 'fa.cpp' 'break'
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
break_re '// update contig' 'fa.cpp' 'break'
commands
    silent
    printf "header update (can be late?):\t"
    #c
end

#b fa.cpp:216
break_re '// update assembly' 'fa.cpp' 'break'
commands
    silent
    printf "boundary update:\t"
    #c
end

break_re '*contxt_idx = kc->kct_l++;//GDB:1' 'fa.cpp' 'break'
commands
    silent
    pkct k
end

break_re '*contxt_idx = kc->kct_l++;//GDB:2' 'fa.cpp' 'break'
commands
    silent
    pkct k
end

break_re '//GDB:move$' 'fa.cpp' 'break'
commands
    silent
    pkct k
end

#break_re 'b.moved = b.sk - thisk + 1;$' 'fa.cpp' 'break'
#commands
#    silent
#    pkct
#end


define reached_boundary
    bt 1
    info locals
end
break_re '// GDB$' 'fa.cpp' 'break' 'reached_boundary'


break_re 'kc->last_uqct = kc->uqct;' 'fa.cpp' 'break'

#########################################################################################
# mapping

b load_kc

#b map_fq_se


#leave:
r

