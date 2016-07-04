
source gdb_py/frame_filter.py
#source /home/roel/.local/lib/python2.7/site-packages/voltron/entry.py
#voltron init
set disassembly-flavor intel

define hook-quit
    set confirm off
end

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

define pdna
    if $argc != 2
        help pdna
    else
        call print_dna($arg0, $arg1, KEY_WIDTH)
    end
end

document pdna
	Prints dna.
	Syntax: pdna <dna (or number)> <delimitor>
end


define pseq
    if $argc != 0
        help pseq
    else
        printf "[\t+%u]:\t", seq.p>>1
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
    if $argc != 3
        help break_re
    else
        python break_re_py($arg0, $arg1, $arg2)
    end
end

document break_re
	Insert specified break in file where match occurs [must be only one].
	Syntax: break_re <match> <file> <break|tbreak>
end

set $dbg=9

define run_until
    if kc->iter != 0
        set $dbg=9
    end
    c
end


# leave this: is for buffers.
#tb key_init.cpp:190
#break_re '_addtoseq(kc->s, seq.t); // kc->s_l grows here' 'key_init.cpp' 'tbreak'
#commands
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
#b fa_read

#b key_init.cpp:190
break_re 'get_kct(kc, seq, 1);' 'key_init.cpp' 'break'
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
        if seq.p != 0
            printf "N stretch of %u\n", corr
        else
            printf "N stretch at start of contig (not inserted)\n"
        end
    end
    c
end

break_re 'kc->ct += kc->kct_l;' 'key_init.cpp' 'tbreak'
commands
    silent
    if $dbg > 6
        print show_mantras(kc, 0, 0, kc->bnd)
    end
    run_until
end


#################################################################
# fa.cpp: print definitions

define dbg_whereat
    if $dbg > 4
        set listsize 1
        list
        set listsize 10
    end
end

define dbg_pseq
    if $dbg > $arg0
        call print_seq(&seq, KEY_WIDTH)
    end
end

define dbg_print
    if $dbg > $arg0
        printf $arg1
    end
end

define dbg_kct
    if $dbg > $arg0
        call print_kct(kc, bnd, b, k)
    end
end

define dbg_mantras
    if $dbg > $arg0
        print show_mantras(kc, b.obnd, b.obnd_l, bnd)
    end
end

#################################################################
# fa.cpp: debug print locations

define break_commands
    commands $arg0
        silent
        dbg_whereat
        run_until
    end
end
python break_re_py2('//P;', 'fa.cpp', 'break', 1)

define break_commands
    commands $arg0
        silent
        dbg_pseq 6
        #dbg_whereat
        run_until
    end
end
python break_re_py2('//O; .* occurance', 'fa.cpp', 'break', 1);

break_re '//~ uniq$' 'fa.cpp' 'break'
commands
    silent
    if *k
        if ~*k & DUP_BIT
            if $dbg > 5
                call print_posseq(b.s, *k, KEY_WIDTH)
            end
            dbg_print 5 "uniq----^^^\n"
        else
            if $dbg > 6
                call print_posseq(b.s, *k, KEY_WIDTH)
            end
        end
    end
    run_until
end

break_re '//~ also update header$' 'fa.cpp' 'break'
commands
    silent
    if $dbg > 0
        printf "stored k offset %u for hdr %u\nnext hdr\n", b.tgtk - kc->kct, bnd->ho - 1
    end
    if $dbg > 1
        printf "2bit sequence offset became %u:\t", b.s + kc->h[bnd->ho]->len - kc->s
        call print_dna(b.s[kc->h[bnd->ho]->len], '.', 4)
        printf "..\n"
    end
    run_until
end


define break_commands
    commands $arg0
        silent
        dbg_whereat
        dbg_mantras 7
        run_until
    end
end
python break_re_py2('//M;', 'fa.cpp', 'break', 1)

define break_commands
    commands $arg0
        silent
        dbg_whereat
        dbg_kct 8
        run_until
    end
end
python break_re_py2('//K;', 'fa.cpp', 'break', 1)

define break_commands
    commands $arg0
        silent
        dbg_whereat
        dbg_kct 8
        dbg_mantras 7
        run_until
    end
end
python break_re_py2('//B;', 'fa.cpp', 'break', 1)

define break_commands
    commands $arg0
        silent
        if $dbg > 6
            pseq
        end
        run_until
    end
end
python break_re_py2('//S;', 'fa.cpp', 'break', 1)


#TODO:
# branch analysis:
# store if/for/(do-)while branch taken info, print if taken first n times. (not taken is difficult)

#TODO:
# parse line and print all variables here. (i.e. before execution)


#########################################################################################
# mapping

b load_kc

#b map_fq_se


#leave:
r

