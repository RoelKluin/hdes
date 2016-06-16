
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

define pkct
    #bt
    if $argc != 0
        call print_kct(kc, b, $arg0)
    else
        call print_kct(kc, b, 0)
    end
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


define run_until
#    if kc->extension == 0 || kc->iter == 0
        c
#    else
#        bt 1
#    end
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

break_re 'kc->uqct += kc->kct_l;' 'key_init.cpp' 'tbreak'
commands
    silent
    print show_mantras(kc, kc->bnd->begin())
    run_until
end


#################################################################
# fa.cpp


define handle_non_uniques
    #bt 2
    #printf "prev:%u\tseq.p:%u\tpend:%u\text-1:%u\t", prev, seq.p, pend, kc->ext + 1
    pseq
    run_until
end

#b fa.cpp:157
break_re 'for (;;) {' 'fa.cpp' 'tbreak'
commands
#    silent
    pkct kc->kct + *contxt_idx
    print show_mantras(kc, b.it)
#    wa seq.dna
#    commands
#        #silent
#        pseq
#        run_until
#    end
    break_re 'for(;;) {' 'fa.cpp' 'break'
    commands
#        wa seq.dna
#        commands
#            #silent
#            pseq
#            run_until
#        end
#        #silent
        handle_non_uniques
    end
    handle_non_uniques
end

#b fa.cpp:251
#break_re '/NB(hdr_end_k(kc, h) >= b.tgtk);' 'fa.cpp' 'break'
#commands
#    #silent
#    printf "uniq at\t%u\n", *k >> 1
#    run_until
#end

#b fa.cpp:201
#
## if the assigned break was 1:
#
#command 1
#pdna seq 3 '\n'
#p seq.p + (s - kc->s)
#end

#b fa.cpp:216
#break next_mantra
#break_re '// update assembly' 'fa.cpp' 'break'
#commands
#    #silent
#    printf "boundary update:\n"
#    run_until
#end

break_re '//GDB:1$' 'fa.cpp' 'break'
commands
#   silent
    pdna seq.dna ','
    run_until
end

break_re '//GDB:moved' 'fa.cpp' 'break'
commands
#   silent
    printf "\nThese were moved to kct end.\n"
    pkct thisk
    run_until
end

break_re '//GDB:2$' 'fa.cpp' 'break'
commands
#    silent
    pkct thisk
    run_until
end

break_re '//GDB:move$' 'fa.cpp' 'break'
commands
#    silent
    printf "^^^---moved up\n"
    pkct k
    run_until
end

#break shrink_mantra
break_re '//GDB:mantra1$' 'fa.cpp' 'break'
commands
#    silent
    print show_mantras(kc, b.it)
    run_until
end

#break reached_boundary
#break_re '//GDB:mantra2$' 'fa.cpp' 'break'
#commands
#    #silent
#    print show_mantras(kc, b.it)
#    run_until
#end

break_re '//GDB:next mantra$' 'fa.cpp' 'break'
commands
#    silent
    printf "next mantra\n"
    print show_mantras(kc, b.it)
    run_until
end


#break reached_boundary
#break_re 'b.moved = b.tgtk - thisk + 1;$' 'fa.cpp' 'break'
#commands
#    #silent
#    pkct
#    run_until
#end

#break next_mantra if b.prev == 0

#define reached_boundary
#    bt 1
#    info locals
#end
#break_re '// GDB$' 'fa.cpp' 'break' 'reached_boundary'

###########################################################################
# ext_uq_iter()

break_re '//GDB:UQ1$' 'fa.cpp' 'break'
commands
#    silent
    call print_posseq(b.s, *k, KEY_WIDTH)
    if ~*k & DUP_BIT
        printf "uniq----^^^\n"
    end
    run_until
end

break_re '//GDB:UQ2$' 'fa.cpp' 'break'
commands
#    silent
    call print_posseq(b.s, *k, KEY_WIDTH)
    if ~*k & DUP_BIT
        printf "uniQ----^^^\n"
    end
    run_until
end

break_re '// check whether last uniq was adjoining end' 'fa.cpp' 'break'
commands
    printf "new hdr\n"
    pkct k
    print show_mantras(kc, b.it)
    run_until
end

break_re '// also update header$' 'fa.cpp' 'break'
commands
#    silent
    printf "stored offset %u for hdr %u\nnext hdr\n", b.tgtk - kc->kct, (*b.it).ho
    if (*b.it).ho != kc->h_l - 1
        printf "2bit sequence offset became %u:\t", b.s + kc->h[(*b.it).ho]->len - kc->s
        call print_dna(b.s[kc->h[(*b.it).ho]->len], '.', 4)
        printf "..\n"
    else
        printf "(looping)\n"
    end
    run_until
end

break_re 'kc->uqct = k - b.tgtk;' 'fa.cpp' 'break'
commands
#    silent
    pkct k
    print show_mantras(kc, b.it)
    run_until
end

break_re 'kc->kct_l = skctl;' 'fa.cpp' 'break'
commands
#    silent
    #b ext_uq_iter
    pkct kc->kct + kc->kct_l
    print show_mantras(kc, b.it)
    run_until
end

break_re '//GDB:BUG$' 'fa.cpp' 'break'
commands
#    silent
    #b ext_uq_iter
    p/x k[-1]
    run_until
end



#########################################################################################
# mapping

b load_kc

#b map_fq_se


#leave:
r

