import subprocess
import gdb
import re


class InScope (gdb.Function):
    """Check if all the given variables or macros are in scope.
       Receives as argument a list of names separated by
       whitespace."""

    def __init__(self):
        super(InScope, self).__init__("in_scope")

    def invoke(self, var):
        vars = set(var.string().split())
        found = set()
        pc = gdb.get_selected_frame().get_pc()
        block = gdb.get_block_for_pc(pc)
        while block:
            for sym in block:
                if (sym.is_argument() or sym.is_constant() or
                        sym.is_function() or sym.is_variable()):
                    sym_name = sym.get_print_name()
                    if sym_name in vars:
                        found.add(sym_name)
            block = block.get_superblock()

        return vars == found

InScope()


# name_in_program = gdb.parse_and_eval('p->name').string()

def break_re_py2(e, f, b, allow_multi=None):
    with open(f, 'r') as o:
        a = [i for (i, l) in enumerate(o, 1) if re.search(e, l)]
        if len(a) > 1 and not allow_multi:
            return
        gdb.execute("printf \""+f+" ~ "+e+" => "+str(len(a))+" matches\n\"")
        for i in a:
            gdb.execute(b+" "+f+":"+str(int(i)))
            gdb.execute("break_commands $bpnum")


def break_re_py(m, f, b, allow_multi=None):
    grepcmd = "grep -n '"+m+"' '"+f+"'|cut -d ':' -f 1"
    lines = subprocess.check_output(grepcmd, shell=True).splitlines()
    if len(lines) > 1 and not allow_multi:
        return
    for n in lines:
        gdb.execute(b+" "+f+":"+str(int(n)))
