import subprocess
import gdb


def break_re_py(m, f, b, cmd=None):
    grepcmd = "grep -n '"+m+"' '"+f+"'|cut -d ':' -f 1"
    for n in subprocess.check_output(grepcmd, shell=True).splitlines():
        gdb.execute(b+" "+f+":"+str(int(n)))
