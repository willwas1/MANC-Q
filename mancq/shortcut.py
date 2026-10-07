"""Puts a MANC-Q shortcut on the Windows desktop:  python -m mancq shortcut

The shortcut starts MANC-Q with the same Python that ran this command (Anaconda, a normal Python, or a virtual
environment), without a black console window. For Anaconda it goes through Anaconda's own launcher (cwp.py),
the same way the Anaconda Navigator shortcuts do, so the environment is set up properly.
"""
import os
import subprocess
import sys


def _find_cwp():
    """Anaconda's cwp.py: in the base installation, which is sys.prefix or two levels up for a named env."""
    p = os.path.abspath(sys.prefix)
    for _ in range(3):
        c = os.path.join(p, "cwp.py")
        if os.path.exists(c):
            return c
        p = os.path.dirname(p)
    return None


def _ps(s):
    return "'" + str(s).replace("'", "''") + "'"


def make_shortcut(name="MANC-Q"):
    if os.name != "nt":
        print("Desktop shortcuts are only made on Windows. On other systems, start MANC-Q with:  python -m mancq")
        return 1
    exe_dir = os.path.dirname(os.path.abspath(sys.executable))
    pythonw = os.path.join(exe_dir, "pythonw.exe")
    if not os.path.exists(pythonw):
        pythonw = os.path.abspath(sys.executable)
    cwp = _find_cwp()
    if cwp:                                    # Anaconda: activate the environment first
        args = f'"{cwp}" "{os.path.abspath(sys.prefix)}" "{pythonw}" -m mancq gui'
        how = "Anaconda"
    else:
        args = "-m mancq gui"
        how = "this Python"
    workdir = os.path.join(os.path.expanduser("~"), "Documents")
    if not os.path.isdir(workdir):
        workdir = os.path.expanduser("~")
    ps = ("$d = [Environment]::GetFolderPath('Desktop'); "
          f"$p = Join-Path $d {_ps(name + '.lnk')}; "
          "$s = (New-Object -ComObject WScript.Shell).CreateShortcut($p); "
          f"$s.TargetPath = {_ps(pythonw)}; $s.Arguments = {_ps(args)}; "
          f"$s.WorkingDirectory = {_ps(workdir)}; $s.Description = 'MANC-Q NMR metabolite quantification'; "
          f"$s.IconLocation = {_ps(pythonw + ',0')}; $s.Save(); Write-Output $p")
    r = subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        print("Could not make the shortcut automatically. Make it by hand instead:\n"
              "  right-click the desktop > New > Shortcut, and for the location paste:\n"
              f'  "{pythonw}" {args}\n'
              "  then name it MANC-Q.")
        if r.stderr.strip():
            print("(Windows said: " + r.stderr.strip().splitlines()[-1] + ")")
        return 1
    print(f"Made {r.stdout.strip()}  (starts MANC-Q with {how}). Double-click it to open MANC-Q.")
    return 0


if __name__ == "__main__":
    sys.exit(make_shortcut())
