from pathlib import Path
import subprocess, re

def run_prodigy(pdb_path,prodigy_exe='prodigy'):
    cmd=[prodigy_exe,str(Path(pdb_path)),'--selection','A','B','-q']
    p=subprocess.run(cmd,check=True,capture_output=True,text=True)
    text=(p.stdout or '').strip()
    # quiet mode: either '<name> -9.373' or just '-9.373'
    nums=re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?',text)
    if not nums: raise RuntimeError(f"Could not parse PRODIGY output: {text!r}")
    return float(nums[-1])
