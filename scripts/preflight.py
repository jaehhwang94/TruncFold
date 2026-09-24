#!/usr/bin/env python
import importlib, shutil, subprocess, sys

def check_import(name):
    try: m=importlib.import_module(name); print(f'[OK] import {name} {getattr(m,"__version__","")}'); return True
    except Exception as e: print(f'[FAIL] import {name}: {e}'); return False

def main():
    ok=True
    for m in ['numpy','pandas','scipy','Bio','truncfold']: ok=check_import(m) and ok
    try: import freesasa; print('[OK] FreeSASA Python module')
    except Exception as e: print('[FAIL] freesasa:',e); ok=False
    exe=shutil.which('prodigy'); print('[OK] prodigy CLI '+exe if exe else '[FAIL] prodigy CLI missing'); ok=bool(exe) and ok
    if shutil.which('nvidia-smi'):
        p=subprocess.run(['nvidia-smi','--query-gpu=name,memory.total','--format=csv,noheader'],capture_output=True,text=True); print('[GPU]',p.stdout.strip())
    else: print('[WARN] nvidia-smi not found in this environment')
    sys.exit(0 if ok else 1)
if __name__=='__main__': main()
