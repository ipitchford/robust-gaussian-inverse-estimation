#!/usr/bin/env python3
"""Run finite checks and the exact gap certificate with source-bound receipts."""
from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path
import json, subprocess, sys, time

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [ROOT / "code/inverse_risk_checks.py", ROOT / "code/analytic_isotropic_gap.py"]

def main():
    parser = ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output)
    if output.exists():
        raise SystemExit(f"refusing to overwrite {output}")
    started = time.time()
    runs = []
    for script in SCRIPTS:
        proc = subprocess.run([sys.executable, str(script)], cwd=ROOT, text=True, capture_output=True)
        if proc.returncode:
            raise SystemExit(proc.stderr or proc.stdout)
        payload = json.loads(proc.stdout)
        if payload.get("status") != "PASS":
            raise SystemExit(f"{script.name}: non-PASS status")
        runs.append({"script":script.relative_to(ROOT).as_posix(),"sha256":sha256(script.read_bytes()).hexdigest(),"result":payload})
    receipt = {"status":"PASS","python":sys.version,"optimized":not __debug__,"elapsed_seconds":time.time()-started,"runs":runs}
    output.write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"status":"PASS","output":str(output),"runs":len(runs)}, indent=2))

if __name__ == "__main__":
    main()
