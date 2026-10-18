# Retained replay receipt

Final normal-wrapper and optimized-wrapper receipts are stored in `certificates/replay-normal.json` and `certificates/replay-optimized.json`. Each binds the executed source hashes and reports both components as `PASS`. In both receipts the verifier subprocesses ran in normal mode; only the second orchestration wrapper used `-O`.

These are retained local runs, not independent reproduction or formal verification. Final archive and public-provider binding are separate publication gates.
