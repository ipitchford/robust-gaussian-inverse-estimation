# Replay

With Python 3.12 and the pinned requirements:

```sh
python code/replay.py --output replay-output.json
python -O code/replay.py --output replay-output-wrapper-optimized.json
```

The output path must be new. The runner hashes both verifier sources, executes them in fresh **normal-mode** subprocesses, and checks the finite checker's summary and the exact certificate's `PASS` status. The second command optimizes only the orchestration wrapper; its receipt records `verifier_subprocess_optimized: false`. It is not evidence that the assertion-based verifier subprocesses ran under `-O`.

The replay validates finite instances and the exact two-band arithmetic. It does not execute the general uniform proof, establish novelty, or bind archive/provider identities.
