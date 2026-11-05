# 17. Command-line tool

The same solver is available headless, for scripting and batch runs:

```sh
cargo run --release -p stellar-cli -- [options]
```

| Option | Meaning |
|---|---|
| `--mass <M>` | stellar mass in M☉ (default 1) |
| `--x <X>` | hydrogen mass fraction (default 0.70) |
| `--z <Z>` | metallicity (default 0.02) |
| `--alpha <α>` | mixing-length parameter (default 1.8) |
| `--f-cno <f>` | CNO fraction of $`Z`$ (default 0.7) |
| `--points <N>` | mesh points (default 800) |
| `-o, --output <file>` | write the profile CSV |
| `--json <file>` | write the summary JSON |
| `--sweep <M1,M2,…>` | solve a ZAMS sequence by continuation and print a table |
| `--evolve` | evolve the star from the ZAMS and print the track |
| `--until <tams\|age\|max>` | stop condition: end of main sequence (default), age, or as far as possible |
| `--max-age <Gyr>` | target age with `--until age`; otherwise a safety limit |
| `--track <file>` | write the evolutionary track (one row per model) as CSV |
| `-v, --verbose` | print every Newton iteration (every model with `--evolve`) |
| `-h, --help` / `-V, --version` | print help / version |

## Examples

```sh
stellar-cli --mass 1 -o sun.csv --json sun.json
stellar-cli --sweep 0.5,1,2,5,10,20
stellar-cli --mass 3 --z 0.008 -v
stellar-cli --mass 1 --evolve --until age --max-age 4.57
stellar-cli --mass 5 --evolve --until max --track track.csv -o final_profile.csv
```

The CLI and the app both use the `stellar-core` library, so identical parameters give the same model to
within the solver tolerance. The app warm-starts from the previous model; the CLI starts cold unless
`--sweep` is used.

---
← Previous: [Constants & units](16-constants.md) · [Contents](README.md) · Next: [References & software](18-references.md) →
