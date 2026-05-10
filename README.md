# pgen-samplebind

A focused tool that does the one thing `plink2 --pmerge` doesn't yet do: bind two or more PFILE/BFILE/EIGENSTRAT datasets that share variants but contain different samples.

Targeted at the ancient-DNA / population-genetics community working with the 1240k SNP capture panel, where this operation is a daily workflow blocker.

## Status

Early development (v0.1.0.dev0). API surface defined, implementation in progress. See `cs-wiki/projects/pgen-samplebind.md` (HLD) and `cs-wiki/projects/pgen-samplebind-lld.md` (LLD) for the full design.

## Quick start

```bash
pip install pgen-samplebind
pgen-samplebind --help
```

## Subcommands

- `pgen-samplebind merge`    — bind inputs into one output PFILE
- `pgen-samplebind validate` — check alignment, no output written
- `pgen-samplebind hash`     — emit canonical variant-set hash
- `pgen-samplebind inspect`  — structured summary of one input

## License

MIT
