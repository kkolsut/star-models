# Embedded physics data

Tables compiled into `stellar-core` with `include_str!` (no file access at run time, so the web
build works). Read by `src/physics/opacity_table.rs`. Plan task 7.1; `docs/physics/03-opacity.md`.

| Directory | Content | Authors and reference | Retrieved from |
|---|---|---|---|
| `opal/GN93hz` | OPAL Rosseland mean opacities, type 1, 126 (X, Z) tables, Grevesse & Noels (1993) mixture, file dated 10/03/95 | C. A. Iglesias & F. J. Rogers, *Updated OPAL opacities*, ApJ 464, 943 (1996) | mirror of A. I. Boothroyd, CITA, `https://www.cita.utoronto.ca/~boothroy/data/kappa/GN93hz.gz` (the LLNL site `opalopacity.llnl.gov` was not reachable from the build environment); md5 equal to the mirror's `OPALmd5sum.txt` |
| `ferguson2005_gn93/*.tron` | Low-temperature Rosseland mean opacities with molecules and grains, 155 (X, Z) files, Grevesse & Noels (1993) mixture (set `f05.g93`) | J. W. Ferguson, D. R. Alexander, F. Allard, T. Barman, J. G. Bodnarik, P. H. Hauschildt, A. Heffner-Wong & A. Tamanai, *Low-temperature opacities*, ApJ 623, 585 (2005) | Wichita State University, `http://webs.wichita.edu/physics/opacity/Data/f05.g93.tar.gz`, via the Internet Archive copy of 2016-07-29 (the live site blocks automated access); all 155 md5 sums equal those of the CITA mirror (`MOLf05_g93_md5sum.txt`) |
| `potekhin/condtab21wd.dat` | Electron thermal conductivity of fully ionised, non-magnetic electron–ion plasmas, ion charge 1–60, log ρ −6…9.75, log T 3…9 ("weakly damped" version of the 2021 update) | A. Y. Potekhin; S. Cassisi, A. Y. Potekhin, A. Pietrinferni, M. Catelan & M. Salaris, ApJ 661, 1094 (2007); S. Cassisi, A. Y. Potekhin, M. Salaris & A. Pietrinferni, A&A 654, A149 (2021); review: A. Y. Potekhin, J. A. Pons & D. Page, Space Sci. Rev. 191, 239 (2015) | `http://www.ioffe.ru/astro/conduct/condtab21wd.dat` (page last updated 2021-05-18), retrieved 2026-11-17 |

Checksums of the files as embedded:

| File | md5 |
|---|---|
| `opal/GN93hz` | `32c29813fda00e37fc2411d3cd6cb1f6` |
| `ferguson2005_gn93/*.tron` (concatenated in `ls` order) | `13d7dad6c4cff4c01cfd595f05870eae` |
| `potekhin/condtab21wd.dat` | `e0ed52d6d26f352414f7cb43d49103cb` |

The files are unchanged copies of the originals.

## Licence status (decision of task 7.1)

None of the three sources states a licence. The OPAL and Ferguson pages offer the tables for
download for use in stellar modelling; the Potekhin page asks users to cite Potekhin, Pons & Page
(2015). The tables are copied unchanged, as is common practice in stellar-evolution codes that
redistribute them (MESA, Boothroyd's OPAL routines, ESTER).

**Decided 2026-11-17:** the tables are embedded, with the full attribution above, in the About page
(§7 Microphysics and §18 References) and in the documentation. The repository's MIT licence does
not cover these files: they remain the work of their authors, and anyone redistributing the
binary should keep this attribution. If an author objects, the run-time "load tables" path is the
fallback (it would not work for the web build); written permission from the three groups has not
been requested yet and is listed as an open item in `docs/execution-plan.md`.
