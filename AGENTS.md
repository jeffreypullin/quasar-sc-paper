# Agent instructions

## Never run Nextflow

Do not invoke `nextflow` in this repository under any circumstances. The user runs
the pipeline themselves. This applies to every subcommand and flag, including ones
that look read-only:

- `nextflow run ...` — no
- `nextflow run ... -preview` — no
- `nextflow run ... -stub-run` — no
- `nextflow inspect`, `nextflow log`, `nextflow clean` — no

`-preview` and `-stub-run` are **not** safe. They do not execute tasks, but they do
append an entry to `.nextflow/history` with a **new session UUID**. Since
`nextflow run -resume` with no explicit session ID resumes the *most recent* run, a
single preview silently redirects the user's next `-resume` to an empty cache and
forces a multi-hour rerun from scratch.

Do not edit, move, or delete `.nextflow/`, `.nextflow.log*`, or the work directory
(`data/nextflow`, set in `conf/csd3.config`) either. Reading them is fine.

## Validating pipeline changes instead

When you change `main.nf`, `modules/*.nf`, or `lib/*.groovy`, verify by inspection
rather than execution:

- Re-read the channel definitions and check operator arities by hand. `combine(..., by: 0)`
  appends the second tuple's remaining elements, so the closure that follows must
  destructure exactly the right number of arguments.
- Check that every process added to a `modules/*.nf` file is also added to the
  matching `include { ... }` list at the top of `main.nf`.
- Create and edit all text files with LF line endings from the outset. Never write
  CRLF line endings. In particular, new `bin/` scripts must be LF-only and
  executable so their shebangs work without any normalization step.
- Test `bin/` scripts standalone on small synthetic inputs. This is encouraged and
  is the right way to check a script works.

Then hand the change back to the user to run.

