# Release notes

One file per release, named for the tag: `release-notes/1.0.586.md`.

The Release workflow publishes the matching file as the body of the GitHub Release page when a tag
is pushed. Commit the file to `master` **before** tagging — the release job checks out the tag, so a
file added afterwards is not seen, and the page will already have been created from generated notes.

If no file exists for the version, the page is still created so that the Releases sidebar and the
`github/v/release` badge move — both read `releases/latest`, which is the Release rather than the
tag, so without a page they keep advertising the previous version indefinitely. That fallback body
is a list of pull request titles under a placeholder line saying so. It is a stopgap, not the
deliverable: edit the page afterwards, or better, write the file before tagging.

## What the file should say

Every line answers *so what?* A pull request title says what its author did; notes say what it means
for someone who already uses mzLib and needs to know whether to upgrade.

> `Pep model stability (#2642)`

against the same change written for a reader:

> **Match-between-runs was not reproducible (#1155)** — the PEP model's training rows were appended
> per thread under a lock, so a peptide whose only evidence was a transferred peak could report its
> real intensity in one run and zero in the next. If you have been averaging repeated runs to smooth
> this out, you can stop.

`1.0.586.md` is a worked example. What it does that a generated list cannot:

- **Leads with what breaks.** Anything that changes answers existing callers already have goes above
  the list, before the reader has to find it halfway down.
- **Orders by impact, not by pull request number.** Documentation-only items go last, labelled as
  such. Nobody minds a boring item marked boring; they mind hunting for the important one.
- **Names what does *not* change.** "Files predating `MSLEVEL` all carry `PEPMASS`, so they read
  exactly as before" stops a reader over-generalising a behavior change into a panic.
- **Gives the numbers.** Counts, sizes, versions, timings.
- **Points at the downstream fix** where there is one, so a breaking change arrives with its remedy.

Reading the diffs is what makes this possible, and it is the part that cannot be generated. A
release where the notes were skipped is a release where nobody found out that variant-database
q-values had been computed against a short decoy set.
