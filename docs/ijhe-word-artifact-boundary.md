# IJHE Word-artifact boundary

## Checked files

The two Word files supplied from the separate `SLABx_LH2` workspace were
checked on 8 October 2026:

- `01_Manuscript.docx`
- `02_Supplementary_Information.docx`

They are a populated SLABx–LH2/PSEP submission package. The manuscript title,
abstract, 210-record FFI receptor field, SLABx version lineage, and PSEP
reproducibility capsule are all specific to that package. They are **not** the
DEGALI IJHE manuscript source in this repository, whose evidence lanes are the
PRESLHY E3.5 no-fit baseline, the six-arc FFI/DNV screen, the Hecht–Panda
boundary, and the separate transient 3-D numerical verification.

Do not upload the supplied Word files as the DEGALI IJHE submission and do not
merge their numerical claims into `MANUSCRIPT_LH2_VALIDATION.md` or
`IJHE_SUPPLEMENTARY_INFORMATION.md`. The two evidence packages may be cited as
separate development lines, but they must retain their own source, observation
operator, and scope boundaries.

## Current DEGALI submission source

The canonical IJHE draft is maintained as the following version-controlled
source set:

- `MANUSCRIPT_LH2_VALIDATION.md`
- `IJHE_SUPPLEMENTARY_INFORMATION.md`
- `IJHE_HIGHLIGHTS.md`
- `outputs/ijhe-highlights-2026-10-09.docx` (standalone portal upload)
- `IJHE_COVER_LETTER_DRAFT.md`
- `IJHE_GRAPHICAL_ABSTRACT.svg`
- `docs/ijhe-figure-table-register.md`

Before upload, the corresponding author must provide the author list,
contributions, competing-interest/funding wording, and approve the final
Word/PDF conversion. The supplied SLABx Word declarations cannot be copied
into the DEGALI manuscript without author confirmation.

An author-neutral draft can now be generated reproducibly with
`tools/build_ijhe_docx.py`. The local draft
`outputs/ijhe-manuscript-draft-2026-10-09.docx` was exported through the
installed Microsoft Word renderer to a 12-page PDF and visually inspected
page by page. It embeds the six main figures next to their captions for the
single-file review copy; canonical SVG/PDF figure artifacts remain separate
upload files. It deliberately retains author placeholders and is a
layout-verified draft, not the final submission file.

## Visual QA status

The bundled LibreOffice render script could not rasterize DOCX files on this
host because `soffice.exe` is unavailable. The supplied Word package was not
edited or overwritten. Microsoft Word was used only to render the generated
author-neutral draft; the final upload still requires author fields, figure
permissions and a final Word/PDF approval.
