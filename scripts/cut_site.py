# pylint: disable=line-too-long
# pylint: disable=invalid-name

"""
Automatic detection of the cut site from an Alleles_frequency_table.

Replaces the hand-maintained 'cut_position' key in parameters.txt. That key was
only ever correct when (a) every table in a run shared one cut position and
(b) CRISPResso2's plot window was not clipped by the amplicon edge. Neither
holds for pooled panels: the Pallaseni 2024 QC-pass panel has 5,583 targets at
TWO different cut positions (39 and 53), so no single constant can be right.

WHAT "CUT POSITION" MEANS HERE
------------------------------
Two different coordinate systems are in play and the old code conflated them:

  amplicon coordinate  index into the ungapped reference amplicon
                       (Original_Sequence, PAM, 3rd_nt, 4th_nt, ...)
  alignment coordinate index into the gapped Aligned_Sequence /
                       Reference_Sequence of ONE row

They coincide only when that row's reference carries no gaps 5' of the cut and
starts at amplicon position 0. detect_cut_position() returns the AMPLICON
coordinate. project_cut_to_alignment() converts it to the alignment coordinate
for a given row.

WHY A ROW'S REFERENCE CAN DIFFER FROM THE AMPLICON
--------------------------------------------------
Substitutions never touch Reference_Sequence -- it is the amplicon with '-'
padding and nothing else. So for every row:

    row.Reference_Sequence.replace('-', '')  is a CONTIGUOUS SLICE of the amplicon

It is the whole amplicon for most rows, and a shorter slice when the read
carries an insertion (which consumes width in a fixed-width plot window) or
runs past the window edge. amplicon_reference() recovers the amplicon as the
longest such slice and records each row's offset into it.
"""

IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}

COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(seq):
    """Reverse complement of a DNA string."""
    return seq.translate(COMPLEMENT)[::-1]


def pam_matches(observed, pattern):
    """IUPAC-aware PAM comparison."""
    if len(observed) != len(pattern):
        return False
    return all(o in IUPAC.get(p, '') for o, p in zip(observed, pattern))


# cut_offset is measured from the 5' end of the protospacer, on the strand the
# protospacer is written on. spacer_len=None means "use len(gRNA_seq)".
#
# SpCas9   : blunt cut 3 bp 5' of a 3' NGG PAM  -> offset = len(spacer) - 3
# LbCas12a : staggered cut, offset 23 from a 5' TTTV PAM. This reproduces the
#            convention the original CRISPRlysis code used (gRNA_pos + 23 with
#            a 'TTT' check immediately 5' of the spacer) so existing Cas12a
#            results do not move.
NUCLEASES = {
    'SpCas9':   {'pam': 'NGG',  'pam_side': "3'", 'cut_offset': lambda n: n - 3},
    'SaCas9':   {'pam': 'NNGRRT', 'pam_side': "3'", 'cut_offset': lambda n: n - 3},
    'LbCas12a': {'pam': 'TTTV', 'pam_side': "5'", 'cut_offset': lambda n: 23},
    'AsCas12a': {'pam': 'TTTV', 'pam_side': "5'", 'cut_offset': lambda n: 23},
}


class CutSiteError(ValueError):
    """Raised when the cut site cannot be determined from a table."""


# ---------------------------------------------------------------------------
# 1. recover the amplicon and each row's offset into it
# ---------------------------------------------------------------------------

def amplicon_reference(reference_sequences):
    """Return (amplicon, {gapped_reference: offset}).

    `reference_sequences` is any iterable of gapped Reference_Sequence strings.
    Raises CutSiteError if some row's reference is not a unique contiguous
    slice of the longest one -- that would mean the table mixes amplicons, or
    that a slice is so short it matches the amplicon in more than one place.
    """
    gapped = list(dict.fromkeys(reference_sequences))       # unique, order kept
    degapped = {g: g.replace('-', '').upper() for g in gapped}
    if not degapped:
        raise CutSiteError('empty allele table: no Reference_Sequence rows')

    # PATCH 2026-11-13 (preprint rerun): CRISPResso's quantification window is not
    # identical for every allele - rows whose alignment reaches further 5' or 3' carry a
    # reference that EXTENDS the longest one instead of being a slice of it (seen in
    # Kurgan 2021: 6 of 496 tables, e.g. APOC3-1_129, where 22 of 56 variants start 2 bp
    # earlier). Assemble the amplicon from such overhanging variants before testing
    # containment, so only genuinely foreign sequences raise.
    amplicon = max(degapped.values(), key=len)
    MIN_OVERLAP = 20
    changed = True
    while changed:
        changed = False
        for d in sorted(degapped.values(), key=len, reverse=True):
            if d in amplicon or len(d) < MIN_OVERLAP:
                continue
            for k in range(min(len(d), len(amplicon)), MIN_OVERLAP - 1, -1):
                if d[-k:] == amplicon[:k]:                  # variant hangs over the 5' end
                    amplicon = d[:-k] + amplicon; changed = True; break
                if d[:k] == amplicon[-k:]:                  # variant hangs over the 3' end
                    amplicon = amplicon + d[k:]; changed = True; break
            if changed:
                break
    offsets = {}
    ambiguous, foreign = [], []
    for g, d in degapped.items():
        hits = amplicon.count(d)
        if hits == 0:
            foreign.append(d)
        elif hits > 1:
            ambiguous.append(d)
        else:
            offsets[g] = amplicon.find(d)
    if foreign:
        raise CutSiteError(
            f'{len(foreign)} reference variant(s) are not a slice of the '
            f'{len(amplicon)} bp amplicon (shortest {min(map(len, foreign))} bp). '
            'The table may mix amplicons.')
    if ambiguous:
        raise CutSiteError(
            f'{len(ambiguous)} reference variant(s) match the amplicon in more '
            f'than one place (shortest {min(map(len, ambiguous))} bp), so their '
            'position is undefined.')
    return amplicon, offsets


# ---------------------------------------------------------------------------
# 2. locate the cut in amplicon coordinates
# ---------------------------------------------------------------------------

def detect_cut_position(reference_sequences, nuclease, gRNA_seq, fallback=None):
    """Cut position in AMPLICON coordinates (index of the first base 3' of the cut).

    Returns (cut_position, info) where info is a dict describing how it was
    found -- write it into the report so every run is auditable.

    `fallback` (normally parameters_dic.get('cut_position')) is used only when
    the nuclease is not in NUCLEASES, i.e. when there is no PAM rule to apply.
    """
    amplicon, _ = amplicon_reference(reference_sequences)

    if nuclease not in NUCLEASES:
        if fallback is None:
            raise CutSiteError(
                f"nuclease {nuclease!r} is not one of {sorted(NUCLEASES)} and no "
                "'cut_position' fallback is set in parameters.txt")
        return int(fallback), {
            'source': 'parameters.txt fallback',
            'nuclease': nuclease, 'amplicon_len': len(amplicon),
            'strand': None, 'protospacer_start': None, 'pam': None,
        }

    spec = NUCLEASES[nuclease]
    guide = str(gRNA_seq).upper()
    if not guide or set(guide) - set('ACGTN'):
        raise CutSiteError(f'gRNA sequence {gRNA_seq!r} is not a plain DNA string')
    n = len(guide)

    hits = []
    for strand, probe in (('+', guide), ('-', revcomp(guide))):
        start = amplicon.find(probe)
        while start != -1:
            hits.append((strand, start))
            start = amplicon.find(probe, start + 1)
    if not hits:
        raise CutSiteError(
            f'gRNA {guide} not found on either strand of the amplicon {amplicon}')
    if len(hits) > 1:
        raise CutSiteError(
            f'gRNA {guide} matches the amplicon {len(hits)} times at {hits}; '
            'the cut site is ambiguous')

    strand, pos = hits[0]
    pam_len = len(spec['pam'])
    if strand == '+':
        lo, hi = ((pos + n, pos + n + pam_len) if spec['pam_side'] == "3'"
                  else (pos - pam_len, pos))
        cut = pos + spec['cut_offset'](n)
    else:
        lo, hi = ((pos - pam_len, pos) if spec['pam_side'] == "3'"
                  else (pos + n, pos + n + pam_len))
        cut = pos + n - spec['cut_offset'](n)

    if lo < 0 or hi > len(amplicon):
        raise CutSiteError(
            f'the {nuclease} PAM would fall outside the amplicon (needs '
            f'[{lo}:{hi}] of {len(amplicon)} bp). Widen the CRISPResso2 window '
            'or supply a longer amplicon.')
    observed = amplicon[lo:hi]
    if strand == '-':
        observed = revcomp(observed)
    if not pam_matches(observed, spec['pam']):
        raise CutSiteError(
            f'observed PAM {observed} does not match {spec["pam"]} for {nuclease} '
            f'(protospacer at {pos} on the {strand} strand). Check that the file '
            'name really ends in the protospacer sequence.')

    if not 0 <= cut <= len(amplicon):
        raise CutSiteError(f'computed cut {cut} lies outside the amplicon')

    return cut, {
        'source': f'{nuclease} PAM rule',
        'nuclease': nuclease, 'amplicon_len': len(amplicon),
        'strand': strand, 'protospacer_start': pos, 'pam': observed,
    }


# ---------------------------------------------------------------------------
# 3. amplicon coordinate -> alignment column, per row
# ---------------------------------------------------------------------------

def project_cut_to_alignment(ref_aln, offset, cut_position, insertion_side='left'):
    """Alignment-column index of `cut_position` for one row, or None if that
    row's reference does not cover the cut.

    Walks the gapped reference counting real bases, starting the counter at the
    row's offset into the amplicon. Substitutions are invisible here by
    construction -- they live in Aligned_Sequence only.

    insertion_side decides where the cut goes when an insertion sits exactly on
    it: 'left' puts the cut before the inserted bases (they read as downstream
    of the cut), 'right' puts it after them.
    """
    seen = offset
    for col, base in enumerate(ref_aln):
        if base == '-':
            continue
        if seen == cut_position:
            if insertion_side == 'left':
                while col > 0 and ref_aln[col - 1] == '-':
                    col -= 1
            return col
        seen += 1
    if seen != cut_position:
        return None
    # The cut sits immediately after the last reference base this row covers.
    # Any trailing gap columns are an insertion at the cut, so honour
    # insertion_side here too rather than always returning the far end.
    col = len(ref_aln)
    if insertion_side == 'left':
        while col > 0 and ref_aln[col - 1] == '-':
            col -= 1
    return col


def cut_position_per_row(reference_sequences, cut_position, insertion_side='left'):
    """{gapped_reference: alignment column or None} for a whole table.

    Computed once per distinct gap pattern rather than once per row -- allele
    tables have far more rows than distinct references.
    """
    _, offsets = amplicon_reference(reference_sequences)
    return {g: project_cut_to_alignment(g, off, cut_position, insertion_side)
            for g, off in offsets.items()}
