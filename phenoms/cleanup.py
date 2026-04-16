"""
Optional preprocessing: renumber a mobile PDB trajectory to match a reference PDB’s
residue numbers via sequence alignment (handles indels, e.g. Δ747–749 vs WT).

PHENOMS bond labels use PDB residue sequence numbers. Pick which structure defines the
target scheme (e.g. renumber WT to mut, or mut to WT). With ``fill_unmapped_mobile=True``,
mobile residues that align to gaps in the reference (longer WT vs deletion in mut, or
N/C overhangs) get consecutive numbers so every mapped standard residue can be rewritten.

No extra dependencies: global alignment is Needleman–Wunsch (pure Python).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple, Union

PathLike = Union[str, Path]

# Standard amino acid 3-letter -> 1-letter (common MD ff names).
THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
    "ASX": "B",
    "GLX": "Z",
    "XAA": "X",
    "UNK": "X",
}


@dataclass(frozen=True)
class ResidueKey:
    """Identifies one residue in a PDB (chain + sequence + insertion code)."""

    chain: str
    resseq: int
    icode: str = " "

    def __str__(self) -> str:
        ic = self.icode if self.icode and self.icode != " " else ""
        return f"{self.chain}:{self.resseq}{ic}"


@dataclass
class RenumberReport:
    """Summary of an align-and-renumber operation."""

    n_reference_residues: int
    n_mobile_residues: int
    n_mapped: int
    n_mismatch_at_aligned_positions: int
    n_mobile_only_gaps_in_reference: int  # mobile residues aligned to gap in ref (before fill)
    n_filled_unmapped_mobile: int  # assigned numbers for those (overhang + ref-gap columns)
    aligned_reference_sequence: str
    aligned_mobile_sequence: str


def _parse_atom_record(line: str) -> Optional[Tuple[ResidueKey, str]]:
    """
    Parse ATOM/HETATM line for residue identity and 1-letter code.
    Returns None if not a standard amino acid we map.
    """
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    if len(line) < 27:
        return None
    resname = line[17:20].strip().upper()
    chain = line[21] if len(line) > 21 else " "
    icode = line[26] if len(line) > 26 else " "
    try:
        resseq = int(line[22:26])
    except ValueError:
        return None
    letter = THREE_TO_ONE.get(resname)
    if letter is None:
        return None
    key = ResidueKey(chain=chain, resseq=resseq, icode=icode if icode != "" else " ")
    return key, letter


def _format_resseq_field(new_num: int) -> str:
    """PDB columns 23–26: right-justified residue sequence number, width 4."""
    s = str(int(new_num))
    if len(s) > 4:
        raise ValueError(f"Residue number {new_num} does not fit in PDB columns 23–26")
    return s.rjust(4)


def ordered_residue_run(
    pdb_path: PathLike,
    *,
    chain_id: Optional[str] = None,
) -> Tuple[List[ResidueKey], str]:
    """
    First occurrence order of protein residues in the PDB (file order), optional chain filter.

    Returns
    -------
    keys : list of ResidueKey
    sequence : str
        One-letter amino acid sequence in the same order.
    """
    path = Path(pdb_path)
    seen: Set[ResidueKey] = set()
    keys: List[ResidueKey] = []
    letters: List[str] = []

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parsed = _parse_atom_record(line)
            if parsed is None:
                continue
            key, letter = parsed
            if chain_id is not None and key.chain != chain_id:
                continue
            if key in seen:
                continue
            seen.add(key)
            keys.append(key)
            letters.append(letter)

    return keys, "".join(letters)


def _unique_chains_in_pdb(pdb_path: PathLike) -> List[str]:
    chains: set[str] = set()
    with Path(pdb_path).open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                chains.add(line[21])
    out = sorted(chains)
    return out


def _default_chain_if_needed(pdb_path: PathLike, chain_id: Optional[str]) -> str:
    if chain_id is not None:
        return chain_id
    chains = [c for c in _unique_chains_in_pdb(pdb_path) if c.strip()]
    if not chains:
        chains = _unique_chains_in_pdb(pdb_path)
    if len(chains) == 1:
        return chains[0]
    if len(chains) == 0:
        raise ValueError(f"No ATOM/HETATM records found: {pdb_path}")
    raise ValueError(
        f"Multiple chains {chains!r} in {pdb_path}; pass chain_id= explicitly."
    )


def needleman_wunsch(
    s1: str,
    s2: str,
    *,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
) -> Tuple[str, str]:
    """
    Global alignment of two strings. Returns (aligned_s1, aligned_s2) with '-' gaps.
    """
    n, m = len(s1), len(s2)
    dp: List[List[int]] = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_diag = dp[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            dp[i][j] = max(
                score_diag,
                dp[i - 1][j] + gap,
                dp[i][j - 1] + gap,
            )
    # Traceback
    i, j = n, m
    a1: List[str] = []
    a2: List[str] = []
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score_diag = dp[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            if dp[i][j] == score_diag:
                a1.append(s1[i - 1])
                a2.append(s2[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i - 1][j] + gap:
            a1.append(s1[i - 1])
            a2.append("-")
            i -= 1
            continue
        a1.append("-")
        a2.append(s2[j - 1])
        j -= 1
    return "".join(reversed(a1)), "".join(reversed(a2))


def _assign_filled_resseqs(
    n: int,
    left_bound: Optional[int],
    right_bound: Optional[int],
) -> List[int]:
    """
    Assign ``n`` consecutive PDB residue numbers for mobile-only alignment columns.

    - N-terminal overhang (no left bound): ``right_bound - n .. right_bound - 1`` (min 1).
    - C-terminal overhang (no right bound): ``left_bound + 1 .. left_bound + n``.
    - Internal (WT extra residues vs mut gap, e.g. deletion): prefer ``left+1 .. left+n``
      while staying below ``right_bound`` when possible; otherwise continue past ``left``.
    """
    if n <= 0:
        return []
    if right_bound is not None and left_bound is None:
        start = max(1, int(right_bound) - n)
        return list(range(start, start + n))
    if right_bound is None and left_bound is not None:
        L = int(left_bound)
        return list(range(L + 1, L + 1 + n))
    if left_bound is not None and right_bound is not None:
        L, R = int(left_bound), int(right_bound)
        span = R - L - 1
        if span >= n:
            return list(range(L + 1, L + 1 + n))
        # Not enough integers strictly between L and R: number consecutively from L+1
        # (covers classic deletion where mut jumps L -> R and WT has residues in between).
        return list(range(L + 1, L + 1 + n))
    # No anchors (should not happen for n > 0 in valid alignments)
    return list(range(1, n + 1))


def build_mobile_to_reference_map(
    ref_keys: Sequence[ResidueKey],
    ref_seq: str,
    mob_keys: Sequence[ResidueKey],
    mob_seq: str,
    *,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
    fill_unmapped_mobile: bool = False,
) -> Tuple[Dict[ResidueKey, int], RenumberReport]:
    """
    Align mobile sequence to reference; map each mobile residue key -> reference resseq.

    Matched columns (both non-gap) assign mobile residue the reference residue number
    at that column.

    If ``fill_unmapped_mobile`` is False, mobile residues aligned to a gap in the reference
    are left unmapped.

    If True, those residues get consecutive numbers: N/C overhangs extend from the first/last
    mapped reference number; internal gaps (e.g. WT residues where mut has a deletion) are
    numbered ``L+1, L+2, ...`` from the last mapped reference residue ``L`` before the gap
    (and fit below the next reference ``R`` when integers are available).
    """
    if len(ref_keys) != len(ref_seq) or len(mob_keys) != len(mob_seq):
        raise ValueError("keys and sequence lengths must match")

    aln_ref, aln_mob = needleman_wunsch(ref_seq, mob_seq, match=match, mismatch=mismatch, gap=gap)

    mapping: Dict[ResidueKey, int] = {}
    i_ref = i_mob = 0
    n_mismatch = 0
    n_mob_only = 0

    for a, b in zip(aln_ref, aln_mob):
        if a != "-" and b != "-":
            rk = ref_keys[i_ref]
            mk = mob_keys[i_mob]
            mapping[mk] = rk.resseq
            if a != b:
                n_mismatch += 1
            i_ref += 1
            i_mob += 1
        elif a != "-" and b == "-":
            i_ref += 1
        else:
            # gap in reference, residue in mobile
            n_mob_only += 1
            i_mob += 1

    n_filled = 0
    if fill_unmapped_mobile and n_mob_only > 0:
        mapping, n_filled = _fill_mobile_gaps_in_map(
            mapping,
            aln_ref,
            aln_mob,
            ref_keys,
            mob_keys,
        )

    report = RenumberReport(
        n_reference_residues=len(ref_keys),
        n_mobile_residues=len(mob_keys),
        n_mapped=len(mapping),
        n_mismatch_at_aligned_positions=n_mismatch,
        n_mobile_only_gaps_in_reference=n_mob_only,
        n_filled_unmapped_mobile=n_filled,
        aligned_reference_sequence=aln_ref,
        aligned_mobile_sequence=aln_mob,
    )
    return mapping, report


def _fill_mobile_gaps_in_map(
    mapping: Dict[ResidueKey, int],
    aln_ref: str,
    aln_mob: str,
    ref_keys: Sequence[ResidueKey],
    mob_keys: Sequence[ResidueKey],
) -> Tuple[Dict[ResidueKey, int], int]:
    """Assign resseq to mobile residues that sit in ref-gap columns; return updated map and count."""
    i_ref = i_mob = 0
    left_rs: Optional[int] = None
    cur_mob_idx: List[int] = []
    filled = 0

    def flush(right_rs: Optional[int]) -> None:
        nonlocal filled, mapping, cur_mob_idx
        if not cur_mob_idx:
            return
        nums = _assign_filled_resseqs(len(cur_mob_idx), left_rs, right_rs)
        for idx, rseq in zip(cur_mob_idx, nums):
            mapping[mob_keys[idx]] = rseq
            filled += 1
        cur_mob_idx = []

    for a, b in zip(aln_ref, aln_mob):
        if a != "-" and b != "-":
            flush(ref_keys[i_ref].resseq)
            left_rs = ref_keys[i_ref].resseq
            i_ref += 1
            i_mob += 1
        elif a != "-" and b == "-":
            i_ref += 1
        else:
            cur_mob_idx.append(i_mob)
            i_mob += 1
    flush(None)

    return mapping, filled


def renumber_pdb_atom_lines(
    lines: Sequence[str],
    mob_to_ref_resseq: Dict[ResidueKey, int],
) -> List[str]:
    """Return new lines with residue sequence numbers rewritten where mapping applies."""
    out: List[str] = []
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")) or len(line) < 27:
            out.append(line)
            continue
        resname = line[17:20].strip().upper()
        if resname not in THREE_TO_ONE:
            out.append(line)
            continue
        chain = line[21] if len(line) > 21 else " "
        icode = line[26] if len(line) > 26 else " "
        try:
            resseq = int(line[22:26])
        except ValueError:
            out.append(line)
            continue
        key = ResidueKey(chain=chain, resseq=resseq, icode=icode if icode != "" else " ")
        new_num = mob_to_ref_resseq.get(key)
        if new_num is None:
            out.append(line)
            continue
        new_field = _format_resseq_field(new_num)
        # Columns 23–26 (0-based 22:26)
        out.append(line[:22] + new_field + line[26:])
    return out


def iter_pdb_models(lines: Sequence[str]) -> Iterator[Tuple[int, int]]:
    """
    Yield (start_idx, end_idx_inclusive) for each MODEL ... ENDMDL block.
    If no MODEL records, yield (0, len-1) for whole file.
    """
    if not any(line.startswith("MODEL") for line in lines):
        yield 0, len(lines) - 1
        return
    start = None
    for i, line in enumerate(lines):
        if line.startswith("MODEL"):
            start = i
        elif line.startswith("ENDMDL") and start is not None:
            yield start, i
            start = None


def renumber_pdb_file(
    mobile_pdb: PathLike,
    output_pdb: PathLike,
    mob_to_ref_resseq: Dict[ResidueKey, int],
) -> None:
    """Rewrite ``mobile_pdb`` to ``output_pdb`` applying residue number mapping (all MODELs)."""
    in_path = Path(mobile_pdb)
    out_path = Path(output_pdb)
    with in_path.open("r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    out_lines: List[str] = []
    last_end = -1
    for start, end in iter_pdb_models(lines):
        out_lines.extend(lines[last_end + 1 : start])
        block = lines[start : end + 1]
        out_lines.extend(renumber_pdb_atom_lines(block, mob_to_ref_resseq))
        last_end = end
    out_lines.extend(lines[last_end + 1 :])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.writelines(out_lines)


def align_and_renumber_pdb(
    reference_pdb: PathLike,
    mobile_pdb: PathLike,
    output_pdb: PathLike,
    *,
    chain_id: Optional[str] = None,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
    fill_unmapped_mobile: bool = False,
) -> RenumberReport:
    """
    Align mobile chain sequence to reference chain sequence and write a new PDB where
    mobile residue numbers match the reference numbering for aligned positions.

    Parameters
    ----------
    reference_pdb
        PDB whose residue *numbers* define the target scheme (e.g. mut, or WT).
    mobile_pdb
        Trajectory or structure to renumber (e.g. WT when reference is mut).
    output_pdb
        Path for the renumbered copy (atoms unchanged; only resSeq fields updated where mapped).
    chain_id
        PDB chain to use on both files. If ``None`` and the file has exactly one non-empty
        chain, that chain is used; otherwise this must be set (e.g. ``'A'``).
    fill_unmapped_mobile
        If True, assign consecutive ``resSeq`` to mobile residues in alignment columns where
        the reference has a gap (WT-only segment vs mut deletion, or N/C overhang on mobile).

    Returns
    -------
    RenumberReport
        Alignment statistics (check ``n_mismatch_at_aligned_positions`` and unmapped counts).

    Notes
    -----
    - If ``fill_unmapped_mobile`` is False, mobile residues opposite a reference gap keep
      their original ``resSeq``.
    - Non-standard residues (not in ``THREE_TO_ONE``) are not mapped and are copied unchanged.
    """
    ref_c = _default_chain_if_needed(reference_pdb, chain_id)
    mob_c = _default_chain_if_needed(mobile_pdb, chain_id)
    ref_keys, ref_seq = ordered_residue_run(reference_pdb, chain_id=ref_c)
    mob_keys, mob_seq = ordered_residue_run(mobile_pdb, chain_id=mob_c)
    if not ref_keys or not mob_keys:
        raise ValueError("Empty residue list for reference or mobile (check chain_id).")

    mapping, report = build_mobile_to_reference_map(
        ref_keys,
        ref_seq,
        mob_keys,
        mob_seq,
        match=match,
        mismatch=mismatch,
        gap=gap,
        fill_unmapped_mobile=fill_unmapped_mobile,
    )
    renumber_pdb_file(mobile_pdb, output_pdb, mapping)
    return report


def renumber_many_to_reference(
    reference_pdb: PathLike,
    mobile_pdbs: Iterable[PathLike],
    output_dir: PathLike,
    *,
    chain_id: Optional[str] = None,
    suffix: str = "_renumbered_to_ref.pdb",
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
    fill_unmapped_mobile: bool = False,
) -> Dict[str, RenumberReport]:
    """
    Renumber multiple mobile PDBs to the same reference; write under ``output_dir``.

    Output names: ``{stem}{suffix}`` for each input stem.
    Returns mapping stem -> RenumberReport.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    reports: Dict[str, RenumberReport] = {}
    for mob in mobile_pdbs:
        p = Path(mob)
        outp = out_dir / f"{p.stem}{suffix}"
        rep = align_and_renumber_pdb(
            reference_pdb,
            p,
            outp,
            chain_id=chain_id,
            match=match,
            mismatch=mismatch,
            gap=gap,
            fill_unmapped_mobile=fill_unmapped_mobile,
        )
        reports[p.stem] = rep
    return reports
