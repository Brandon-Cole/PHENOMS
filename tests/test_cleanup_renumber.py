"""Tests for phenoms.cleanup sequence alignment and PDB renumbering."""

from pathlib import Path

from phenoms.cleanup import (
    ResidueKey,
    align_and_renumber_pdb,
    build_mobile_to_reference_map,
    needleman_wunsch,
    renumber_pdb_atom_lines,
)


def test_needleman_wunsch_deletion():
    a, b = needleman_wunsch("ACDEFG", "ADEFG", match=2, mismatch=-1, gap=-2)
    assert len(a) == len(b)
    assert "-" in a or "-" in b


def test_build_map_simple():
    ref_keys = [
        ResidueKey(" ", 10),
        ResidueKey(" ", 11),
        ResidueKey(" ", 12),
        ResidueKey(" ", 13),
    ]
    mob_keys = [
        ResidueKey(" ", 1),
        ResidueKey(" ", 2),
        ResidueKey(" ", 3),
    ]
    mapping, rep = build_mobile_to_reference_map(ref_keys, "ABCD", mob_keys, "ACD")
    assert mapping[ResidueKey(" ", 1)] == 10
    assert mapping[ResidueKey(" ", 2)] == 12
    assert mapping[ResidueKey(" ", 3)] == 13
    assert rep.n_mapped == 3
    assert rep.n_filled_unmapped_mobile == 0


def test_build_map_fill_deletion_in_reference():
    # Reference = mut (deletion): A C D — B missing vs mobile WT A B C D
    ref_keys = [ResidueKey(" ", 100), ResidueKey(" ", 102), ResidueKey(" ", 103)]
    mob_keys = [ResidueKey(" ", 1), ResidueKey(" ", 2), ResidueKey(" ", 3), ResidueKey(" ", 4)]
    mapping, rep = build_mobile_to_reference_map(
        ref_keys,
        "ACD",
        mob_keys,
        "ABCD",
        fill_unmapped_mobile=True,
    )
    assert mapping[ResidueKey(" ", 1)] == 100
    assert mapping[ResidueKey(" ", 2)] == 101  # WT-only vs mut gap
    assert mapping[ResidueKey(" ", 3)] == 102
    assert mapping[ResidueKey(" ", 4)] == 103
    assert rep.n_filled_unmapped_mobile == 1


def test_renumber_pdb_atom_lines():
    mob_to_ref = {ResidueKey("A", 5): 100}
    line = "ATOM      1  N   ALA A   5      10.0   20.0   30.0  1.00  0.00           N\n"
    # Note: chain at col 22 = 'A', resseq 22:26
    out = renumber_pdb_atom_lines([line], mob_to_ref)
    assert " 100" in out[0][22:26] or out[0][22:26].strip() == "100"


def test_align_and_renumber_minimal_pdb(tmp_path: Path):
    ref = tmp_path / "ref.pdb"
    mob = tmp_path / "mob.pdb"
    out = tmp_path / "out.pdb"
    ref.write_text(
        "MODEL        1\n"
        "ATOM      1  N   ALA A  10       0.0    0.0    0.0  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A  10       1.0    0.0    0.0  1.00  0.00           C\n"
        "ATOM      3  N   GLY A  11       2.0    0.0    0.0  1.00  0.00           N\n"
        "ENDMDL\n"
        "END\n",
        encoding="utf-8",
    )
    mob.write_text(
        "MODEL        1\n"
        "ATOM      1  N   ALA A  99       0.0    0.0    0.0  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A  99       1.0    0.0    0.0  1.00  0.00           C\n"
        "ATOM      3  N   GLY A 100       2.0    0.0    0.0  1.00  0.00           N\n"
        "ENDMDL\n"
        "END\n",
        encoding="utf-8",
    )
    rep = align_and_renumber_pdb(ref, mob, out, chain_id="A")
    assert rep.n_mapped == 2
    text = out.read_text(encoding="utf-8")
    assert "  10" in text  # ALA renumbered to 10
    assert "  11" in text  # GLY renumbered to 11
