"""
Write PDB with B-factors set from per-residue values (e.g. difference, variance) for visualization.
"""

import re


def _residue_number_from_donor_label(donor_residue):
    """Extract integer residue number from 'ALA1' or 'GLY42'."""
    m = re.search(r"\d+", str(donor_residue))
    return int(m.group()) if m else None


def write_pdb_bfactors(
    pdb_path,
    residue_to_value,
    output_path,
    value_column=None,
    model_index=0,
):
    """
    Write a PDB file with B-factors set to the given per-residue values (e.g. difference).
    Use for PyMOL/Chimera visualization of differential protection or fluctuation.

    Parameters
    ----------
    pdb_path : str
        Reference PDB file path.
    residue_to_value : dict or pd.DataFrame
        If dict: mapping residue_number (int) -> value (float).
        If DataFrame: must have 'Residue Number' and a value column (or 'Difference').
    output_path : str
        Output PDB path.
    value_column : str or None
        If residue_to_value is a DataFrame, column to use for B-factor.
        If None, uses ``Difference_clipped`` when present, else ``Difference``.
    model_index : int
        Only modify a single static model/frame (0-based). Useful for multi-model PDBs.
    """
    try:
        from Bio.PDB import PDBParser, PDBIO
    except ImportError:
        # Fallback: update B-factor column directly in fixed-width PDB text.
        # B-factor is columns 61-66 (1-indexed).
        if hasattr(residue_to_value, "columns"):
            col = value_column
            if col is None:
                col = (
                    "Difference_clipped"
                    if "Difference_clipped" in residue_to_value.columns
                    else "Difference"
                )
            if "Residue Number" not in residue_to_value.columns or col not in residue_to_value.columns:
                raise ValueError(f"DataFrame must have 'Residue Number' and '{col}'")
            mapping = dict(
                zip(
                    residue_to_value["Residue Number"].astype(int),
                    residue_to_value[col].astype(float),
                )
            )
        else:
            mapping = dict(residue_to_value)

        with open(pdb_path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        has_model_records = any(line.startswith("MODEL") for line in lines)

        # If the PDB contains many MODEL blocks (e.g., 500 frames), emit only the selected
        # static model. Otherwise the viewer will still show all frames.
        if has_model_records:
            model_blocks: list[tuple[int, int]] = []
            current_start = None
            for i, line in enumerate(lines):
                if line.startswith("MODEL"):
                    current_start = i
                elif line.startswith("ENDMDL") and current_start is not None:
                    model_blocks.append((current_start, i))
                    current_start = None

            if not model_blocks:
                raise ValueError(f"MODEL records found but no ENDMDL blocks parsed: {pdb_path}")

            model_index_int = int(model_index)
            if model_index_int < 0 or model_index_int >= len(model_blocks):
                raise ValueError(
                    f"model_index={model_index_int} out of range (PDB has {len(model_blocks)} models): {pdb_path}"
                )

            start, end = model_blocks[model_index_int]
            prefix_lines = lines[:start]
            selected_lines = lines[start : end + 1]

            # Keep only footer lines after the selected model, up to the next MODEL.
            suffix_lines: list[str] = []
            for line in lines[end + 1 :]:
                if line.startswith("MODEL"):
                    break
                suffix_lines.append(line)

            out_lines: list[str] = []
            out_lines.extend(prefix_lines)

            for line in selected_lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Residue sequence number is columns 23-26 (1-indexed).
                    res_str = line[22:26].strip()
                    try:
                        res_num = int(res_str)
                    except ValueError:
                        out_lines.append(line)
                        continue

                    value = mapping.get(res_num, 0.0)
                    try:
                        value_f = float(value)
                    except (TypeError, ValueError):
                        value_f = 0.0

                    bfac_str = f"{value_f:6.2f}"  # width=6 to match PDB column.
                    # Replace columns 61-66: indices 60..65 inclusive -> [60:66]
                    if len(line) >= 66:
                        out_lines.append(line[:60] + bfac_str + line[66:])
                    else:
                        out_lines.append(line.rstrip("\n") + bfac_str + "\n")
                else:
                    out_lines.append(line)

            out_lines.extend(suffix_lines)

            # Ensure PDB terminator exists.
            if not any(l.startswith("END") for l in out_lines[-10:]):
                out_lines.append("END\n")

            with open(output_path, "w", encoding="utf-8") as f:
                f.writelines(out_lines)
            return

        # Single-model PDB: update B-factors in-place across the whole file.
        out_lines: list[str] = []
        for line in lines:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                out_lines.append(line)
                continue

            # Residue sequence number is columns 23-26 (1-indexed).
            res_str = line[22:26].strip()
            try:
                res_num = int(res_str)
            except ValueError:
                out_lines.append(line)
                continue

            value = mapping.get(res_num, 0.0)
            try:
                value_f = float(value)
            except (TypeError, ValueError):
                value_f = 0.0

            bfac_str = f"{value_f:6.2f}"  # width=6 to match PDB column.
            # Replace columns 61-66: indices 60..65 inclusive -> [60:66]
            if len(line) >= 66:
                out_lines.append(line[:60] + bfac_str + line[66:])
            else:
                out_lines.append(line.rstrip("\n") + bfac_str + "\n")

        with open(output_path, "w", encoding="utf-8") as f:
            f.writelines(out_lines)
        return

    if hasattr(residue_to_value, "columns"):
        # DataFrame
        col = value_column
        if col is None:
            col = (
                "Difference_clipped"
                if "Difference_clipped" in residue_to_value.columns
                else "Difference"
            )
        if "Residue Number" not in residue_to_value.columns or col not in residue_to_value.columns:
            raise ValueError(f"DataFrame must have 'Residue Number' and '{col}'")
        mapping = dict(zip(residue_to_value["Residue Number"].astype(int), residue_to_value[col].astype(float)))
    else:
        mapping = dict(residue_to_value)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    models = list(structure.get_models())
    if not models:
        raise ValueError(f"No models found in PDB: {pdb_path}")
    if int(model_index) < 0 or int(model_index) >= len(models):
        raise ValueError(
            f"model_index={model_index} out of range (PDB has {len(models)} models): {pdb_path}"
        )

    selected_model = models[int(model_index)]
    for chain in selected_model:
        for residue in chain:
            res_id = residue.get_id()
            res_num = res_id[1]
            value = mapping.get(res_num, 0.0)
            for atom in residue:
                atom.set_bfactor(float(value))
    io = PDBIO()
    # Emit only the selected model to produce a static-frame PDB.
    io.set_structure(selected_model)
    io.save(output_path)


def residue_value_map_from_comparison(comparison_df, value_column="Difference_clipped"):
    """Build dict residue_number -> value from comparison DataFrame for write_pdb_bfactors."""
    if "Residue Number" not in comparison_df.columns or value_column not in comparison_df.columns:
        raise ValueError(f"comparison_df must have 'Residue Number' and '{value_column}'")
    return dict(zip(comparison_df["Residue Number"].astype(int), comparison_df[value_column].astype(float)))
