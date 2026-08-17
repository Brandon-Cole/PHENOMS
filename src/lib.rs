//! Baker–Hubbard: r(H···A) < 0.25 nm, ∠D–H–A > 120° (at H: vectors d−h and a−h).
//! Batch over all frames; returns (frame, donor, hydrogen, acceptor) for each bond.

use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray1};
use pyo3::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

/// Baker-Hubbard criteria: r(H..A) < 0.25 nm, angle D-H-A > 120°.
const R_HA_MAX_NM: f64 = 0.25;
const R_HA_MAX_NM_SQ: f64 = R_HA_MAX_NM * R_HA_MAX_NM;
/// Use spatial pruning only when the brute-force search space is large enough to amortize grid setup.
const SPATIAL_PAIR_PRODUCT_THRESHOLD: usize = 100_000;

#[inline(always)]
fn cell_index(coord: f64, cell_size: f64) -> i32 {
    (coord / cell_size).floor() as i32
}

type CellKey = (i32, i32, i32);

/// Baker–Hubbard: r(H···A) < 0.25 nm, angle D–H–A > 120°.
/// Angle at H is between vectors H→D (d − h) and H→A (a − h), matching MDTraj.
#[inline(always)]
fn is_hbond_xyz(xyz: &[f64], d_ptr: usize, h_ptr: usize, a_ptr: usize) -> bool {
    let hax = xyz[a_ptr] - xyz[h_ptr];
    let hay = xyz[a_ptr + 1] - xyz[h_ptr + 1];
    let haz = xyz[a_ptr + 2] - xyz[h_ptr + 2];
    let r_ha_sq = hax * hax + hay * hay + haz * haz;
    if r_ha_sq >= R_HA_MAX_NM_SQ || r_ha_sq <= 0.0 {
        return false;
    }

    let hdx = xyz[d_ptr] - xyz[h_ptr];
    let hdy = xyz[d_ptr + 1] - xyz[h_ptr + 1];
    let hdz = xyz[d_ptr + 2] - xyz[h_ptr + 2];
    let r_hd_sq = hdx * hdx + hdy * hdy + hdz * hdz;
    if r_hd_sq <= 0.0 {
        return false;
    }

    // ∠D–H–A > 120°  <=>  dot(hd,ha) < -0.5*|hd||ha|
    let dot = hdx * hax + hdy * hay + hdz * haz;
    dot < 0.0 && (dot * dot) > 0.25 * r_hd_sq * r_ha_sq
}

fn build_acceptor_grid(
    xyz: &[f64],
    base: usize,
    acceptors: &[usize],
    cell_size: f64,
) -> HashMap<CellKey, Vec<usize>> {
    let mut grid: HashMap<CellKey, Vec<usize>> = HashMap::new();
    for &acc in acceptors {
        let a_ptr = base + acc * 3;
        let key = (
            cell_index(xyz[a_ptr], cell_size),
            cell_index(xyz[a_ptr + 1], cell_size),
            cell_index(xyz[a_ptr + 2], cell_size),
        );
        grid.entry(key).or_default().push(acc);
    }
    grid
}

fn candidate_acceptors_near_h(
    grid: &HashMap<CellKey, Vec<usize>>,
    xyz: &[f64],
    h_ptr: usize,
    cell_size: f64,
) -> Vec<usize> {
    let cx = cell_index(xyz[h_ptr], cell_size);
    let cy = cell_index(xyz[h_ptr + 1], cell_size);
    let cz = cell_index(xyz[h_ptr + 2], cell_size);
    let mut out = Vec::new();
    for dx in -1..=1 {
        for dy in -1..=1 {
            for dz in -1..=1 {
                if let Some(accs) = grid.get(&(cx + dx, cy + dy, cz + dz)) {
                    out.extend(accs.iter().copied());
                }
            }
        }
    }
    out
}

fn collect_bonds_brute_frame(
    xyz: &[f64],
    frame: usize,
    base: usize,
    valid_pairs: &[(usize, usize)],
    valid_acceptors: &[usize],
    out: &mut Vec<u32>,
) {
    for &(d, h) in valid_pairs {
        let d_ptr = base + d * 3;
        let h_ptr = base + h * 3;
        for &acc in valid_acceptors {
            if acc == d {
                continue;
            }
            let a_ptr = base + acc * 3;
            if is_hbond_xyz(xyz, d_ptr, h_ptr, a_ptr) {
                out.push(frame as u32);
                out.push(d as u32);
                out.push(h as u32);
                out.push(acc as u32);
            }
        }
    }
}

fn collect_bonds_spatial_frame(
    xyz: &[f64],
    frame: usize,
    base: usize,
    valid_pairs: &[(usize, usize)],
    valid_acceptors: &[usize],
    out: &mut Vec<u32>,
) {
    let grid = build_acceptor_grid(xyz, base, valid_acceptors, R_HA_MAX_NM);
    for &(d, h) in valid_pairs {
        let d_ptr = base + d * 3;
        let h_ptr = base + h * 3;
        let candidates = candidate_acceptors_near_h(&grid, xyz, h_ptr, R_HA_MAX_NM);
        for acc in candidates {
            if acc == d {
                continue;
            }
            let a_ptr = base + acc * 3;
            if is_hbond_xyz(xyz, d_ptr, h_ptr, a_ptr) {
                out.push(frame as u32);
                out.push(d as u32);
                out.push(h as u32);
                out.push(acc as u32);
            }
        }
    }
}

fn validate_inputs(
    xyz: &[f64],
    n_frames: usize,
    n_atoms: usize,
    donor_h_pairs: Vec<(usize, usize)>,
    acceptor_indices: Vec<usize>,
) -> Result<(Vec<(usize, usize)>, Vec<usize>), PyErr> {
    if xyz.len() != n_frames * n_atoms * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "xyz_flat length must equal n_frames * n_atoms * 3",
        ));
    }
    let valid_pairs: Vec<(usize, usize)> = donor_h_pairs
        .into_iter()
        .filter(|&(d, h)| d < n_atoms && h < n_atoms)
        .collect();
    let valid_acceptors: Vec<usize> = acceptor_indices
        .into_iter()
        .filter(|&a| a < n_atoms)
        .collect();
    Ok((valid_pairs, valid_acceptors))
}

fn bonds_to_pyarray(py: Python<'_>, out: Vec<u32>) -> PyResult<Py<PyArray2<u32>>> {
    let n_rows = out.len() / 4;
    let arr = if n_rows == 0 {
        PyArray2::zeros_bound(py, [0, 4], false)
    } else {
        let a2 = Array2::from_shape_vec((n_rows, 4), out)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("bond buffer shape: {e}")))?;
        a2.into_pyarray_bound(py)
    };
    Ok(arr.unbind())
}

fn run_kernel<F>(
    py: Python<'_>,
    xyz: &[f64],
    n_frames: usize,
    n_atoms: usize,
    valid_pairs: &[(usize, usize)],
    valid_acceptors: &[usize],
    mut collect_frame: F,
) -> PyResult<Py<PyArray2<u32>>>
where
    F: FnMut(&[f64], usize, usize, &[(usize, usize)], &[usize], &mut Vec<u32>),
{
    let stride = n_atoms * 3;
    let mut out: Vec<u32> = Vec::new();
    out.reserve(n_frames.saturating_mul(valid_pairs.len() / 4).saturating_mul(4));

    for frame in 0..n_frames {
        let base = frame * stride;
        collect_frame(xyz, frame, base, valid_pairs, valid_acceptors, &mut out);
    }
    bonds_to_pyarray(py, out)
}

fn run_kernel_threaded<F>(
    py: Python<'_>,
    xyz: &[f64],
    n_frames: usize,
    n_atoms: usize,
    valid_pairs: &[(usize, usize)],
    valid_acceptors: &[usize],
    n_threads: usize,
    collect_frame: F,
) -> PyResult<Py<PyArray2<u32>>>
where
    F: Fn(&[f64], usize, usize, &[(usize, usize)], &[usize], &mut Vec<u32>) + Sync,
{
    let stride = n_atoms * 3;
    let threads = n_threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("rayon pool build failed: {e}")))?;

    let out: Vec<u32> = pool.install(|| {
        (0..n_frames)
            .into_par_iter()
            .flat_map_iter(|frame| {
                let base = frame * stride;
                let mut frame_out: Vec<u32> = Vec::new();
                collect_frame(xyz, frame, base, valid_pairs, valid_acceptors, &mut frame_out);
                frame_out
            })
            .collect()
    });
    bonds_to_pyarray(py, out)
}

/// Brute-force Baker–Hubbard (reference implementation for regression tests).
#[pyfunction]
fn run_baker_hubbard_brute(
    py: Python<'_>,
    xyz_flat: PyReadonlyArray1<f64>,
    n_frames: usize,
    n_atoms: usize,
    donor_h_pairs: Vec<(usize, usize)>,
    acceptor_indices: Vec<usize>,
) -> PyResult<Py<PyArray2<u32>>> {
    let xyz = xyz_flat.as_slice()?;
    let (valid_pairs, valid_acceptors) =
        validate_inputs(xyz, n_frames, n_atoms, donor_h_pairs, acceptor_indices)?;
    run_kernel(py, xyz, n_frames, n_atoms, &valid_pairs, &valid_acceptors, collect_bonds_brute_frame)
}

/// Run Baker-Hubbard over all frames. Uses spatial pruning for large search spaces
/// and brute force for smaller ones (same geometry either way).
#[pyfunction]
fn run_baker_hubbard(
    py: Python<'_>,
    xyz_flat: PyReadonlyArray1<f64>,
    n_frames: usize,
    n_atoms: usize,
    donor_h_pairs: Vec<(usize, usize)>,
    acceptor_indices: Vec<usize>,
) -> PyResult<Py<PyArray2<u32>>> {
    let xyz = xyz_flat.as_slice()?;
    let (valid_pairs, valid_acceptors) =
        validate_inputs(xyz, n_frames, n_atoms, donor_h_pairs, acceptor_indices)?;
    let search_space = valid_pairs.len().saturating_mul(valid_acceptors.len());
    if search_space >= SPATIAL_PAIR_PRODUCT_THRESHOLD {
        run_kernel(
            py,
            xyz,
            n_frames,
            n_atoms,
            &valid_pairs,
            &valid_acceptors,
            collect_bonds_spatial_frame,
        )
    } else {
        run_kernel(
            py,
            xyz,
            n_frames,
            n_atoms,
            &valid_pairs,
            &valid_acceptors,
            collect_bonds_brute_frame,
        )
    }
}

/// Threaded Baker–Hubbard (adaptive spatial vs brute per frame batch).
#[pyfunction]
fn run_baker_hubbard_with_threads(
    py: Python<'_>,
    xyz_flat: PyReadonlyArray1<f64>,
    n_frames: usize,
    n_atoms: usize,
    donor_h_pairs: Vec<(usize, usize)>,
    acceptor_indices: Vec<usize>,
    n_threads: usize,
) -> PyResult<Py<PyArray2<u32>>> {
    let xyz = xyz_flat.as_slice()?;
    let (valid_pairs, valid_acceptors) =
        validate_inputs(xyz, n_frames, n_atoms, donor_h_pairs, acceptor_indices)?;
    let search_space = valid_pairs.len().saturating_mul(valid_acceptors.len());
    let collect = if search_space >= SPATIAL_PAIR_PRODUCT_THRESHOLD {
        collect_bonds_spatial_frame as fn(&[f64], usize, usize, &[(usize, usize)], &[usize], &mut Vec<u32>)
    } else {
        collect_bonds_brute_frame
    };
    run_kernel_threaded(
        py,
        xyz,
        n_frames,
        n_atoms,
        &valid_pairs,
        &valid_acceptors,
        n_threads,
        collect,
    )
}

#[pymodule]
fn phenoms_hbond_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_baker_hubbard, m)?)?;
    m.add_function(wrap_pyfunction!(run_baker_hubbard_brute, m)?)?;
    m.add_function(wrap_pyfunction!(run_baker_hubbard_with_threads, m)?)?;
    Ok(())
}
