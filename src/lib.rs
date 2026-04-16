//! Baker–Hubbard: r(H···A) < 0.25 nm, ∠D–H–A > 120° (at H: vectors d−h and a−h).
//! Batch over all frames; returns (frame, donor, hydrogen, acceptor) for each bond.

use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray1};
use pyo3::prelude::*;
use rayon::prelude::*;

/// Baker-Hubbard criteria: r(H..A) < 0.25 nm, angle D-H-A > 120°.
const R_HA_MAX_NM: f64 = 0.25;
const R_HA_MAX_NM_SQ: f64 = R_HA_MAX_NM * R_HA_MAX_NM;

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
    // Compare squared magnitudes to avoid sqrt:
    // dot^2 > 0.25 * r_hd_sq * r_ha_sq and dot < 0.
    let dot = hdx * hax + hdy * hay + hdz * haz;
    dot < 0.0 && (dot * dot) > 0.25 * r_hd_sq * r_ha_sq
}

/// Run Baker-Hubbard over all frames. Coordinates in nm, flat row-major (n_frames * n_atoms * 3).
///
/// Returns (n_bonds, 4) array: each row is (frame, donor, hydrogen, acceptor).
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
    if xyz.len() != n_frames * n_atoms * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "xyz_flat length must equal n_frames * n_atoms * 3",
        ));
    }

    let mut out: Vec<u32> = Vec::new();
    // Heuristic reserve to reduce reallocations.
    out.reserve(n_frames.saturating_mul(donor_h_pairs.len() / 4).saturating_mul(4));
    let stride = n_atoms * 3;

    let valid_pairs: Vec<(usize, usize)> = donor_h_pairs
        .into_iter()
        .filter(|&(d, h)| d < n_atoms && h < n_atoms)
        .collect();
    let valid_acceptors: Vec<usize> = acceptor_indices
        .into_iter()
        .filter(|&a| a < n_atoms)
        .collect();

    for frame in 0..n_frames {
        let base = frame * stride;
        for &(d, h) in &valid_pairs {
            let d_ptr = base + d * 3;
            let h_ptr = base + h * 3;
            for &acc in &valid_acceptors {
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

/// Threaded variant of Baker–Hubbard over all frames.
/// Uses rayon and honors `n_threads` (>=1).
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
    if xyz.len() != n_frames * n_atoms * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "xyz_flat length must equal n_frames * n_atoms * 3",
        ));
    }

    let stride = n_atoms * 3;
    let valid_pairs: Vec<(usize, usize)> = donor_h_pairs
        .into_iter()
        .filter(|&(d, h)| d < n_atoms && h < n_atoms)
        .collect();
    let valid_acceptors: Vec<usize> = acceptor_indices
        .into_iter()
        .filter(|&a| a < n_atoms)
        .collect();

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
                for &(d, h) in &valid_pairs {
                    let d_ptr = base + d * 3;
                    let h_ptr = base + h * 3;
                    for &acc in &valid_acceptors {
                        if acc == d {
                            continue;
                        }
                        let a_ptr = base + acc * 3;
                        if is_hbond_xyz(xyz, d_ptr, h_ptr, a_ptr) {
                            frame_out.push(frame as u32);
                            frame_out.push(d as u32);
                            frame_out.push(h as u32);
                            frame_out.push(acc as u32);
                        }
                    }
                }
                frame_out
            })
            .collect()
    });

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

#[pymodule]
fn phenoms_hbond_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_baker_hubbard, m)?)?;
    m.add_function(wrap_pyfunction!(run_baker_hubbard_with_threads, m)?)?;
    Ok(())
}
