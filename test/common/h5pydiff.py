#!/usr/bin/env python3
"""Compare two HDF5 files for structural and value equivalence."""

import argparse
import sys

import h5py
import numpy as np


def _attr_dtype(obj, name):
    return h5py.h5a.open(obj.id, name.encode()).dtype


def _is_complex_compound(dtype):
    """Return True if dtype is an H5T_COMPOUND with 're' and 'im' float fields."""
    return (
        dtype.names is not None
        and set(dtype.names) == {"re", "im"}
        and np.issubdtype(dtype["re"], np.floating)
        and np.issubdtype(dtype["im"], np.floating)
    )


def _threshold_msg(rel_ok, abs_ok, rel_tol, abs_tol):
    """Build a failure reason string describing which thresholds were breached."""
    if rel_tol is not None and abs_tol is not None:
        if not rel_ok and not abs_ok:
            return "relative and absolute thresholds breached"
        if not rel_ok:
            return "relative threshold breached"
        return "absolute threshold breached"
    if rel_tol is not None:
        return "relative threshold breached"
    return "absolute threshold breached"


def _error_details(n_fail, n_total, max_abs, max_rel):
    """Build the count + error-value suffix appended to a failure reason."""
    parts = []
    if n_total > 1:
        parts.append(f"{n_fail}/{n_total} elements")
    if max_abs is not None:
        parts.append(f"max abs err = {max_abs:.3e}")
    if max_rel is not None and not np.isinf(max_rel):
        parts.append(f"max rel err = {max_rel:.3e}")
    return "; ".join(parts)


def _check_real(a1, a2, rel_tol, abs_tol):
    """Compare real-valued arrays. Returns (ok, reason_or_None)."""
    if rel_tol is None and abs_tol is None:
        return bool(np.array_equal(a1, a2, equal_nan=True)), None

    a1f = np.asarray(a1, dtype=float)
    a2f = np.asarray(a2, dtype=float)
    diff = np.abs(a1f - a2f)
    ref = np.maximum(np.abs(a1f), np.abs(a2f))

    with np.errstate(invalid="ignore", divide="ignore"):
        rel_diff = np.where(ref > 0, diff / ref, np.where(diff == 0, 0.0, np.inf))

    rel_ok = True
    abs_ok = True
    if rel_tol is not None:
        rel_ok = bool(np.all(rel_diff <= rel_tol))
    if abs_tol is not None:
        abs_ok = bool(np.all(diff <= abs_tol))

    if rel_ok and abs_ok:
        return True, None

    n_total = diff.size
    if rel_tol is not None and abs_tol is not None:
        if not rel_ok and not abs_ok:
            n_fail = int(np.sum((rel_diff > rel_tol) | (diff > abs_tol)))
        elif not rel_ok:
            n_fail = int(np.sum(rel_diff > rel_tol))
        else:
            n_fail = int(np.sum(diff > abs_tol))
    elif rel_tol is not None:
        n_fail = int(np.sum(rel_diff > rel_tol))
    else:
        n_fail = int(np.sum(diff > abs_tol))

    details = _error_details(
        n_fail, n_total,
        max_abs=float(np.max(diff)),
        max_rel=float(np.max(rel_diff)),
    )
    reason = _threshold_msg(rel_ok, abs_ok, rel_tol, abs_tol)
    if details:
        reason += "; " + details
    return False, reason


def _check_complex(a1, a2, rel_tol, abs_tol):
    """Compare complex compound (re/im) arrays. Returns (ok, reason_or_None).

    Relative tolerance is scaled per-element by the mean magnitude of the two
    complex values, so each component difference must satisfy:
        |re1 - re2| <= rel_tol * norm   and   |im1 - im2| <= rel_tol * norm
    where norm = 0.5 * (|z1| + |z2|).

    For error reporting, per-element error is max(re_diff, im_diff); relative
    error is that quantity divided by norm.
    """
    if rel_tol is None and abs_tol is None:
        return bool(np.array_equal(a1, a2)), None

    re_diff = np.abs(a1["re"] - a2["re"])
    im_diff = np.abs(a1["im"] - a2["im"])
    norm = 0.5 * (
        np.sqrt(a1["re"] ** 2 + a1["im"] ** 2)
        + np.sqrt(a2["re"] ** 2 + a2["im"] ** 2)
    )
    comp_diff = np.maximum(re_diff, im_diff)  # per-element absolute error

    with np.errstate(invalid="ignore", divide="ignore"):
        comp_rel = np.where(norm > 0, comp_diff / norm,
                            np.where(comp_diff == 0, 0.0, np.inf))

    rel_ok = True
    abs_ok = True
    if rel_tol is not None:
        rel_ok = bool(np.all(re_diff <= rel_tol * norm) and np.all(im_diff <= rel_tol * norm))
    if abs_tol is not None:
        abs_ok = bool(np.all(re_diff <= abs_tol) and np.all(im_diff <= abs_tol))

    if rel_ok and abs_ok:
        return True, None

    n_total = comp_diff.size
    if rel_tol is not None and abs_tol is not None:
        if not rel_ok and not abs_ok:
            n_fail = int(np.sum((comp_rel > rel_tol) | (comp_diff > abs_tol)))
        elif not rel_ok:
            n_fail = int(np.sum(comp_rel > rel_tol))
        else:
            n_fail = int(np.sum(comp_diff > abs_tol))
    elif rel_tol is not None:
        n_fail = int(np.sum(comp_rel > rel_tol))
    else:
        n_fail = int(np.sum(comp_diff > abs_tol))

    details = _error_details(
        n_fail, n_total,
        max_abs=float(np.max(comp_diff)),
        max_rel=float(np.max(comp_rel)),
    )
    reason = _threshold_msg(rel_ok, abs_ok, rel_tol, abs_tol)
    if details:
        reason += "; " + details
    return False, reason


def _values_match(v1, v2, dtype, rel_tol, abs_tol):
    """Compare two values. Returns (ok, reason_or_None)."""
    a1 = np.asarray(v1)
    a2 = np.asarray(v2)
    if a1.shape != a2.shape:
        return False, f"shape {a1.shape} vs {a2.shape}"

    if _is_complex_compound(dtype):
        return _check_complex(a1, a2, rel_tol, abs_tol)
    if np.issubdtype(dtype, np.floating):
        return _check_real(a1, a2, rel_tol, abs_tol)
    # integers, strings, and other types: always exact
    ok = bool(np.array_equal(a1, a2))
    return ok, None


def _child_path(parent, name):
    return f"{parent.rstrip('/')}/{name}"


def compare_attrs(path, obj1, obj2, rel_tol, abs_tol, diffs):
    keys1 = set(obj1.attrs.keys())
    keys2 = set(obj2.attrs.keys())

    for k in sorted(keys1 - keys2):
        diffs.append(f"ATTR  {path}@{k}: only in file 1")
    for k in sorted(keys2 - keys1):
        diffs.append(f"ATTR  {path}@{k}: only in file 2")

    for k in sorted(keys1 & keys2):
        desc = f"{path}@{k}"
        dtype1 = _attr_dtype(obj1, k)
        dtype2 = _attr_dtype(obj2, k)
        if dtype1 != dtype2:
            diffs.append(f"ATTR  {desc}: type mismatch ({dtype1} vs {dtype2})")
        v1, v2 = obj1.attrs[k], obj2.attrs[k]
        ok, reason = _values_match(v1, v2, dtype1, rel_tol, abs_tol)
        if not ok:
            msg = f"ATTR  {desc}: value mismatch"
            if reason:
                msg += f" ({reason})"
            diffs.append(msg)


def compare_datasets(path, ds1, ds2, rel_tol, abs_tol, diffs):
    dtype1, dtype2 = ds1.dtype, ds2.dtype
    if dtype1 != dtype2:
        diffs.append(f"DSET  {path}: type mismatch ({dtype1} vs {dtype2})")

    if ds1.shape != ds2.shape:
        diffs.append(f"DSET  {path}: shape mismatch ({ds1.shape} vs {ds2.shape})")
    else:
        v1, v2 = ds1[()], ds2[()]
        ok, reason = _values_match(v1, v2, dtype1, rel_tol, abs_tol)
        if not ok:
            msg = f"DSET  {path}: value mismatch"
            if reason:
                msg += f" ({reason})"
            diffs.append(msg)

    compare_attrs(path, ds1, ds2, rel_tol, abs_tol, diffs)


def compare_groups(path, g1, g2, rel_tol, abs_tol, diffs):
    compare_attrs(path, g1, g2, rel_tol, abs_tol, diffs)

    keys1 = set(g1.keys())
    keys2 = set(g2.keys())

    for k in sorted(keys1 - keys2):
        kind = "GROUP" if isinstance(g1[k], h5py.Group) else "DSET "
        diffs.append(f"{kind} {_child_path(path, k)}: only in file 1")
    for k in sorted(keys2 - keys1):
        kind = "GROUP" if isinstance(g2[k], h5py.Group) else "DSET "
        diffs.append(f"{kind} {_child_path(path, k)}: only in file 2")

    for k in sorted(keys1 & keys2):
        child = _child_path(path, k)
        obj1, obj2 = g1[k], g2[k]

        if type(obj1) != type(obj2):
            diffs.append(f"      {child}: type mismatch (group vs dataset)")
            continue

        if isinstance(obj1, h5py.Group):
            compare_groups(child, obj1, obj2, rel_tol, abs_tol, diffs)
        elif isinstance(obj1, h5py.Dataset):
            compare_datasets(child, obj1, obj2, rel_tol, abs_tol, diffs)


def main():
    parser = argparse.ArgumentParser(
        description="Compare two HDF5 files for structural and value equivalence."
    )
    parser.add_argument("file1", help="First HDF5 file")
    parser.add_argument("file2", help="Second HDF5 file")
    parser.add_argument(
        "--relative", "-r",
        type=float,
        default=None,
        metavar="REL",
        help="Relative tolerance for float/complex comparisons. "
             "For complex, each component is compared against REL * element_norm.",
    )
    parser.add_argument(
        "--absolute", "-a",
        type=float,
        default=None,
        metavar="ABS",
        help="Absolute tolerance for float/complex comparisons.",
    )
    args = parser.parse_args()

    diffs = []
    with h5py.File(args.file1, "r") as f1, h5py.File(args.file2, "r") as f2:
        compare_groups("/", f1, f2, args.relative, args.absolute, diffs)

    for d in diffs:
        print(d)

    return 1 if diffs else 0


if __name__ == "__main__":
    sys.exit(main())
