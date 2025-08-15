import h5py
import numpy as np
from typing import Optional
from py3dgespar import run_gespar3d_stats


def load_volume(path: str, crop: Optional[int] = None) -> np.ndarray:
    """Load a knee MRI volume and optionally center-crop to a cube.

    Parameters
    ----------
    path: str
        Path to the HDF5 file containing ``kspace`` data from fastMRI.
    crop: Optional[int]

        If provided, the volume is center-cropped to ``crop``^3 voxels.

    Returns
    -------
    np.ndarray
        Complex-valued volume in image space.
    """
    with h5py.File(path, "r") as f:
        kspace = f["kspace"][()]
    # Convert each 2D slice to image space
    volume = np.fft.ifft2(kspace, norm="ortho")
    if crop is not None:
        s, h, w = volume.shape
        dim = min(crop, s, h, w)
        s0 = (s - dim) // 2
        h0 = (h - dim) // 2
        w0 = (w - dim) // 2
        volume = volume[s0:s0+dim, h0:h0+dim, w0:w0+dim]
    return volume


def main() -> None:
    data_path = (
        "testdata/fastMRI/knee_singlecoil_test~/singlecoil_test/file1000343.h5"
    )
    volume = load_volume(data_path, crop=16)
    dimlen = volume.shape[0]
    x = volume.ravel()
    k = 20
    m = x.size
    max_t = 1000
    snr = np.inf
    nmse, supp, runtime, iterations = run_gespar3d_stats(
        x, dimlen, k, m, max_t, snr, verbose=True
    )
    print(f"NMSE: {nmse}")
    print(f"Support matches: {supp}")
    print(f"Runtime: {runtime}")
    print(f"Iterations: {iterations}")


if __name__ == "__main__":
    main()
