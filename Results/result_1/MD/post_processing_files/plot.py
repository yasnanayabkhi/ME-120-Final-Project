#!/usr/bin/env python3
"""
Turn the XVG outputs from GROMACS MD into publication-ready figures.


Reads:  ./*.xvg
Writes: ./figures/Dacomitinib_FtsZ_<analysis>.png


Run from this directory:
   conda activate md_plot
   python plot.py


Pure XVG parsing (no dependencies beyond matplotlib/numpy).
"""


from __future__ import annotations


import argparse
from pathlib import Path


import matplotlib.pyplot as plt
import numpy as np


SYSTEM = "Dacomitinib_FtsZ"
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)


# Global style — no gridlines, clean look
plt.rcParams.update({
    "axes.grid": False,
    "font.family": "sans-serif",
    "axes.spines.top": False,
    "axes.spines.right": False,
})




def load_xvg(path: Path) -> tuple[np.ndarray, list[str]]:
    """
    Read an XVG file as an Nx2 (or NxK) array of floats.
    Also extract legend/label information from @ lines.
    Returns: (data array, legend list)
    """
    rows = []
    legends = []


    for line in path.read_text().splitlines():
        line = line.strip()


        if line.startswith("@") and "legend" in line:
            parts = line.split('"')
            if len(parts) >= 2:
                legends.append(parts[1])


        if not line or line.startswith(("#", "@")):
            continue


        try:
            rows.append([float(x) for x in line.split()])
        except ValueError:
            continue


    return np.asarray(rows), legends




def save_fig(fig, name: str) -> None:
    """Save figure to figures/ directory."""
    out = FIG / f"{SYSTEM}_{name}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.relative_to(HERE)}")




def plot_rmsd() -> None:
    """Plot RMSD trajectories for FtsZ and Dacomitinib with mean lines."""
    files = {
        "FtsZ Backbone": "rmsd_FtsZ",
        "Dacomitinib": "rmsd_LIG",
    }


    fig, ax = plt.subplots(figsize=(7, 4))
    found_any = False
    colors = ["#1f77b4", "#ff7f0e"]


    for (label, fname), color in zip(files.items(), colors):
        path = HERE / f"{fname}.xvg"
        if not path.exists():
            continue


        d, _ = load_xvg(path)
        ax.plot(d[:, 0], d[:, 1], label=label, linewidth=1.5, color=color)
        mean = np.mean(d[:, 1])
        ax.axhline(mean, color=color, linestyle='--', linewidth=1.2,
                   alpha=0.7, label=f"{label} mean ({mean:.3f} nm)")
        found_any = True


    if not found_any:
        return


    ax.set_xlabel("Time (ps)", fontsize=11, fontweight="bold")
    ax.set_ylabel("RMSD (nm)", fontsize=11, fontweight="bold")
    ax.set_title("FtsZ and Dacomitinib Binding Pose Stability", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    save_fig(fig, "rmsd")




def plot_rmsf() -> None:
    """Plot per-residue flexibility of FtsZ."""
    path = HERE / "rmsf.xvg"
    if not path.exists():
        return


    d, _ = load_xvg(path)


    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(d[:, 0], d[:, 1], linewidth=1.5, color="#1f77b4")
    mean = np.mean(d[:, 1])
    ax.axhline(mean, color='red', linestyle='--', linewidth=1.2,
               alpha=0.7, label=f"Mean RMSF ({mean:.3f} nm)")
    ax.set_xlabel("FtsZ Residue Number", fontsize=11, fontweight="bold")
    ax.set_ylabel("RMSF (nm)", fontsize=11, fontweight="bold")
    ax.set_title("FtsZ Per-Residue Flexibility", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    save_fig(fig, "rmsf")




def plot_hbond() -> None:
    """Plot hydrogen bond dynamics between FtsZ and Dacomitinib with mean line."""
    path = HERE / "hbond_protein_lig.xvg"
    if not path.exists():
        return


    d, _ = load_xvg(path)


    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(d[:, 0], d[:, 1], label="H-bonds", linewidth=1.5, color="#2ca02c")
    mean = np.mean(d[:, 1])
    ax.axhline(mean, color='red', linestyle='--', linewidth=1.2,
               alpha=0.7, label=f"Mean ({mean:.2f})")
    ax.set_xlabel("Time (ps)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Number of H-Bonds", fontsize=11, fontweight="bold")
    ax.set_title("FtsZ-Dacomitinib H-Bond Formation Dynamics", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    save_fig(fig, "hbond_protein_lig")




def plot_interaction_energy() -> None:
    """Plot interaction energy between FtsZ and Dacomitinib with mean lines."""
    path = HERE / "interaction_energy.xvg"
    if not path.exists():
        return


    d, _ = load_xvg(path)


    fig, ax = plt.subplots(figsize=(7, 4))


    if d.shape[1] >= 3:
        ax.plot(d[:, 0], d[:, 1], label="Coulombic (SR)", linewidth=1.5, color="#ff7f0e")
        ax.plot(d[:, 0], d[:, 2], label="LJ (SR)", linewidth=1.5, color="#d62728")
        total = d[:, 1] + d[:, 2]
        ax.plot(d[:, 0], total, label="Total", linewidth=1.5,
                linestyle='--', color="#1f77b4", alpha=0.7)
        mean_total = np.mean(total)
        ax.axhline(mean_total, color="#1f77b4", linestyle=':', linewidth=1.2,
                   alpha=0.7, label=f"Total mean ({mean_total:.1f} kJ/mol)")


    ax.set_xlabel("Time (ps)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Energy (kJ/mol)", fontsize=11, fontweight="bold")
    ax.set_title("FtsZ-Dacomitinib Interaction Energy", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    save_fig(fig, "interaction_energy")




def plot_mindist() -> None:
    """Plot minimum contact distance between FtsZ and Dacomitinib."""
    path = HERE / "mindist.xvg"
    if not path.exists():
        return


    d, _ = load_xvg(path)


    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(d[:, 0], d[:, 1], linewidth=1.5, color="#9467bd")
    mean = np.mean(d[:, 1])
    ax.axhline(mean, color='red', linestyle='--', linewidth=1.2,
               alpha=0.7, label=f"Mean ({mean:.3f} nm)")
    ax.set_xlabel("Time (ps)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Distance (nm)", fontsize=11, fontweight="bold")
    ax.set_title("FtsZ-Dacomitinib Minimum Contact Distance", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    save_fig(fig, "mindist")




def plot_sasa() -> None:
    """Plot solvent-accessible surface area of Dacomitinib."""
    path = HERE / "sasa.xvg"
    if not path.exists():
        return


    d, _ = load_xvg(path)


    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(d[:, 0], d[:, 1], linewidth=1.5, color="#8c564b")
    mean = np.mean(d[:, 1])
    ax.axhline(mean, color='red', linestyle='--', linewidth=1.2,
               alpha=0.7, label=f"Mean ({mean:.3f} nm²)")
    ax.set_xlabel("Time (ps)", fontsize=11, fontweight="bold")
    ax.set_ylabel("SASA (nm²)", fontsize=11, fontweight="bold")
    ax.set_title("Dacomitinib Burial (Solvent-Accessible Surface Area)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    save_fig(fig, "sasa")




def plot_distance() -> None:
    """Plot center of mass distance between FtsZ and Dacomitinib."""
    path = HERE / "distance.xvg"
    if not path.exists():
        return


    d, _ = load_xvg(path)


    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(d[:, 0], d[:, 1], linewidth=1.5, color="#e377c2")
    mean = np.mean(d[:, 1])
    ax.axhline(mean, color='red', linestyle='--', linewidth=1.2,
               alpha=0.7, label=f"Mean ({mean:.3f} nm)")
    ax.set_xlabel("Time (ps)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Distance (nm)", fontsize=11, fontweight="bold")
    ax.set_title("FtsZ-Dacomitinib Center of Mass Distance", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    save_fig(fig, "distance")




def plot_gyrate() -> None:
    """Plot Radius of Gyration of FtsZ."""
    path = HERE / "gyrate.xvg"
    if not path.exists():
        return


    d, _ = load_xvg(path)


    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(d[:, 0], d[:, 1], linewidth=1.5, color="#1f77b4")
    mean = np.mean(d[:, 1])
    ax.axhline(mean, color='red', linestyle='--', linewidth=1.2,
               alpha=0.7, label=f"Mean ({mean:.3f} nm)")
    ax.set_xlabel("Time (ps)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Radius of Gyration (nm)", fontsize=11, fontweight="bold")
    ax.set_title("FtsZ Radius of Gyration", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    save_fig(fig, "gyrate")




def plot_rmsd_dist() -> None:
    """Plot RMSD distribution from clustering analysis of Dacomitinib poses."""
    path = HERE / "rmsd-dist.xvg"
    if not path.exists():
        return


    d, _ = load_xvg(path)


    fig, ax = plt.subplots(figsize=(7, 4))
    if d.shape[1] >= 2:
        ax.plot(d[:, 0], d[:, 1], linewidth=1.5, color="#bcbd22",
                marker='o', markersize=3)
    ax.set_xlabel("RMSD (nm)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Count", fontsize=11, fontweight="bold")
    ax.set_title("Dacomitinib Binding Pose Clustering: RMSD Distribution",
                 fontsize=12, fontweight="bold")
    save_fig(fig, "rmsd_dist")




def main() -> None:
    """Generate all available plots."""
    argparse.ArgumentParser(
        description="Plot GROMACS XVG files for Dacomitinib-FtsZ MD simulation"
    ).parse_args()


    print(f"\nGenerating plots from {HERE}/\n")


    plot_rmsd()
    plot_rmsf()
    plot_hbond()
    plot_interaction_energy()
    plot_mindist()
    plot_sasa()
    plot_distance()
    plot_gyrate()
    plot_rmsd_dist()


    print(f"\nDone! Figures saved to {FIG.relative_to(HERE)}/\n")




if __name__ == "__main__":
    main()
