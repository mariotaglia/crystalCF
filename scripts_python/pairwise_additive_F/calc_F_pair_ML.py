# =====================================================================
# Pairwise additive free energy using Context-Aware Neural Network
# Compatible with TorchContextRegressor (2-phase + asymptotic training)
# =====================================================================

import numpy as np
import os
import joblib
import torch
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# =====================================================================
# 0. LOAD TRAINED MODEL
# =====================================================================

DEVICE = torch.device("cpu")

MODEL_DIR = os.path.expanduser(
    "~/develop/crystalCF/scripts_python/pairwise_additive_F"
)
MODEL_FILE = os.path.join(MODEL_DIR, "torch_context_model.pkl")

try:
    DEVICE = torch.device("cpu")
    checkpoint = joblib.load(MODEL_FILE)

    FEATURES = checkpoint["features"]
    F_SCALE = checkpoint["F_scale"]

    idx_d = FEATURES.index("D")
    idx_ctx = list(range(idx_d))

    hidden_dim = checkpoint.get("hidden_dim", 64)

    # --- Arquitectura EXACTA a la de entrenamiento ---
    import torch.nn as nn

    class ExponentialAsymptoticNet(nn.Module):
        def __init__(self, context_dim, hidden_dim=64):
            super().__init__()

            self.encoder = nn.Sequential(
                nn.Linear(context_dim, hidden_dim),
                nn.SiLU(),
                nn.Linear(hidden_dim, hidden_dim),
                nn.SiLU(),
            )

            self.D0 = nn.Linear(hidden_dim, 1)
            self.A = nn.Linear(hidden_dim, 4)
            self.B = nn.Linear(hidden_dim, 4)
            self.softplus = nn.Softplus()
            self.p = nn.Linear(hidden_dim, 1)

        def forward(self, ctx, d):
            z = self.encoder(ctx)

            D0 = self.softplus(self.D0(z))
            A = self.A(z)
            B = self.softplus(self.B(z))
            p = 2.0 + 10.0 * torch.sigmoid(self.p(z))
            x = d - D0


            out = (
                A[:, 0:1] * torch.exp(-B[:, 0:1] * x**2) +
                A[:, 1:2] * torch.exp(-B[:, 1:2] * torch.pow(torch.abs(x), p)) +
                A[:, 2:3] * torch.exp(-B[:, 2:3] * x**4) +
                A[:, 3:4] * torch.exp(-B[:, 3:4] * x**5)
            )

            return out

    MODEL = ExponentialAsymptoticNet(context_dim=len(idx_ctx))
    MODEL.load_state_dict(checkpoint["model_state"])
    MODEL.to(DEVICE)
    MODEL.eval()

except Exception as e:
    raise RuntimeError(f"Error cargando el modelo ML: {e}")

# =====================================================================
# 1. FEATURE ENGINEERING
# =====================================================================
def add_features(df):
    df = df.copy()
    df['r_sum'] = df['rA'] + df['rB']
    df['l_tot'] = df['lA'] + df['lB']
    df['steric'] = (df['rA']**2 * df['covA']) + (df['rB']**2 * df['covB'])
    df['soft_A'] = df['lA'] / df['rA']
    df['soft_B'] = df['lB'] / df['rB']
    df['vol_A'] = df['rA']**3
    df['vol_B'] = df['rB']**3
    df['cov_mean'] = (df['covA'] + df['covB']) / 2
    df['geometric_mean_radius'] = np.sqrt(df['rA'] * df['rB'])
    df['harmonic_mean_radius'] = 2 * df['rA'] * df['rB'] / (df['rA'] + df['rB'])
    return df

# =====================================================================
# 2. BUILD FEATURE VECTOR
# =====================================================================

def build_features(rA, rB, lA, lB, covA, covB, D):
    D = np.atleast_1d(D)

    df = pd.DataFrame({
        "rA": rA,
        "rB": rB,
        "lA": lA,
        "lB": lB,
        "covA": covA,
        "covB": covB,
        "D": D,
    })

    df = add_features(df)

    missing = set(FEATURES) - set(df.columns)
    if missing:
        raise RuntimeError(f"Missing required features: {missing}")

    return df[FEATURES].values

# =====================================================================
# 3. PREDICTION FUNCTION
# =====================================================================

def predict_F(rA, rB, lA, lB, covA, covB, D_vals):
    X = build_features(rA, rB, lA, lB, covA, covB, D_vals)

    ctx = torch.tensor(X[:, idx_ctx], dtype=torch.float32, device=DEVICE)
    d = torch.tensor(X[:, [idx_d]], dtype=torch.float32, device=DEVICE)

    with torch.no_grad():
        F_scaled = MODEL(ctx, d).cpu().numpy().ravel()

    return F_scaled * F_SCALE


# =====================================================================
# 4. READ DEFINITIONS.txt
# =====================================================================

if not os.path.isfile("DEFINITIONS.txt"):
    raise RuntimeError("No se encontró DEFINITIONS.txt")

NNN = 0
radii = []

with open("DEFINITIONS.txt") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    line = line.strip()
    if line == "! number of particles":
        NNN = int(lines[i+1].split()[0])
    elif line == "!particle semiaxis x y z in nm":
        j = i + 1
        while j < len(lines) and lines[j].strip() and not lines[j].startswith("!"):
            vals = [float(x) for x in lines[j].split()]
            radii.append(vals[0])
            j += 1
        break

# =====================================================================
# 5. GENERATE distances.dat USING crystalCF
# =====================================================================

with open("DEFINITIONS.txt", "r") as f:
    lines = f.readlines()

with open("DEFINITIONS.txt_clt", "w") as f:
    for line in lines:
        if not any(k in line for k in ["dumpcluster", "cutoffcluster", "cluster_same"]):
            f.write(line)
    f.write("dumpcluster 2\n")
    f.write("cutoffcluster 15.0\n")
    f.write("cluster_same 0\n")

os.system("mv DEFINITIONS.txt DEFINITIONS.txt_tmp")
os.system("mv DEFINITIONS.txt_clt DEFINITIONS.txt")
os.system("~/develop/crystalCF/crystalCF")
os.system("mv DEFINITIONS.txt_tmp DEFINITIONS.txt")

# =====================================================================
# 6. READ distances.dat
# =====================================================================

if not os.path.isfile("distances.dat"):
    raise RuntimeError("No se encontró distances.dat")

with open("distances.dat") as f:
    lines = [l.strip() for l in f if l.strip()]

numdists = int(lines[0])

dists   = np.zeros(numdists)
weights = np.zeros(numdists)
partA   = np.zeros(numdists, dtype=int)
partB   = np.zeros(numdists, dtype=int)

idx = 1
for i in range(numdists):
    dists[i]   = float(lines[idx]); idx += 1
    weights[i] = float(lines[idx]); idx += 1
    partA[i]   = int(lines[idx]);   idx += 1
    partB[i]   = int(lines[idx]);   idx += 1

# =====================================================================
# 8. PAIRWISE FREE ENERGY
# =====================================================================

FIXED_LIGAND_LENGTH = 12.0
FIXED_COVERAGE = 5.85

Fpair = 0.0

print("Distance | F_ML | Weight")

for i in range(numdists):
    rA = radii[partA[i]-1]
    rB = radii[partB[i]-1]
    d  = dists[i]
    w  = weights[i]

    F = predict_F(
        rA, rB,
        FIXED_LIGAND_LENGTH, FIXED_LIGAND_LENGTH,
        FIXED_COVERAGE, FIXED_COVERAGE,
        d
    )[0]

    print(f"{d:8.4f} | {F:10.4f} | {w:6.3f}")
    Fpair += w * F/2

print("\n===================================")
print(f"F pares = {Fpair:.6f}")
print("===================================")

with open("F_pairwise.dat", "w") as f:
    f.write(f"1 {Fpair}\n")
