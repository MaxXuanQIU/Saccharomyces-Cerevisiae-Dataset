# This file is used to inspect the content and structure of processed protein data, helping to understand the components of the dataset.
import pickle
import os
import numpy as np

# Prioritize reading data with ESM features, fallback if not available
pkl_candidates = [
    'data/processed_data/processed_protein_data_with_esm.pkl',
    'data/processed_data/processed_protein_data.pkl'
]

pkl_path = None
for p in pkl_candidates:
    if os.path.exists(p):
        pkl_path = p
        break

if pkl_path is None:
    raise FileNotFoundError("No available pkl data file found, please confirm if files exist in the processed_data directory.")

print(f"Loading data file: {pkl_path}")
with open(pkl_path, 'rb') as f:
    data = pickle.load(f)

# List all keys
print("Keys contained in data:", list(data.keys()))

def safe_len(x):
    try:
        return len(x)
    except Exception:
        return "N/A"

def shape_of(x):
    if isinstance(x, np.ndarray):
        return tuple(x.shape)
    return "N/A"

# Describe the length and shape information of each key
print("\n=== Scale Information for Each Key ===")
for k, v in data.items():
    info_len = safe_len(v)
    info_shape = shape_of(v)
    v_type = type(v).__name__
    extra = ""
    # More detailed description for common keys
    if k == 'protein_features':
        if isinstance(v, dict) and len(v) > 0:
            any_vec = next(iter(v.values()))
            if isinstance(any_vec, np.ndarray):
                extra = f", feature vector shape={any_vec.shape}, dtype={any_vec.dtype}"
        elif isinstance(v, np.ndarray):
            extra = f", shape={v.shape}, dtype={v.dtype}"
    elif k == 'labels' and isinstance(v, np.ndarray):
        extra = f", shape={v.shape}, dtype={v.dtype}, total positive samples={int(v.sum())}"
    print(f"- {k}: type={v_type}, length={info_len}, shape={info_shape}{extra}")

# Print the first data entry as an example (using the first protein as example)
print("\n=== Example: First Protein Entry ===")
proteins = data.get('proteins', [])
go_terms = data.get('go_terms', [])
labels = data.get('labels', None)
protein_features = data.get('protein_features', None)

if isinstance(proteins, list) and len(proteins) > 0:
    first_pid = proteins[0]
else:
    first_pid = None

if first_pid is None:
    print("Could not find example protein (proteins is empty).")
else:
    print(f"Protein ID: {first_pid}")

    # Feature example
    feat_vec = None
    if isinstance(protein_features, dict):
        feat_vec = protein_features.get(first_pid, None)

    if isinstance(feat_vec, np.ndarray):
        preview = np.array2string(feat_vec[:8], precision=4, separator=', ')
        print(f"Features: shape={feat_vec.shape}, dtype={feat_vec.dtype}, first 8 items={preview}")
    else:
        print("Features: Could not find feature vector for this protein or type not supported for preview.")

    # Label example - get labels for the first protein
    if isinstance(labels, np.ndarray) and labels.ndim == 2:
        # Get label vector for the first protein (row 0)
        first_protein_labels = labels[0, :]
        pos_idx = np.where(first_protein_labels > 0.5)[0].tolist()
        pos_go = [go_terms[i] for i in pos_idx] if isinstance(go_terms, list) and len(go_terms) > 0 else pos_idx
        print(f"Labels: positive GO count={len(pos_idx)}, first few GOs={pos_go[:10]}")
    else:
        print("Labels: Could not find labels for this protein or label matrix unavailable.")

    # Network edge example
    edges = data.get('network_edges', [])
    if isinstance(edges, list) and len(edges) > 0:
        print(f"Network edge examples (first 3): {edges[:3]}")
    else:
        print("Network edges: unavailable or empty.")