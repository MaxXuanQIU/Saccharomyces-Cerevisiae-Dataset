# esm_feature_extractor.py
import os
import torch
import pickle
import numpy as np

from .common import (
    PATHS, ensure_dirs,
    load_esm_model,
    parse_fasta_gene_symbol_to_seq,
    build_string_to_gene_symbol_from_gaf
)

def extract_esm_features():
    ensure_dirs()

    # Read processed data
    if not os.path.exists(PATHS['PROCESSED_PKL']):
        raise FileNotFoundError(f"File not found: {PATHS['PROCESSED_PKL']}, please run data_processor.py first")

    with open(PATHS['PROCESSED_PKL'], 'rb') as f:
        data = pickle.load(f)

    proteins = list(data.get('proteins', []))
    if not proteins:
        print("No proteins available for feature extraction (proteins list is empty).")
        return

    # Build STRING -> gene symbol mapping (for finding sequences from FASTA)
    print("Building STRING->gene symbol mapping...")
    string_to_symbol = build_string_to_gene_symbol_from_gaf(PATHS['GAF'], target_string_ids=set(proteins))
    print(f"  Number of proteins mapped to gene symbols: {len(string_to_symbol)} / {len(proteins)}")

    # Build gene symbol -> sequence mapping
    print("Parsing FASTA (gene symbol -> sequence)...")
    if not os.path.exists(PATHS['FASTA']):
        raise FileNotFoundError(f"File not found: {PATHS['FASTA']}, please check FASTA path.")
    gene_to_seq = parse_fasta_gene_symbol_to_seq(PATHS['FASTA'])
    print(f"  Number of available gene symbols in FASTA: {len(gene_to_seq)}")

    # Combine into STRING -> sequence mapping
    protein_to_sequence = {}
    missing_seq = []
    for sid in proteins:
        sym = string_to_symbol.get(sid)
        if sym:
            seq = gene_to_seq.get(sym.upper())
            if seq:
                protein_to_sequence[sid] = seq
            else:
                missing_seq.append(sid)
        else:
            missing_seq.append(sid)

    print(f"Number of proteins with sequences found: {len(protein_to_sequence)} / {len(proteins)}")
    if missing_seq:
        print(f"Warning: {len(missing_seq)} proteins have no sequences found, will be skipped.")

    # Load ESM model
    print("Loading ESM model...")
    model, alphabet = load_esm_model()
    batch_converter = alphabet.get_batch_converter()
    model.eval()

    # Device and precision settings
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    use_amp = device.type == "cuda"

    # Extract features
    protein_features = {}
    seq_items = list(protein_to_sequence.items())
    if not seq_items:
        print("No available sequences, terminating. Please check if FASTA is complete and GAF to STRING ID mapping is correct.")
        return

    # Batch size can be adjusted; 8~16 for GPU depending on memory
    batch_size = 8 if device.type == "cuda" else 2

    for i in range(0, len(seq_items), batch_size):
        batch = seq_items[i:i + batch_size]
        batch_labels, batch_strs, batch_tokens = batch_converter([(k, v) for k, v in batch])
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            if use_amp:
                with torch.cuda.amp.autocast(dtype=torch.float16):
                    results = model(batch_tokens, repr_layers=[33], return_contacts=False)
            else:
                results = model(batch_tokens, repr_layers=[33], return_contacts=False)

            token_reps = results["representations"][33]

        for j, (sid, _seq) in enumerate(batch):
            # Remove start/end special tokens
            rep = token_reps[j, 1:-1].mean(0).detach().cpu().float().numpy()
            protein_features[sid] = rep
        print(f"Progress: {min(i + batch_size, len(seq_items))}/{len(seq_items)}")

    # Merge and save
    # data['protein_features'] = protein_features
    # with open(PATHS['PROCESSED_PKL_ESM'], 'wb') as f:
    #     pickle.dump(data, f)

    # Alignment: uniform filtering, drop samples missing sequences/features to ensure consistency ===
    print("Aligning data structure to match samples with extracted ESM features...")
    proteins_all = list(data.get('proteins', []))
    # Remove dependency on protein_to_idx
    labels_all = data.get('labels', None)
    edges_all = data.get('network_edges', [])

    keep_set = set(protein_features.keys())
    dropped = [p for p in proteins_all if p not in keep_set]
    print(f"Will drop samples without sequences/features: {len(dropped)} samples")

    # 1) Regenerate proteins (maintain original order subsequence)
    proteins_new = [p for p in proteins_all if p in keep_set]

    # 2) Filter label rows: use old index list selected_rows
    labels_new = labels_all
    if isinstance(labels_all, np.ndarray) and labels_all.ndim == 2:
        selected_rows = np.array([i for i, p in enumerate(proteins_all) if p in keep_set], dtype=int)
        labels_new = labels_all[selected_rows, :]
        assert labels_new.shape[0] == len(proteins_new), "Number of filtered label rows inconsistent with protein count"

    # 3) Filter network_edges (maintain original logic)
    edges_new = []
    if isinstance(edges_all, list) and len(edges_all) > 0 and isinstance(edges_all[0], (list, tuple)) and len(edges_all[0]) >= 2:
        e0 = edges_all[0]
        if isinstance(e0[0], str):
            keep = keep_set
            edges_new = [(u, v) for (u, v, *rest) in edges_all if (u in keep and v in keep)]
        else:
            # Old index -> new index mapping (temporary, not persistent)
            selected_rows = [i for i, p in enumerate(proteins_all) if p in keep_set]
            old_to_new_idx = {old_i: new_i for new_i, old_i in enumerate(selected_rows)}
            for e in edges_all:
                u, v = e[0], e[1]
                if u in old_to_new_idx and v in old_to_new_idx:
                    edges_new.append((old_to_new_idx[u], old_to_new_idx[v]))
    else:
        edges_new = edges_all # Keep as is if structure cannot be determined

    # 4) Filter protein_features
    protein_features_new = {p: protein_features[p] for p in proteins_new}

    # Write back
    data['proteins'] = proteins_new
    data['labels'] = labels_new
    data['network_edges'] = edges_new
    data['protein_features'] = protein_features_new

    with open(PATHS['PROCESSED_PKL_ESM'], 'wb') as f:
        pickle.dump(data, f)

    # Print aligned data size
    print(f"Number of proteins after alignment: {len(data['proteins'])}")
    if isinstance(data['labels'], np.ndarray):
        print(f"Label matrix shape after alignment: {data['labels'].shape}")
    
    print(f"ESM feature extraction completed and aligned! Saved to {PATHS['PROCESSED_PKL_ESM']}")

if __name__ == "__main__":
    extract_esm_features()