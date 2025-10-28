import pandas as pd
import numpy as np
from collections import defaultdict
import pickle
import re

from .common import (
    PATHS, ensure_dirs, parse_fasta, SYS_NAME_PAT, TAXON_PREFIX
)

def main():
    print("Starting protein data processing...")
    ensure_dirs()

    # 1. Read protein network data
    print("1. Reading network data...")
    network_df = pd.read_csv(PATHS['NETWORK'], sep=' ')
    print(f"   Network contains {len(network_df)} interactions")
    print(f"   Unique protein count: {len(set(network_df['protein1']).union(set(network_df['protein2'])))}")
    
    # 2. Read unique STRING protein list
    print("2. Reading protein list...")
    proteins_df = pd.read_csv(PATHS['PROTEINS'], sep='\t')
    string_proteins = set(proteins_df['#string_protein_id'])
    print(f"   Unique proteins in STRING database: {len(string_proteins)}")
    
    # 3. Read GO annotation data
    print("3. Reading GO annotation data...")
    go_annotations = []
    with open(PATHS['GAF'], 'r', encoding='utf-8') as f:
        for line in f:
            if not line.startswith('!'):  # Skip comment lines
                parts = line.strip().split('\t')
                if len(parts) >= 17:  # Ensure enough columns
                    go_annotations.append(parts)
    
    # Create GO annotation DataFrame
    go_columns = [
        'DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 
        'GO_ID', 'DB_Reference', 'Evidence_Code', 'With_From', 
        'Aspect', 'DB_Object_Name', 'DB_Object_Synonym', 
        'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 
        'Annotation_Extension', 'Gene_Product_Form_ID'
    ]
    go_df = pd.DataFrame(go_annotations, columns=go_columns)
    print(f"   Original GO annotation count: {len(go_df)}")
    
    # 4. Read FASTA sequence data (backup for future ESM features)
    print("4. Reading FASTA sequence data...")
    sequences = parse_fasta(PATHS['FASTA'])
    print(f"   Parsed {len(sequences)} protein sequences")
    
    # 5. Data cleaning and integration
    print("5. Data cleaning and integration...")
    
    # Filter high-quality experimental evidence
    experimental_evidence = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'}
    high_quality_go = go_df[go_df['Evidence_Code'].isin(experimental_evidence)].copy()
    print(f"   High-quality experimental evidence annotations: {len(high_quality_go)}")
    
    # Select top 50 most common GO terms
    top_go_terms = high_quality_go['GO_ID'].value_counts().head(50).index.tolist()
    print(f"   Selected top 50 most common GO terms")
    
    # 5.1 Map GO annotations to STRING protein IDs (4932.YXXXXXX)
    print("   Mapping GO annotations to STRING protein IDs...")
    
    def row_to_string_ids(row):
        ids = set()
        # Synonyms usually contain systematic names (e.g., FUN24|TSV115|YAL001C)
        syn = row.get('DB_Object_Synonym', '')
        if isinstance(syn, str) and syn:
            for token in syn.split('|'):
                tok = token.strip()
                if SYS_NAME_PAT.match(tok):
                    ids.add(TAXON_PREFIX + tok)
        sym = row.get('DB_Object_Symbol', '')
        if isinstance(sym, str) and SYS_NAME_PAT.match(sym.strip()):
            ids.add(TAXON_PREFIX + sym.strip())
        return ids
    
    # Only keep annotations from top_go_terms and build STRING ID -> GO set mapping
    mapped_rows = 0
    protein_to_go = defaultdict(set)  # Key: STRING protein ID (e.g., 4932.YAL001C)
    for _, row in high_quality_go.iterrows():
        go_id = row['GO_ID']
        if go_id not in top_go_terms:
            continue
        string_ids = row_to_string_ids(row)
        if string_ids:
            mapped_rows += 1
            for sid in string_ids:
                protein_to_go[sid].add(go_id)
    print(f"   GO annotation rows mapped to STRING via systematic names: {mapped_rows}")
    print(f"   STRING proteins with GO labels: {len(protein_to_go)}")
    
    # Collect all proteins appearing in the network (STRING IDs)
    all_network_proteins = set(network_df['protein1']).union(set(network_df['protein2']))
    
    # 6. Create final dataset
    print("6. Creating final dataset...")

    # Find proteins with both network connections and GO annotations (intersection by STRING ID)
    valid_proteins = sorted(all_network_proteins.intersection(protein_to_go.keys()))
    print(f"   Valid protein count (with network and GO annotations): {len(valid_proteins)}")
    
    # Create protein feature matrix placeholder (will be replaced with ESM features later)
    protein_features = {}
    for protein in valid_proteins:
        # Create random features for now, will be replaced with real ESM features later
        protein_features[protein] = np.random.randn(1280).astype(np.float32)  # Typical ESM-2 dimension
    
    # Create label matrix
    labels = np.zeros((len(valid_proteins), len(top_go_terms)), dtype=np.float32)
    for p_idx, protein in enumerate(valid_proteins):
        for go_term in protein_to_go.get(protein, ()):
            if go_term in top_go_terms:
                g_idx = top_go_terms.index(go_term)
                labels[p_idx, g_idx] = 1.0
    
    print(f"   Label matrix shape: {labels.shape}")
    print(f"   Positive sample count: {np.sum(labels)}")

    # 7. Save processed data
    print("7. Saving data...")
    output_data = {
        'proteins': valid_proteins,
        'protein_features': protein_features,  # Placeholder, will be replaced with real ESM features
        'go_terms': top_go_terms,
        'labels': labels,
        'network_edges': network_df[['protein1', 'protein2']].values.tolist(),
        #'protein_to_go': dict(protein_to_go)
    }
    with open(PATHS['PROCESSED_PKL'], 'wb') as f:
        pickle.dump(output_data, f)
    
    # Save statistics
    with open(PATHS['STATS_TXT'], 'w', encoding='utf-8') as f:
        f.write("=== Protein Data Preprocessing Statistics ===\n\n")
        f.write(f"Network interactions count: {len(network_df)}\n")
        f.write(f"STRING unique proteins: {len(string_proteins)}\n")
        f.write(f"High-quality GO annotations: {len(high_quality_go)}\n")
        f.write(f"Selected GO terms count: {len(top_go_terms)}\n")
        f.write(f"Valid proteins count: {len(valid_proteins)}\n")
        f.write(f"Label matrix shape: {labels.shape}\n")
        f.write(f"Total positive samples: {np.sum(labels)}\n")
        avg = (np.sum(labels) / len(valid_proteins)) if len(valid_proteins) > 0 else 0.0
        f.write(f"Average GO terms per protein: {avg:.2f}\n\n")
        f.write("Top 10 most common GO terms:\n")
        for go_term in top_go_terms[:10]:
            go_term_index = top_go_terms.index(go_term) if len(valid_proteins) > 0 else 0
            count = int(np.sum(labels[:, go_term_index])) if len(valid_proteins) > 0 else 0
            f.write(f"  {go_term}: {count} proteins\n")
    
    print("\n=== Data processing completed ===")
    print("Output files:")
    print(f"  - {PATHS['PROCESSED_PKL']}")
    print(f"  - {PATHS['STATS_TXT']}")
    print("\nNext step: Run ESM feature extraction script")

if __name__ == "__main__":
    main()