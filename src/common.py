import os
import re

# Project paths
ROOT_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(ROOT_DIR, 'data')
ORIG_DIR = os.path.join(DATA_DIR, 'original_data')
PROC_DIR = os.path.join(DATA_DIR, 'processed_data')

PATHS = {
    'NETWORK': os.path.join(ORIG_DIR, 'protein_network_data.txt'),
    'PROTEINS': os.path.join(ORIG_DIR, 'unique_STRING_proteins.txt'),
    'GAF': os.path.join(ORIG_DIR, 'sgd.gaf'),
    'FASTA': os.path.join(ORIG_DIR, 'idmapping_2025_10_27.fasta'),
    'PROCESSED_PKL': os.path.join(PROC_DIR, 'processed_protein_data.pkl'),
    'PROCESSED_PKL_ESM': os.path.join(PROC_DIR, 'processed_protein_data_with_esm.pkl'),
    'STATS_TXT': os.path.join(PROC_DIR, 'data_statistics.txt'),
}

def ensure_dirs():
    os.makedirs(PROC_DIR, exist_ok=True)

# Shared patterns (yeast systematic name regex: YAL001C, YAL016C-A, etc.)
SYS_NAME_PAT = re.compile(r'^Y[A-P][LR][0-9]{3}[CW](?:-[A-Z])?$')
TAXON_PREFIX = '4932.'

def parse_fasta(filename):
    """Parse FASTA into {id: sequence} using header id components."""
    sequences = {}
    current_id = None
    current_seq = []
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                header_parts = line[1:].split('|')
                if len(header_parts) > 1:
                    current_id = header_parts[1]
                else:
                    current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        sequences[current_id] = ''.join(current_seq)
    return sequences

def parse_fasta_gene_symbol_to_seq(path):
    """Map gene symbol (GN=...) to sequence from UniProt FASTA."""
    gene_to_seq = {}
    current_gene = None
    current_seq = []
    gn_pat = re.compile(r'\bGN=([A-Za-z0-9_\-\.]+)\b')
    with open(path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.rstrip('\n')
            if line.startswith('>'):
                if current_gene and current_seq:
                    seq = ''.join(current_seq).replace(' ', '')
                    if seq:
                        gene_to_seq.setdefault(current_gene.upper(), seq)
                current_gene = None
                current_seq = []
                m = gn_pat.search(line)
                if m:
                    current_gene = m.group(1)
            else:
                if line and not line.startswith('!'):
                    current_seq.append(line.strip())
    if current_gene and current_seq:
        seq = ''.join(current_seq).replace(' ', '')
        if seq:
            gene_to_seq.setdefault(current_gene.upper(), seq)
    return gene_to_seq

def build_string_to_gene_symbol_from_gaf(gaf_path, target_string_ids=None):
    """
    Build mapping: STRING ID (4932.Y*) -> gene symbol from GAF.
    Uses system name tokens in synonyms to derive STRING IDs.
    """
    string_to_symbol = {}
    with open(gaf_path, 'r', encoding='utf-8') as f:
        for raw in f:
            if raw.startswith('!'):
                continue
            parts = raw.rstrip('\n').split('\t')
            if len(parts) < 11:
                continue
            db_obj_symbol = parts[2].strip()
            synonyms = parts[10].strip()
            if not synonyms:
                continue
            for tok in synonyms.split('|'):
                tok = tok.strip()
                if SYS_NAME_PAT.match(tok):
                    sid = TAXON_PREFIX + tok
                    if (target_string_ids is None) or (sid in target_string_ids):
                        string_to_symbol.setdefault(sid, db_obj_symbol)
    return string_to_symbol

def load_esm_model():
    """Load ESM-2 model, supporting fair-esm installs."""
    try:
        import esm  # fair-esm
        if hasattr(esm, "pretrained"):
            return esm.pretrained.esm2_t33_650M_UR50D()
        from esm import pretrained
        return pretrained.esm2_t33_650M_UR50D()
    except Exception as e:
        raise RuntimeError(
            "Failed to load ESM model. Please run:\n"
            "  pip uninstall -y esm\n"
            "  pip install -U fair-esm\n"
            "and restart your terminal before retrying. Original error: {}".format(e)
        )