### Data files

- `protein_network_data.txt` — Protein interaction network from STRING.  
    Source: [STRING Database - Saccharomyces cerevisiae](https://cn.string-db.org/cgi/download?species=Saccharomyces+cerevisiae)
- `list_of_STRING_proteins.txt` — Protein list exported from STRING.  
    Source: [STRING Database - Saccharomyces cerevisiae](https://cn.string-db.org/cgi/download?species=Saccharomyces+cerevisiae)
- `unique_STRING_proteins.txt` — Unique protein IDs (first column extracted from `list_of_STRING_proteins.txt`).  
- `sgd.gaf` — Gene Ontology (GO) annotation file from Gene Ontology.  
    Source: [Gene Ontology Downloads](https://current.geneontology.org/products/pages/downloads.html)
- `idmapping_2025_10_27.fasta` — UniProt mapping results (used `unique_STRING_proteins.txt`; 6,156 IDs mapped to 6,156 results; 444 IDs not mapped).  
    Source: [UniProt ID Mapping](https://www.uniprot.org/id-mapping)