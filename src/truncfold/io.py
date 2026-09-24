from pathlib import Path
AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")

def _validate(seq, path):
    if not seq:
        raise ValueError(f"Empty sequence in {path}")
    bad = sorted(set(seq) - AA)
    if bad:
        raise ValueError(f"Unsupported amino-acid symbols in {path}: {bad}")
    return seq

def read_fasta(path):
    path = Path(path)
    records, header, parts = [], None, []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line: continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, _validate("".join(parts), path)))
            header, parts = (line[1:].strip() or path.stem), []
        else:
            parts.append(line.replace(" ", "").upper())
    if header is not None:
        records.append((header, _validate("".join(parts), path)))
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    return records

def read_single_fasta(path):
    recs = read_fasta(path)
    if len(recs) != 1:
        raise ValueError(f"Expected one FASTA record in {path}, found {len(recs)}")
    return recs[0]

def write_fasta(records, path, width=80):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8') as f:
        for name,seq in records:
            f.write(f">{name}\n")
            for i in range(0,len(seq),width): f.write(seq[i:i+width]+"\n")
