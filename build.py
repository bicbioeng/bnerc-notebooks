"""Builds the BNERC Colab notebook templates. Stdlib only.

The site copies the selected records to the clipboard and opens one of these
notebooks in Colab; the user pastes them into RECORDS and runs all cells.
No per-click GitHub writes, no tokens in the browser (same model as
bicbioeng/gedb-notebooks).
"""
import json, pathlib

def md(*lines): return {"cell_type": "markdown", "metadata": {}, "source": "\n".join(lines)}
def code(src, form=False):
    c = {"cell_type": "code", "metadata": {"cellView": "form"} if form else {}, "source": src.strip("\n"), "outputs": [], "execution_count": None}
    return c
def nb(cells):
    return {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python"}, "colab": {"provenance": []}}, "nbformat": 4, "nbformat_minor": 5}

INTRO = "Opened from **BNERC**, the Bionitrogen Fixation Research Database (BicBioEng lab, University of South Dakota)."

PARAMS = '''
#@title Records to analyze
#@markdown Paste what BNERC copied for you: `accession:source` pairs separated by commas. The source can be left out; it's inferred from the accession.
RECORDS = "GCA_000009705:cyanobacteria_ena, GCA_000020025:cyanobacteria_ena"  #@param {type:"string"}
#@markdown NCBI asks for a contact email with each download request.
NCBI_EMAIL = "readstech@gmail.com"  #@param {type:"string"}

def _parse(text):
    out = []
    for part in text.replace("\\n", ",").split(","):
        part = part.strip()
        if not part:
            continue
        acc, _, src = part.partition(":")
        acc, src = acc.strip(), src.strip()
        if not src:
            src = "cyanobacteria_ncbi" if acc.upper().startswith(("GCA_", "GCF_")) else "cyanorak"
        out.append((acc, src))
    return out

records = _parse(RECORDS)
accessions = [a for a, _ in records]
assert records, "Paste at least one record into RECORDS above."
print(f"{len(records)} record(s):", ", ".join(f"{a} ({s})" for a, s in records))
'''

FETCH = '''
#@title Download genomes as FASTA (.fna)
import glob, os, shutil, subprocess, urllib.request, zipfile
from Bio import Entrez
Entrez.email = NCBI_EMAIL

def _datasets_cli():
    exe = shutil.which("datasets") or os.path.abspath("bin/datasets")
    if not os.path.exists(exe):
        os.makedirs("bin", exist_ok=True)
        urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets", exe)
        os.chmod(exe, 0o755)
    return exe

def fetch(acc, source):
    dest = f"{acc}.fna"
    if os.path.exists(dest):
        return print(f"✅ {dest} already downloaded")
    try:
        if source in ("cyanobacteria_ncbi", "cyanobacteria_gtdb"):
            # Assembly accessions (GCA_/GCF_): NCBI Datasets CLI
            subprocess.run([_datasets_cli(), "download", "genome", "accession", acc, "--include", "genome", "--filename", f"{acc}.zip"], check=True)
            with zipfile.ZipFile(f"{acc}.zip") as z:
                z.extractall(f"genome_{acc}")
            hits = glob.glob(f"genome_{acc}/ncbi_dataset/data/*/*.fna")
            if not hits:
                return print(f"❌ No .fna in the NCBI package for {acc}")
            shutil.move(hits[0], dest)
        elif source == "cyanobacteria_ena":
            urllib.request.urlretrieve(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{acc}?download=true", dest)
        else:
            # cyanorak / cyanobacteria_ncbi_nucelotide: NCBI Nucleotide via Entrez
            with Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text") as h:
                open(dest, "w").write(h.read())
        print(f"✅ {dest} ({os.path.getsize(dest) // 1024} KB)")
    except Exception as e:
        print(f"❌ {acc} ({source}): {e}")

for acc, src in records:
    fetch(acc, src)
missing = [a for a in accessions if not os.path.exists(f"{a}.fna") or os.path.getsize(f"{a}.fna") == 0]
if missing:
    print("⚠️ Not downloaded:", ", ".join(missing))
'''

comparative = nb([
    md("# Comparative genomics with BLAST", INTRO, "",
       "Downloads each selected genome, builds one BLAST database from all of them, then searches it with the first genome as the query and plots the top hits.",
       "", "**How to run:** paste the records BNERC copied into `RECORDS` below, then choose **Runtime → Run all**."),
    code(PARAMS, form=True),
    md("## 1. Install tools"),
    code("!pip -q install biopython pandas matplotlib\n!apt-get -qq install -y ncbi-blast+ > /dev/null && blastn -version | head -1"),
    md("## 2. Download sequences"),
    code(FETCH, form=True),
    md("## 3. Build one BLAST database"),
    code('''
with open("all_sequences.fna", "w") as out:
    for acc in accessions:
        if os.path.exists(f"{acc}.fna"):
            out.write(open(f"{acc}.fna").read())
!makeblastdb -in all_sequences.fna -dbtype nucl -out combined_db
'''),
    md("## 4. BLAST the first genome against all of them"),
    code('''
query = accessions[0]
print("Query:", query)
!blastn -query {query}.fna -db combined_db -out blast_results.txt -outfmt "6 qseqid sseqid length pident bitscore evalue" -max_target_seqs 50
'''),
    md("## 5. Top hits"),
    code('''
import pandas as pd
import matplotlib.pyplot as plt
cols = ["qseqid", "sseqid", "length", "pident", "bitscore", "evalue"]
df = pd.read_csv("blast_results.txt", sep="\\t", names=cols).sort_values("bitscore", ascending=False).head(10)
display(df)
plt.figure(figsize=(10, 5))
plt.barh(df["sseqid"], df["bitscore"], color="#0E8C8C")
plt.title(f"Top BLAST hits for {query}")
plt.xlabel("Bit score")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()
'''),
])

pangenome = nb([
    md("# Pangenome analysis with PPanGGOLiN", INTRO, "",
       "Downloads each selected genome, then annotates, clusters and partitions them into core, shell and cloud genes with PPanGGOLiN, and draws the U-curve and tile plot. Use at least 3 genomes.",
       "", "**How to run:** paste the records BNERC copied into `RECORDS` below, then choose **Runtime → Run all**. Installing PPanGGOLiN takes a few minutes."),
    code(PARAMS, form=True),
    md("## 1. Install PPanGGOLiN (via Miniconda)"),
    code('''
!wget -nc -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
!bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/local > /dev/null
!conda install -y -q -c conda-forge -c bioconda ppanggolin > /dev/null
!pip -q install biopython
!ppanggolin --version
'''),
    md("## 2. Download sequences"),
    code(FETCH, form=True),
    md("## 3. Genome list for PPanGGOLiN"),
    code('''
os.makedirs("genomes_fna", exist_ok=True)
lines = []
for acc in accessions:
    if os.path.exists(f"{acc}.fna"):
        shutil.copy(f"{acc}.fna", f"genomes_fna/{acc}.fna")
        lines.append(f"{acc}\\t{os.path.abspath(f'genomes_fna/{acc}.fna')}")
open("genomes_list.tsv", "w").write("\\n".join(lines) + "\\n")
print(len(lines), "genomes listed")
assert len(lines) >= 2, "Pangenome analysis needs at least 2 downloaded genomes (3+ recommended)."
'''),
    md("## 4. Annotate, cluster, build the graph, partition"),
    code('''
!ppanggolin annotate --fasta genomes_list.tsv --output annotated --force
!ppanggolin cluster -p annotated/pangenome.h5 --identity 0.5 --coverage 0.8 --force
!ppanggolin graph -p annotated/pangenome.h5 --force
!ppanggolin partition -p annotated/pangenome.h5 --force
!ppanggolin info -p annotated/pangenome.h5 --content
'''),
    md("## 5. Export and plots"),
    code('''
os.makedirs("output", exist_ok=True)
!ppanggolin write_pangenome -p annotated/pangenome.h5 --csv --output output/presence_absence --force
!ppanggolin write_pangenome -p annotated/pangenome.h5 --light_gexf --output output/pangenome_graph --force
!ppanggolin draw -p annotated/pangenome.h5 --ucurve --output output/fig_ucurve --force
!ppanggolin draw -p annotated/pangenome.h5 --tile_plot --nocloud --output output/fig_tile --force
from IPython.display import HTML, display
for f in glob.glob("output/fig_*/*.html"):
    print(f)
    display(HTML(open(f).read()))
'''),
])

GEO_PARAMS = '''
#@title GEO series to analyze
#@markdown Paste the GEO accessions BNERC copied for you (e.g. `GSE12345, GSE67890`).
RECORDS = "GSE10000"  #@param {type:"string"}
series = [s.split(":")[0].strip() for s in RECORDS.replace("\\n", ",").split(",") if s.strip()]
assert series, "Paste at least one GEO accession into RECORDS above."
print(series)
'''
geo = nb([
    md("# Expression data from NCBI GEO", INTRO, "",
       "Links each series to GEO2R for interactive differential expression, and downloads its series matrix into a table.",
       "", "**How to run:** paste the GEO accessions BNERC copied into `RECORDS`, then choose **Runtime → Run all**."),
    code(GEO_PARAMS, form=True),
    md("## 1. Open in GEO2R"),
    code('''
from IPython.display import Markdown, display
display(Markdown("\\n".join(f"- [{s} in GEO2R](https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc={s})" for s in series)))
'''),
    md("## 2. Download series matrices"),
    code('''
import gzip, io, urllib.request
import pandas as pd
tables = {}
for s in series:
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{s[:-3]}nnn/{s}/matrix/{s}_series_matrix.txt.gz"
    try:
        raw = gzip.decompress(urllib.request.urlopen(url).read()).decode("utf-8", "replace")
        body = raw.split("!series_matrix_table_begin")[1].split("!series_matrix_table_end")[0]
        tables[s] = pd.read_csv(io.StringIO(body), sep="\\t", index_col=0)
        print(f"✅ {s}: {tables[s].shape[0]} rows × {tables[s].shape[1]} samples")
    except Exception as e:
        print(f"❌ {s}: {e}")
for s, t in tables.items():
    display(t.head())
'''),
])

out = pathlib.Path(__file__).parent
for name, n in [("comparative-genomics", comparative), ("pangenome", pangenome), ("geo-expression", geo)]:
    (out / f"{name}.ipynb").write_text(json.dumps(n, indent=1, ensure_ascii=False) + "\n")
    print("wrote", name)
