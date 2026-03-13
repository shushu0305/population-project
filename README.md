# Haplogroup Discover 🧬

Given an mtDNA or Y-DNA haplogroup, search the AADR ancient genomics dataset for the earliest matching samples, infer a candidate geographic origin, and find notable historical individuals who share the same lineage.

---

## Features

- **Ancient sample search** — filter AADR samples by haplogroup clade (exact + upstream + downstream)
- **Origin inference** — rank candidate countries based on the earliest matching samples
- **VIP matching** — link haplogroups to notable historical figures
- **Streamlit UI** — interactive web interface with Wikipedia integration
- **CLI** — scriptable command-line interface for batch use

---

## Quick start

### 1. Clone the repository

```bash
git clone https://github.com/yourname/haplogroup-discover.git
cd haplogroup-discover
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### 3. Prepare the AADR data file

Tree files and the VIP file are already included in the repository. **You only need to download the AADR file manually:**

1. Go to the [Reich Lab AADR page](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data)
2. Download the annotation file (`.xlsx` or `.tsv`)
3. Place it in the `data/aadr/` directory

```bash
# To use a custom path, copy the environment variable template
cp .env.example .env
# Edit .env and set AADR_PATH to your file location
```

### 4. Launch the web interface

```bash
cd src
streamlit run app.py
```

### 5. Or use the CLI

```bash
cd src

# Basic usage (paths are read automatically from config.py / .env)
python main.py --system mt --target U5b2c

# With explicit paths and saved output
python main.py \
  --system y \
  --target R-M269 \
  --aadr ../data/aadr/AADR\ Annotations\ 2025.xlsx \
  --early-n 10 \
  --save
```

---

## Project structure

```
haplogroup-discover/
├── README.md
├── requirements.txt
├── .env.example
├── .gitignore
│
├── src/
│   ├── app.py              # Streamlit web entry point
│   ├── main.py             # CLI entry point
│   ├── config.py           # Path configuration (supports env var overrides)
│   ├── service.py          # Analysis pipeline orchestration
│   ├── analysis_origin.py  # Origin-focused analyzer
│   ├── analysis_vip.py     # VIP matching analyzer
│   ├── tree_parser.py      # Phylogenetic tree parser
│   ├── aadr_parser.py      # AADR annotation loader
│   ├── vip_parser.py       # VIP haplogroup list loader
│   ├── vip_matcher.py      # Haplogroup relationship matching
│   ├── y_mapper.py         # Y-DNA label resolution utilities
│   └── export_y_aliases.py # Y-DNA alias export tool
│
└── data/
    ├── aadr/               # ⚠ Not included — download separately (see above)
    ├── trees/              # Phylogenetic tree files (included)
    ├── vip/                # VIP haplogroup list (included)
    └── y_aliases.tsv       # Y-DNA label mapping table (included)
```

---

## Data sources

| File | Source | Included in repo |
|------|--------|-----------------|
| AADR annotation | [Reich Lab](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data) | ❌ Download separately |
| mtDNA tree | [PhyloTree Build 17](https://www.phylotree.org/) | ✅ |
| Y-DNA tree | [ISOGG 2016](https://isogg.org/) | ✅ |
| VIP haplogroup list | Custom-curated | ✅ |

---

## VIP file format

The VIP Excel file should have one sheet per system (`mtDNA`, `Y`), with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `Name` / `Individual` / `VIP` | ✅ | Person's name |
| `Haplogroup` / `mtDNA` / `HG` | ✅ | Haplogroup label |
| `System` | optional | `mt` or `y` (inferred from sheet name if absent) |
| `Source` | optional | Citation or reference |
| `Note` | optional | Free-text note |

---

## License

MIT
