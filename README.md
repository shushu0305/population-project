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
git clone https://github.com/shushu0305/population-project.git
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### 3. Launch the web interface

```bash
cd src
streamlit run app.py
```

### 4. Or use the CLI

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


---

