from __future__ import annotations

import re
from html import escape
from urllib.parse import quote

import pandas as pd
import requests
import streamlit as st
import plotly.express as px

from service import run_haplogroup_analysis

# ---------- Internal config ----------
AADR_PATH = "../data/aadr/AADR Annotations 2025.xlsx"
VIP_PATH = "../data/vip/VIPHaplogroups.xlsx"
MT_TREE_PATH = "../data/trees/mt_phyloTree_b17_Tree2.txt"
Y_TREE_PATH = "../data/trees/chrY_hGrpTree_isogg2016.txt"
EARLY_N = 5
PREVIEW_VIP_COUNT = 8

st.set_page_config(page_title="Haplogroup Discover", page_icon="🧬", layout="wide")

# ---------- Session state ----------
if "analysis_result" not in st.session_state:
    st.session_state.analysis_result = None
if "analysis_target" not in st.session_state:
    st.session_state.analysis_target = None
if "analysis_system" not in st.session_state:
    st.session_state.analysis_system = None
if "selected_vip_idx" not in st.session_state:
    st.session_state.selected_vip_idx = None
if "vip_show_all" not in st.session_state:
    st.session_state.vip_show_all = False

# ---------- Styles ----------
st.markdown(
    """
<style>
:root {
    --ink: #111827;
    --muted: #6b7280;
    --line: #e5e7eb;
    --paper: #f8fafc;
    --card: #ffffff;
    --gold: #b7791f;
    --gold-bg: #fff8eb;
    --blue: #2563eb;
    --blue-bg: #eff6ff;
    --green: #2f855a;
    --green-bg: #ecfdf5;
}
html, body, [class*="css"] {
    background: var(--paper);
    color: var(--ink);
}
#MainMenu, footer, header {visibility: hidden;}
.block-container {max-width: 1180px; padding-top: 2rem; padding-bottom: 3rem;}

.hero-title {font-size: 2.2rem; font-weight: 800; margin-bottom: 0.25rem;}
.hero-sub {color: var(--muted); margin-bottom: 1.5rem;}

.metric-card {
    border: 1px solid var(--line);
    border-radius: 16px;
    padding: 14px 16px;
    background: var(--card);
    min-height: 98px;
}
.metric-label {font-size: 12px; color: var(--muted); text-transform: uppercase; letter-spacing: .06em;}
.metric-value {font-size: 28px; font-weight: 700; margin-top: 8px; line-height: 1.1;}

.hero-card, .plain-card {
    border: 1px solid var(--line);
    border-radius: 18px;
    background: var(--card);
    padding: 20px;
    box-shadow: 0 1px 2px rgba(0,0,0,0.03);
}
.eyebrow {font-size: 12px; color: var(--muted); text-transform: uppercase; letter-spacing: .08em; margin-bottom: 8px;}
.hero-sample-id {font-size: 32px; font-weight: 800; margin-bottom: 10px;}
.kv-row {display:flex; gap: 8px; margin-bottom: 6px; font-size: 15px;}
.kv-key {width: 110px; color: var(--muted); flex: 0 0 auto;}

.note-box {
    border-left: 4px solid #f59e0b;
    background: #fffaf0;
    border-radius: 0 12px 12px 0;
    padding: 14px 16px;
    color: #7c5a10;
}

.vip-list-card {
    border: 1px solid var(--line);
    border-radius: 16px;
    background: var(--card);
    padding: 12px 14px;
    margin-bottom: 10px;
}
.vip-list-card.selected {
    border-color: #f3c77a;
    background: #fffaf1;
    box-shadow: 0 2px 8px rgba(183,121,31,0.08);
}
.vip-name-row {
    display: flex;
    align-items: center;
    gap: 10px;
    color: var(--ink);
    font-weight: 700;
    font-size: 15px;
}
.vip-icon {
    width: 28px;
    height: 28px;
    border-radius: 999px;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    background: var(--gold-bg);
    color: var(--gold);
    font-size: 14px;
    flex: 0 0 auto;
}
.vip-sub {margin-top: 6px; color: var(--muted); font-size: 13px;}
.vip-hg-mini {
    margin-top: 6px;
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    color: #7c5a10;
    font-size: 12px;
}

.vip-panel-title {font-size: 30px; font-weight: 800; line-height: 1.1; margin-bottom: 6px;}
.vip-panel-hg {
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    color: #7c5a10;
    font-size: 15px;
    margin-bottom: 14px;
}
.vip-bio {font-size: 15px; line-height: 1.75; color: #374151;}
.meta-grid {display:grid; grid-template-columns: 160px 1fr; gap: 8px 12px; margin-top: 14px; font-size: 14px;}
.meta-k {color: var(--muted);}
.meta-v {color: var(--ink);}

.branch-tree {
    border: 1px solid var(--line);
    border-radius: 18px;
    background: linear-gradient(180deg, #ffffff 0%, #faf8f4 100%);
    padding: 18px 18px 20px;
    margin-top: 16px;
}
.branch-title {font-size: 12px; color: var(--muted); text-transform: uppercase; letter-spacing: .08em; margin-bottom: 14px;}
.branch-wrap {display:flex; flex-direction:column; align-items:center;}
.branch-node {
    display:inline-flex;
    align-items:center;
    gap:8px;
    padding:9px 14px;
    border-radius:999px;
    border:1px solid var(--line);
    background:#fff;
    font-size:14px;
    box-shadow: 0 1px 2px rgba(0,0,0,0.04);
}
.branch-node.ancestor {color: var(--muted);}
.branch-node.query {background: var(--gold-bg); border-color: #f3c77a; color: #8a6218; font-weight: 700;}
.branch-node.vip {background: var(--blue-bg); border-color: #bfdbfe; color: #1d4ed8; font-weight: 700;}
.branch-node.same {background: var(--green-bg); border-color: #a7f3d0; color: var(--green); font-weight: 700;}
.branch-trunk {width:2px; height:22px; background:#d1d5db;}
.branch-split {
    width:min(84%, 460px);
    height:24px;
    border-top:2px solid #d1d5db;
    border-left:2px solid #d1d5db;
    border-right:2px solid #d1d5db;
    border-radius:16px 16px 0 0;
}
.branch-children {width:100%; display:grid; grid-template-columns:1fr 1fr; gap:20px;}
.branch-child {display:flex; flex-direction:column; align-items:center;}
.branch-child-line {width:2px; height:18px; background:#d1d5db;}
.branch-label {font-size:12px; color: var(--muted); margin-top: 8px; text-align:center;}

.show-more-wrap {margin: 4px 0 14px;}

button[kind="secondary"] {
    border-radius: 999px !important;
}
</style>
""",
    unsafe_allow_html=True,
)

# ---------- Helpers ----------
def safe_text(v) -> str:
    try:
        if pd.isna(v):
            return "—"
    except Exception:
        pass
    return str(v) if v is not None else "—"




def get_query_name_labels(result) -> dict:
    system_label = "Y-DNA" if result.system == "y" else "mtDNA"
    terminal_name = safe_text(result.target)
    tree_name = safe_text(getattr(result, "resolved_target_for_tree", None))

    if result.system != "y":
        return {
            "system_label": system_label,
            "terminal_name": terminal_name,
            "isogg_name": terminal_name,
            "display_name": terminal_name,
        }

    if tree_name in {"—", "", terminal_name}:
        tree_name = terminal_name

    return {
        "system_label": system_label,
        "terminal_name": terminal_name,
        "isogg_name": tree_name,
        "display_name": f"{terminal_name} · {tree_name}",
    }


def format_y_haplogroup_pair(terminal_name: object, isogg_name: object) -> str:
    terminal_name = safe_text(terminal_name)
    isogg_name = safe_text(isogg_name)
    if isogg_name in {"—", "", terminal_name}:
        return terminal_name
    return f"Terminal: {terminal_name} | ISOGG: {isogg_name}"


def add_haplogroup_display_columns(df: pd.DataFrame, system: str) -> pd.DataFrame:
    if df is None or df.empty:
        return df

    out = df.copy()

    if system == "y":
        if "y_haplogroup" in out.columns or "y_haplogroup_isogg" in out.columns:
            out["haplogroup_display"] = [
                format_y_haplogroup_pair(
                    row.get("y_haplogroup"),
                    row.get("y_haplogroup_isogg"),
                )
                for _, row in out.iterrows()
            ]
    else:
        if "mt_haplogroup" in out.columns:
            out["haplogroup_display"] = out["mt_haplogroup"].fillna("—").astype(str)

    return out



COUNTRY_NAME_FIXES = {
    "United States of America": "United States",
    "Russian Federation": "Russia",
    "Czech Republic": "Czechia",
    "Türkiye": "Turkey",
}


def get_core_sample_columns(df: pd.DataFrame, system: str) -> pd.DataFrame:
    if df is None or df.empty:
        return df
    df = add_haplogroup_display_columns(df, system)
    core_cols = [
        "sample_id",
        "haplogroup_display",
        "date_mean_bp",
        "political_entity",
        "locality",
        "publication",
    ]
    available = [c for c in core_cols if c in df.columns]
    if not available:
        return df
    return df[available].copy()


def render_origin_map(country_summary: pd.DataFrame):
    if country_summary is None or country_summary.empty:
        st.info("No geographic data available for map visualization.")
        return

    df = country_summary.copy()
    required_cols = {"political_entity", "sample_count", "oldest_bp"}
    missing = required_cols - set(df.columns)
    if missing:
        st.warning(f"Map could not be generated. Missing columns: {', '.join(sorted(missing))}")
        return

    df["political_entity"] = df["political_entity"].astype(str).replace(COUNTRY_NAME_FIXES)
    df["sample_count"] = pd.to_numeric(df["sample_count"], errors="coerce").fillna(0)
    df["oldest_bp"] = pd.to_numeric(df["oldest_bp"], errors="coerce")

    fig = px.choropleth(
        df,
        locations="political_entity",
        locationmode="country names",
        color="sample_count",
        hover_name="political_entity",
        hover_data={
            "sample_count": True,
            "oldest_bp": True,
            "political_entity": False,
        },
        color_continuous_scale="YlOrBr",
    )

    fig.update_layout(
        margin=dict(l=0, r=0, t=10, b=0),
        coloraxis_colorbar_title="Sample count",
        geo=dict(
            showframe=False,
            showcoastlines=True,
            projection_type="equirectangular",
            bgcolor="rgba(0,0,0,0)",
        ),
    )

    st.plotly_chart(fig, use_container_width=True)


@st.cache_data(show_spinner=False)
def cached_run(system: str, target: str):
    return run_haplogroup_analysis(
        system=system,
        target=target,
        aadr_path=AADR_PATH,
        vip_path=VIP_PATH,
        vip_sheet="Y" if system == "y" else "mtDNA",
        mt_tree_path=MT_TREE_PATH,
        y_tree_path=Y_TREE_PATH,
        early_n=EARLY_N,
    )


@st.cache_data(show_spinner=False)
def fetch_wikipedia_intro(vip_name: str) -> dict:
    headers = {"User-Agent": "HaplogroupDiscover/1.0 (Streamlit app)"}

    def candidate_names(name: str) -> list[str]:
        cands = [name.strip()]
        for cleaned in [
            re.sub(r"\s*\([^)]*\)", "", name).strip(),
            re.sub(r"\s*\[[^]]*\]", "", name).strip(),
            re.sub(r"\s+\d{3,4}[^A-Za-z]*\d{0,4}\s*$", "", name).strip(),
        ]:
            if cleaned and cleaned not in cands:
                cands.append(cleaned)
        return cands

    def parse_summary(title: str):
        url = f"https://en.wikipedia.org/api/rest_v1/page/summary/{quote(title, safe='')}"
        resp = requests.get(url, headers=headers, timeout=8)
        if resp.status_code != 200:
            return None
        data = resp.json()
        extract = (data.get("extract") or "").strip()
        if not extract or data.get("type") == "disambiguation":
            return None
        return {
            "title": data.get("title") or title,
            "extract": extract,
            "page_url": ((data.get("content_urls") or {}).get("desktop") or {}).get("page"),
            "thumbnail": ((data.get("thumbnail") or {}).get("source")),
        }

    for cand in candidate_names(vip_name):
        result = parse_summary(cand)
        if result:
            return result

    search_url = "https://en.wikipedia.org/w/api.php"
    params = {
        "action": "query",
        "list": "search",
        "srsearch": vip_name,
        "format": "json",
        "srlimit": 1,
    }
    resp = requests.get(search_url, params=params, headers=headers, timeout=8)
    if resp.status_code == 200:
        hits = ((resp.json() or {}).get("query") or {}).get("search") or []
        if hits:
            title = hits[0].get("title")
            result = parse_summary(title)
            if result:
                return result

    return {
        "title": vip_name,
        "extract": "No Wikipedia summary was found for this person.",
        "page_url": None,
        "thumbnail": None,
    }


def render_metric_card(label: str, value: str):
    st.markdown(
        f"""
        <div class="metric-card">
            <div class="metric-label">{escape(label)}</div>
            <div class="metric-value">{escape(value)}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def render_branch_tree(
    common_ancestor: str,
    query_hg: str,
    vip_hg: str,
    vip_name: str,
    query_terminal: str | None = None,
    query_isogg: str | None = None,
):
    common_ancestor = safe_text(common_ancestor)
    query_hg = safe_text(query_hg)
    vip_hg = safe_text(vip_hg)
    vip_name = safe_text(vip_name)
    query_terminal = safe_text(query_terminal or query_hg)
    query_isogg = safe_text(query_isogg or query_hg)
    query_label = (
        escape(query_terminal)
        if query_terminal == query_isogg
        else f"{escape(query_terminal)}<br><span style=\"font-size:12px;color:#6b7280;\">ISOGG: {escape(query_isogg)}</span>"
    )

    if query_hg == vip_hg:
        html = f"""
        <div class="branch-tree">
            <div class="branch-title">Haplogroup lineage</div>
            <div class="branch-wrap">
                <div class="branch-node same"><span>◎</span><span>{escape(query_terminal)}</span></div>
                <div class="branch-label">Terminal: {escape(query_terminal)}<br>ISOGG: {escape(query_isogg)}<br>Query haplogroup and {escape(vip_name)} are the same clade</div>
            </div>
        </div>
        """
    else:
        html = f"""
        <div class="branch-tree">
            <div class="branch-title">Haplogroup lineage</div>
            <div class="branch-wrap">
                <div class="branch-node ancestor"><span>⟡</span><span>{escape(common_ancestor)}</span></div>
                <div class="branch-trunk"></div>
                <div class="branch-split"></div>
                <div class="branch-children">
                    <div class="branch-child">
                        <div class="branch-child-line"></div>
                        <div class="branch-node query"><span>⌘</span><span>{escape(query_terminal)}</span></div>
                        <div class="branch-label">Terminal: {escape(query_terminal)}<br>ISOGG: {escape(query_isogg)}</div>
                    </div>
                    <div class="branch-child">
                        <div class="branch-child-line"></div>
                        <div class="branch-node vip"><span>✦</span><span>{escape(vip_hg)}</span></div>
                        <div class="branch-label">{escape(vip_name)}</div>
                    </div>
                </div>
            </div>
        </div>
        """
    st.markdown(html, unsafe_allow_html=True)


# ---------- Header ----------
st.markdown('<div class="hero-title"> Haplogroup Discover</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="hero-sub">Enter an mtDNA or Y-DNA haplogroup to explore the earliest ancient samples, inferred origin, and notable individuals.</div>',
    unsafe_allow_html=True,
)

# ---------- Input ----------
with st.sidebar:
    st.header("Query")
    system = st.selectbox("System", ["y", "mt"], format_func=lambda x: "Y-DNA" if x == "y" else "mtDNA")
    target = st.text_input("Haplogroup", value="R-M269")
    run_btn = st.button("Run analysis", use_container_width=True)

# keep last result visible even after reruns
should_run = run_btn or (
    st.session_state.analysis_result is not None
    and st.session_state.analysis_target == target.strip()
    and st.session_state.analysis_system == system
)

if not should_run:
    st.markdown(
        '<div class="note-box">Enter a haplogroup in the sidebar and click <b>Run analysis</b> to explore ancient samples, inferred origin, and notable individuals.</div>',
        unsafe_allow_html=True,
    )
    st.stop()

if run_btn or st.session_state.analysis_result is None or st.session_state.analysis_target != target.strip() or st.session_state.analysis_system != system:
    try:
        with st.spinner("Running analysis..."):
            result = cached_run(system=system, target=target.strip())
        st.session_state.analysis_result = result
        st.session_state.analysis_target = target.strip()
        st.session_state.analysis_system = system
        st.session_state.selected_vip_idx = 0
        st.session_state.vip_show_all = False
    except Exception as e:
        st.error(f"Analysis failed: {e}")
        st.stop()
else:
    result = st.session_state.analysis_result

matched_samples = result.matched_samples
oldest_sample = result.oldest_sample
early_samples = result.early_samples
country_summary = result.country_summary
vip_matches = result.vip_matches
system_label = "Y-DNA" if result.system == "y" else "mtDNA"
name_info = get_query_name_labels(result)

# ---------- Hero/result section ----------
if result.system == "y":
    st.subheader(f"{name_info['terminal_name']} · {system_label}")
    st.markdown(
        f"""
        <div class="plain-card" style="margin-bottom:16px;">
            <div class="eyebrow">Haplogroup naming</div>
            <div class="kv-row"><span class="kv-key">Terminal name</span><span>{escape(name_info['terminal_name'])}</span></div>
            <div class="kv-row"><span class="kv-key">ISOGG name</span><span>{escape(name_info['isogg_name'])}</span></div>
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.subheader(f"{result.target} · {system_label}")
hero_left, _ = st.columns([1.2, 1])

with hero_left:
    if not oldest_sample.empty:
        row = oldest_sample.iloc[0]
        st.markdown(
            f"""
            <div class="hero-card">
                <div class="eyebrow">Earliest observed ancient sample</div>
                <div class="hero-sample-id">{escape(safe_text(row.get('sample_id')))}</div>
                {
                    f'<div class="kv-row"><span class="kv-key">Terminal name</span><span>{escape(safe_text(row.get("y_haplogroup")))}</span></div>'
                    f'<div class="kv-row"><span class="kv-key">ISOGG name</span><span>{escape(safe_text(row.get("y_haplogroup_isogg")))}</span></div>'
                    if result.system == "y"
                    else f'<div class="kv-row"><span class="kv-key">Haplogroup</span><span>{escape(safe_text(row.get("mt_haplogroup")))}</span></div>'
                }
                <div class="kv-row"><span class="kv-key">Date</span><span>{escape(safe_text(row.get('date_mean_bp')))} BP</span></div>
                <div class="kv-row"><span class="kv-key">Country</span><span>{escape(safe_text(row.get('political_entity')))}</span></div>
                <div class="kv-row"><span class="kv-key">Locality</span><span>{escape(safe_text(row.get('locality')))}</span></div>
                <div class="kv-row"><span class="kv-key">Publication</span><span>{escape(safe_text(row.get('publication')))}</span></div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.warning("No matching ancient samples were found.")

st.markdown("### Summary")
c1, c2, c3, c4, c5 = st.columns(5)
with c1:
    render_metric_card("Direct upstream", safe_text(result.upstream_label))
with c2:
    oldest_bp = safe_text(oldest_sample.iloc[0]["date_mean_bp"]) if not oldest_sample.empty else "—"
    render_metric_card("Oldest BP", oldest_bp)
with c3:
    render_metric_card("Candidate origin", safe_text(result.candidate_origin_country))
with c4:
    render_metric_card("Matched samples", str(len(matched_samples)))
with c5:
    render_metric_card("VIP matches", str(len(vip_matches)))

# ---------- Tabs ----------
tab1, tab2, tab3, tab4 = st.tabs(["Overview", "Ancient Samples", "Origin", "Notable People"])

with tab1:
    left, right = st.columns([1, 1.4])
    with left:
        st.markdown("#### Query info")
        st.write("**Target:**", result.target)
        st.write("**System:**", system_label)
        if result.system == "y":
            st.write("**Terminal name:**", name_info["terminal_name"])
            st.write("**ISOGG name:**", name_info["isogg_name"])
        else:
            st.write("**Haplogroup name:**", name_info["terminal_name"])
        st.write("**Resolved target for tree:**", safe_text(result.resolved_target_for_tree))
        st.write("**Direct upstream:**", safe_text(result.upstream_label))
        st.write("**Included labels count:**", len(result.included_labels))
        if st.checkbox("Show included labels", value=False):
            st.code(", ".join(result.included_labels[:200]))
    with right:
        st.markdown("#### Analysis summary")
        st.markdown(
            """
            - **Matched samples**: all ancient samples matching the query or included clades  
            - **Oldest sample**: the earliest matching sample in the dataset  
            - **Early samples**: the earliest subset used for summary display  
            - **Country summary**: country-level aggregation based on early samples  
            - **Candidate origin country**: the current rule-based inferred country
            """
        )

with tab2:
    oldest_sample_display = get_core_sample_columns(oldest_sample, result.system)
    early_samples_display = get_core_sample_columns(early_samples, result.system)
    matched_samples_display = get_core_sample_columns(matched_samples, result.system)

    st.markdown("#### Oldest sample")
    st.dataframe(oldest_sample_display, use_container_width=True, hide_index=True)

    st.markdown("#### Early samples")
    st.dataframe(early_samples_display, use_container_width=True, hide_index=True)

    st.markdown("#### All matched samples")
    st.dataframe(matched_samples_display, use_container_width=True, hide_index=True)

    if not matched_samples.empty:
        st.download_button(
            "Download matched samples CSV",
            data=matched_samples.to_csv(index=False).encode("utf-8-sig"),
            file_name=f"{result.target}_matched_samples.csv",
            mime="text/csv",
        )

with tab3:
    st.markdown("#### Candidate origin")
    if result.candidate_origin_country:
        st.success(f"Candidate origin country: **{result.candidate_origin_country}**")
    else:
        st.warning("Insufficient data to infer a candidate origin country.")

    st.info(
        "Origin is inferred from the geographic distribution of the earliest matched ancient samples. "
        "The map below shows where those early samples are concentrated."
    )

    st.markdown("#### Geographic distribution of early samples")
    render_origin_map(country_summary)

    st.markdown("#### Country summary")
    st.dataframe(country_summary, use_container_width=True, hide_index=True)

    if not country_summary.empty:
        st.caption(
            "Darker shading indicates countries with more early matched samples. "
            "Hover over a country to view sample count and oldest BP."
        )

with tab4:
    st.markdown("#### Notable People")

    if vip_matches.empty:
        st.warning("No notable individuals matched this haplogroup.")
    else:
        keyword = st.text_input("Search by name", value="", key="vip_search")
        filtered = vip_matches.copy()

        if keyword.strip():
            filtered = filtered[
                filtered["vip_name"].astype(str).str.contains(keyword, case=False, na=False)
            ]
            st.session_state.vip_show_all = True

        display_df = filtered.reset_index(drop=True)
        total = len(display_df)

        if total == 0:
            st.info("No matching people were found.")
        else:
            if st.session_state.selected_vip_idx is None or st.session_state.selected_vip_idx not in display_df.index:
                st.session_state.selected_vip_idx = int(display_df.index[0])

            show_all = st.session_state.vip_show_all or total <= PREVIEW_VIP_COUNT
            shown_df = display_df if show_all else display_df.head(PREVIEW_VIP_COUNT)

            if st.session_state.selected_vip_idx not in shown_df.index:
                shown_df = (
                    pd.concat([shown_df, display_df.loc[[st.session_state.selected_vip_idx]]])
                    .drop_duplicates()
                    .sort_index()
                )

            shown_text = f"{total} result(s)"
            if not show_all and total > PREVIEW_VIP_COUNT:
                shown_text += f" — showing first {PREVIEW_VIP_COUNT}"
            st.caption(shown_text)

            if total > PREVIEW_VIP_COUNT:
                st.markdown('<div class="show-more-wrap">', unsafe_allow_html=True)
                if not show_all:
                    if st.button(f"Show remaining {total - PREVIEW_VIP_COUNT}", key="show_more_vips"):
                        st.session_state.vip_show_all = True
                        st.rerun()
                else:
                    if st.button("Show fewer", key="show_less_vips"):
                        st.session_state.vip_show_all = False
                        st.rerun()
                st.markdown("</div>", unsafe_allow_html=True)

            list_col, detail_col = st.columns([1, 1.45], gap="large")
            pending_select = None

            with list_col:
                for idx, row in shown_df.iterrows():
                    vip_name = safe_text(row.get("vip_name"))
                    vip_hg = safe_text(row.get("vip_haplogroup"))
                    common_ancestor = safe_text(row.get("common_ancestor"))
                    is_selected = st.session_state.selected_vip_idx == idx
                    card_class = "vip-list-card selected" if is_selected else "vip-list-card"
                    icon = "✦" if is_selected else "◌"

                    st.markdown(f'<div class="{card_class}">', unsafe_allow_html=True)
                    clicked = st.button(f"{icon}  {vip_name}", key=f"vip_pick_{idx}", use_container_width=True)
                    st.markdown(
                        f'<div class="vip-sub">Common ancestor: {escape(common_ancestor)}</div>'
                        f'<div class="vip-hg-mini">{escape(vip_hg)}</div>'
                        f"</div>",
                        unsafe_allow_html=True,
                    )
                    if clicked and not is_selected:
                        pending_select = idx

            if pending_select is not None:
                st.session_state.selected_vip_idx = pending_select
                st.rerun()

            with detail_col:
                sel_idx = st.session_state.selected_vip_idx
                if sel_idx is None or sel_idx not in display_df.index:
                    st.info("Select a person from the list to view details.")
                else:
                    row = display_df.loc[sel_idx]
                    vip_name = safe_text(row.get("vip_name"))
                    vip_hg = safe_text(row.get("vip_haplogroup"))
                    common_ancestor = safe_text(row.get("common_ancestor"))
                    relation = safe_text(row.get("relation"))
                    source = safe_text(row.get("source"))
                    note = safe_text(row.get("note"))

                    wiki = fetch_wikipedia_intro(vip_name)
                    wiki_title = safe_text(wiki.get("title"))
                    wiki_extract = safe_text(wiki.get("extract"))
                    wiki_url = wiki.get("page_url")

                    if result.system == "y":
                        st.markdown(
                            f'''
                            <div class="plain-card" style="margin-bottom:16px;">
                                <div class="eyebrow">Query haplogroup naming</div>
                                <div class="kv-row"><span class="kv-key">Terminal name</span><span>{escape(name_info["terminal_name"])}</span></div>
                                <div class="kv-row"><span class="kv-key">ISOGG name</span><span>{escape(name_info["isogg_name"])}</span></div>
                            </div>
                            ''',
                            unsafe_allow_html=True,
                        )

                    st.markdown('<div class="plain-card">', unsafe_allow_html=True)
                    if wiki_url:
                        st.markdown(
                            f'<div class="vip-panel-title"><a href="{escape(wiki_url)}" target="_blank" style="color:inherit;text-decoration:none;">{escape(wiki_title)}</a></div>',
                            unsafe_allow_html=True,
                        )
                    else:
                        st.markdown(f'<div class="vip-panel-title">{escape(wiki_title or vip_name)}</div>', unsafe_allow_html=True)

                    st.markdown(f'<div class="vip-panel-hg">{escape(vip_hg)}</div>', unsafe_allow_html=True)
                    st.markdown(f'<div class="vip-bio">{escape(wiki_extract)}</div>', unsafe_allow_html=True)

                    st.markdown(
                        f"""
                        <div class="meta-grid">
                            <div class="meta-k">Relation to query</div><div class="meta-v">{escape(relation)}</div>
                            <div class="meta-k">Common ancestor</div><div class="meta-v">{escape(common_ancestor)}</div>
                            <div class="meta-k">Source</div><div class="meta-v">{escape(source)}</div>
                            <div class="meta-k">Note</div><div class="meta-v">{escape(note)}</div>
                        </div>
                        """,
                        unsafe_allow_html=True,
                    )
                    st.markdown("</div>", unsafe_allow_html=True)

                    render_branch_tree(common_ancestor, result.target, vip_hg, vip_name, name_info["terminal_name"], name_info["isogg_name"])

            st.download_button(
                "Download VIP results CSV",
                data=filtered.to_csv(index=False).encode("utf-8-sig"),
                file_name=f"{result.target}_vip_matches.csv",
                mime="text/csv",
            )