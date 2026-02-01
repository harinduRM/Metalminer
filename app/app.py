import os
import pandas as pd
import streamlit as st
from ccdc import io

from metalminer.Manual import manual_disorder_html, parse_manual_input
from metalminer.core import Config, process_structure_core, run_pipeline
from metalminer.ligands import analyze_ligands
from metalminer.oxidation import calculate_bvs
from metalminer.utils import generate_enriched_xyz_string, make_xyz_viewer_html
from metalminer.validation import validate_geometry


def parse_list(text_value):
    if not text_value:
        return []
    items = []
    for part in text_value.split(","):
        part = part.strip()
        if part:
            items.append(part)
    return items


def parse_ids_from_upload(uploaded_file):
    if not uploaded_file:
        return []
    content = uploaded_file.getvalue().decode("utf-8")
    lines = [line.strip() for line in content.splitlines() if line.strip()]
    return lines

def trigger_rerun():
    if hasattr(st, "rerun"):
        st.rerun()
    st.experimental_rerun()

def init_manual_refinement_state(pkl_path, csv_path, config):
    pkl_mtime = os.path.getmtime(pkl_path)
    active_key = f"{pkl_path}:{csv_path}:{pkl_mtime}"
    if st.session_state.get("manual_active_key") == active_key:
        return
    df = pd.read_pickle(pkl_path)
    indices = df[(df["VALIDATION_FAILED"] == True) | (df["HAD_RDKIT_ISSUE"] == True)].index.tolist()
    st.session_state["manual_active_key"] = active_key
    st.session_state["manual_df"] = df
    st.session_state["manual_indices"] = indices
    st.session_state["manual_cursor"] = 0
    st.session_state["manual_stage"] = "edit"
    st.session_state["manual_preview"] = None
    st.session_state["manual_command_text"] = ""
    st.session_state["manual_pkl_path"] = pkl_path
    st.session_state["manual_csv_path"] = csv_path
    st.session_state["manual_num_metal_layers"] = int(config.num_metal_layers)
    st.session_state["manual_hydrogen_method"] = config.Hydrogen_Addition_method
    st.session_state["manual_extraction_cutoffs"] = config.Extraction_Cut_off_distances
    st.session_state["manual_extraction_method"] = config.EXTRACTION_METHOD
    st.session_state["manual_geometric_radius"] = float(config.Geometric_radius)

def save_manual_outputs():
    df = st.session_state["manual_df"]
    df.to_pickle(st.session_state["manual_pkl_path"])
    df.to_csv(st.session_state["manual_csv_path"], index=False)

def clear_manual_refinement_state():
    for key in [
        "manual_active_key",
        "manual_df",
        "manual_indices",
        "manual_cursor",
        "manual_stage",
        "manual_preview",
        "manual_command_text",
        "manual_pkl_path",
        "manual_csv_path",
        "manual_num_metal_layers",
        "manual_hydrogen_method",
        "manual_extraction_cutoffs",
        "manual_extraction_method",
        "manual_geometric_radius",
    ]:
        st.session_state.pop(key, None)

def advance_manual_cursor():
    indices = st.session_state["manual_indices"]
    cursor = st.session_state["manual_cursor"]
    if indices:
        indices.pop(cursor)
    if cursor >= len(indices):
        cursor = len(indices)
    st.session_state["manual_cursor"] = cursor
    st.session_state["manual_stage"] = "edit"
    st.session_state["manual_preview"] = None
    st.session_state["manual_command_text"] = ""

def render_manual_refinement_ui():
    st.subheader("Manual Refinement")
    indices = st.session_state.get("manual_indices", [])
    if not indices:
        st.success("No entries require manual refinement.")
        return
    cursor = st.session_state.get("manual_cursor", 0)
    if cursor >= len(indices):
        st.success("Manual refinement complete.")
        save_manual_outputs()
        return

    df = st.session_state["manual_df"]
    idx = indices[cursor]
    row = df.loc[idx]
    ref = row["CSD ID"]
    site = row["Site Label"]
    target_metal = row["Metal"]
    num_metal_layers = st.session_state["manual_num_metal_layers"]
    extraction_cutoffs = st.session_state["manual_extraction_cutoffs"]
    hydrogen_method = st.session_state["manual_hydrogen_method"]
    extraction_method = st.session_state.get("manual_extraction_method", "Topological")
    geometric_radius = st.session_state.get("manual_geometric_radius", 3.8)

    st.write(f"Entry: {ref} | Site: {site} | Metal: {target_metal}")
    st.write(
        f"Validation Failed: {row['VALIDATION_FAILED']} | RDKit Issue: {row['HAD_RDKIT_ISSUE']}"
    )

    entry_reader = io.EntryReader("CSD")
    entry = entry_reader.entry(ref)
    if not entry:
        st.warning(f"Could not load entry {ref} from the CSD.")
        if st.button("Skip entry", key=f"manual_skip_{idx}"):
            advance_manual_cursor()
            trigger_rerun()
        return

    stage = st.session_state.get("manual_stage", "edit")
    if stage == "edit":
        raw_mol, _, _, _, _, _ = process_structure_core(
            entry,
            target_metal,
            site,
            num_metal_layers,
            perform_correction=False,
            EXTRACTION_METHOD=extraction_method,
            Extraction_Cut_off_distances=extraction_cutoffs,
            Geometric_radius=geometric_radius,
        )
        if raw_mol is not None:
            raw_html = manual_disorder_html(raw_mol, site)
            st.components.v1.html(raw_html, height=520, scrolling=False)
        else:
            st.warning("Unable to render the raw structure for this entry.")

        st.session_state["manual_command_text"] = st.text_area(
            "Manual refinement command",
            value=st.session_state.get("manual_command_text", ""),
            help="Example: Atom_delete = [1,2], Atom_average = [[3,4]]",
        )
        col_apply, col_skip, col_delete = st.columns(3)
        if col_apply.button("Apply command", key=f"manual_apply_{idx}"):
            command_text = st.session_state.get("manual_command_text", "")
            action, commands = parse_manual_input(command_text)
            if action == "Delete":
                df.drop(idx, inplace=True)
                save_manual_outputs()
                advance_manual_cursor()
                trigger_rerun()
            mol, metal, _, _, bridge, val_f = process_structure_core(
                entry,
                target_metal,
                site,
                num_metal_layers,
                manual_commands=commands,
                perform_correction=True,
                EXTRACTION_METHOD=extraction_method,
                Extraction_Cut_off_distances=extraction_cutoffs,
                Geometric_radius=geometric_radius,
            )
            if mol is None:
                st.error("Manual refinement failed to generate a structure for preview.")
                return
            lig_data, has_b, _, v_hs = analyze_ligands(
                mol, metal, entry.chemical_name, bridge, hydrogen_method
            )
            is_actually_valid = validate_geometry(mol, metal, virtual_hydrogens=v_hs, check_rdkit=True)
            val_f = not is_actually_valid
            xyz_preview = generate_enriched_xyz_string(mol, v_hs)
            st.session_state["manual_preview"] = {
                "lig_data": lig_data,
                "has_b": has_b,
                "val_f": val_f,
                "xyz_preview": xyz_preview,
                "metal": metal,
            }
            st.session_state["manual_stage"] = "preview"
            trigger_rerun()

        if col_skip.button("Skip entry", key=f"manual_skip_{idx}"):
            advance_manual_cursor()
            trigger_rerun()
        if col_delete.button("Delete entry", key=f"manual_delete_{idx}"):
            df.drop(idx, inplace=True)
            save_manual_outputs()
            advance_manual_cursor()
            trigger_rerun()

    if stage == "preview":
        preview = st.session_state.get("manual_preview")
        if not preview:
            st.warning("No preview data available.")
            st.session_state["manual_stage"] = "edit"
            trigger_rerun()
        preview_html = make_xyz_viewer_html(preview["xyz_preview"], width=450, height=450)
        st.components.v1.html(preview_html, height=480, scrolling=False)
        st.write(f"Validation Result: {'FAILED' if preview['val_f'] else 'PASSED'}")
        col_accept, col_back = st.columns(2)
        if col_accept.button("Accept changes", key=f"manual_accept_{idx}"):
            lig_data = preview["lig_data"]
            df.at[idx, "Coordination Number"] = sum(len(x["Connecting_atom"]) for x in lig_data)
            df.at[idx, "Ligand_count"] = sum(x["number"] for x in lig_data)
            df.at[idx, "AN_n_Ligand_Info"] = str(lig_data)
            df.at[idx, "VALIDATION_FAILED"] = preview["val_f"]
            df.at[idx, "HAD_RDKIT_ISSUE"] = any(
                item.get("HAD_RDKIT_ISSUE", False) for item in lig_data
            )
            df.at[idx, "Bridging_Ligand"] = preview["has_b"]
            df.at[idx, "Total_Ligand_Charge"] = sum(
                item["Charge"] * item["number"] for item in lig_data
            )
            df.at[idx, "AN_n_XYZ Coordinates"] = preview["xyz_preview"]
            if row["OS_Method"] == "BVS (Geometry)":
                updated_bvs = calculate_bvs(preview["metal"], preview["metal"].neighbours)
                df.at[idx, "Oxidation State"] = f"+{round(updated_bvs)}"
            save_manual_outputs()
            advance_manual_cursor()
            trigger_rerun()
        if col_back.button("Back to edit", key=f"manual_back_{idx}"):
            st.session_state["manual_stage"] = "edit"
            st.session_state["manual_preview"] = None
            trigger_rerun()

st.set_page_config(page_title="MetalMiner", layout="wide")
st.title("MetalMiner Local GUI")
st.caption("Streamlit wrapper for the MetalMiner pipeline.")

default_config = Config()

with st.form("metalminer_form"):
    col_left, col_right = st.columns(2)

    with col_left:
        target_metals = st.text_input(
            "Target metals (comma-separated)",
            value=",".join(default_config.TARGET_METAL_list),
            help="Elements to search for. Example: Pu, U, Th",
        )
        target_csd_ids_text = st.text_area(
            "Target CSD IDs (optional, comma or newline-separated)",
            value="",
            help="Example: ABABEF, CSDABC (one per line or comma-separated).",
        )
        target_csd_ids_file = st.file_uploader(
            "Or upload a text/CSV file with CSD IDs",
            type=["txt", "csv"],
            help="One CSD ID per line in a .txt or .csv file.",
        )
        extraction_method = st.selectbox(
            "Extraction method",
            options=["Topological", "Geometric"],
            index=0 if default_config.EXTRACTION_METHOD == "Topological" else 1,
            help="Topological uses bond cutoffs; Geometric uses a radial search.",
        )
        hydrogen_method = st.selectbox(
            "Hydrogen addition method",
            options=["None", "Geometric", "RDkit"],
            index=["None", "Geometric", "RDkit"].index(default_config.Hydrogen_Addition_method)
            if default_config.Hydrogen_Addition_method in ["None", "Geometric", "RDkit"]
            else 1,
            help="Example: Geometric adds hydrogens by geometry; RDkit uses RDKit.",
        )
        disorder_resolve_method = st.selectbox(
            "Disorder resolve method",
            options=["Major_componet", "Minor_componet", "Average", "Hybrid"],
            index=3 if default_config.DISORDER_RESOLVE_METHOD == "Hybrid" else 0,
            help="How to resolve disorder. Example: Hybrid for mixed strategies.",
        )
        oxidation_filter = st.text_input(
            "Oxidation states filter (comma-separated or 'all')",
            value=",".join(default_config.OXIDATION_STATES_FILTER),
            help="Example: all or +2,+3",
        )
        metalloligands = st.text_input(
            "Metalloligands (comma-separated)",
            value=",".join(default_config.metalloligands),
            help="Example: Cr, V, Mo, W",
        )
        r_factor_limit = st.number_input(
            "R-factor limit",
            value=float(default_config.R_FACTOR_LIMIT),
            help="Example: 10.0 (structures above this are skipped).",
        )

    with col_right:
        process_limit = st.number_input(
            "Process limit",
            min_value=0,
            value=int(default_config.PROCESS_LIMIT),
            help="Example: 1000 (0 means no limit).",
        )
        site_timeout_seconds = st.number_input(
            "Per-site timeout (seconds)",
            min_value=0,
            value=int(default_config.SITE_TIMEOUT_SECONDS or 0),
            help="0 disables timeouts. Example: 600 to skip sites taking longer than 10 minutes.",
        )
        retry_timeouts = st.checkbox(
            "Retry timed-out sites on resume",
            value=default_config.RETRY_TIMEOUTS,
            help="If unchecked, timed-out sites are treated as completed when resuming.",
        )
        visualize = st.checkbox(
            "Visualize (py3Dmol)",
            value=default_config.VISUALIZE,
            help="Show 3D views during processing.",
        )
        visualization_limit = st.number_input(
            "Visualization limit",
            min_value=1,
            value=int(default_config.VISUALIZATION_LIMIT),
            help="Example: 15 (max number of structures to render).",
        )
        num_metal_layers = st.number_input(
            "Number of metal layers",
            min_value=1,
            value=int(default_config.num_metal_layers),
            help="Example: 1 or 2 (how many metal shells to include).",
        )
        limit_nm_nm = st.number_input(
            "Extraction cutoff: LIMIT_NM_NM",
            value=float(default_config.Extraction_Cut_off_distances["LIMIT_NM_NM"]),
            help="Only used when the extraction method is Topological.",
        )
        limit_m_nm = st.number_input(
            "Extraction cutoff: LIMIT_M_NM",
            value=float(default_config.Extraction_Cut_off_distances["LIMIT_M_NM"]),
            help="Only used when the extraction method is Topological.",
        )
        limit_h_x = st.number_input(
            "Extraction cutoff: LIMIT_H_X",
            value=float(default_config.Extraction_Cut_off_distances["LIMIT_H_X"]),
            help="Only used when the extraction method is Topological.",
        )
        geometric_radius = st.number_input(
            "Geometric radius",
            min_value=0.1,
            value=float(default_config.Geometric_radius),
            help="Only used when the extraction method is Geometric.",
        )

        filter_polymeric = st.checkbox(
            "Filter polymeric",
            value=default_config.FILTER_POLYMERIC,
            help="Exclude polymeric structures from results.",
        )
        filter_powder = st.checkbox(
            "Filter powder",
            value=default_config.FILTER_POWDER,
            help="Exclude powder diffraction structures.",
        )
        filter_all_disorder = st.checkbox(
            "Filter all disorder",
            value=default_config.FILTER_ALL_DISORDER,
            help="Exclude structures with global disorder.",
        )
        filter_primary_disorder = st.checkbox(
            "Filter primary disorder",
            value=default_config.FILTER_PRIMARY_DISORDER,
            help="Exclude structures with primary sphere disorder.",
        )
        correct_primary_disorder = st.checkbox(
            "Correct primary disorder",
            value=default_config.CORRECT_PRIMARY_DISORDER,
            help="Attempt correction for primary sphere disorder.",
        )
        get_abstract = st.checkbox(
            "Fetch abstract (OpenAlex/Gemini)",
            value=default_config.GET_ABSTRACT,
            help="Fetch abstracts via OpenAlex and Gemini (if configured).",
        )
        edit_manual = st.checkbox(
            "Run manual refinement step",
            value=False,
            help="Manual refinement opens an interactive review panel in the app.",
        )

    submit = st.form_submit_button("Run")

console_placeholder = st.empty()
viewer_placeholder = st.empty()
visuals_placeholder = st.empty()

if submit:
    st.session_state["log_lines"] = []
    st.session_state["visualization_html"] = []

    def emit(text):
        st.session_state["log_lines"].append(text)
        console_placeholder.text_area(
            "Console Output",
            value="\n".join(st.session_state["log_lines"]),
            height=300,
        )

    def render_visualizations():
        html_list = st.session_state.get("visualization_html", [])
        if not html_list:
            return
        with visuals_placeholder.container():
            st.subheader("Visualizations")
            for idx, html in enumerate(html_list, start=1):
                st.markdown(f"View {idx}")
                st.components.v1.html(html, height=450, scrolling=False)

    def on_visualization(html, idx):
        st.session_state["visualization_html"].append(html)
        render_visualizations()

    target_csd_ids = parse_list(target_csd_ids_text)
    target_csd_ids.extend(parse_ids_from_upload(target_csd_ids_file))
    target_csd_ids = target_csd_ids or None

    manual_refinement_requested = edit_manual
    if not manual_refinement_requested:
        clear_manual_refinement_state()
    config = Config(
        TARGET_METAL_list=parse_list(target_metals),
        target_csd_ids=target_csd_ids,
        EXTRACTION_METHOD=extraction_method,
        R_FACTOR_LIMIT=r_factor_limit,
        PROCESS_LIMIT=int(process_limit),
        SITE_TIMEOUT_SECONDS=int(site_timeout_seconds) if int(site_timeout_seconds) > 0 else None,
        RETRY_TIMEOUTS=retry_timeouts,
        VISUALIZE=visualize,
        VISUALIZATION_LIMIT=int(visualization_limit),
        FILTER_POLYMERIC=filter_polymeric,
        FILTER_POWDER=filter_powder,
        FILTER_ALL_DISORDER=filter_all_disorder,
        FILTER_PRIMARY_DISORDER=filter_primary_disorder,
        CORRECT_PRIMARY_DISORDER=correct_primary_disorder,
        Hydrogen_Addition_method=hydrogen_method,
        DISORDER_RESOLVE_METHOD=disorder_resolve_method,
        num_metal_layers=int(num_metal_layers),
        OXIDATION_STATES_FILTER=parse_list(oxidation_filter) or ["all"],
        GET_ABSTRACT=get_abstract,
        Geometric_radius=float(geometric_radius),
        Edit_manual=False,
        Extraction_Cut_off_distances={
            "LIMIT_NM_NM": float(limit_nm_nm),
            "LIMIT_M_NM": float(limit_m_nm),
            "LIMIT_H_X": float(limit_h_x),
        },
        metalloligands=parse_list(metalloligands),
    )

    result = run_pipeline(
        config,
        emit=emit,
        display_views=False,
        on_visualization=on_visualization if visualize else None,
    )

    if visualize:
        if st.session_state["visualization_html"]:
            render_visualizations()
        else:
            viewer_placeholder.info("No py3Dmol visualizations produced for this run.")

    if manual_refinement_requested:
        pkl_paths = result.artifacts.get("pkl_paths", [])
        csv_paths = result.artifacts.get("csv_paths", [])
        if pkl_paths and csv_paths:
            init_manual_refinement_state(pkl_paths[-1], csv_paths[-1], config)
        else:
            clear_manual_refinement_state()
            st.warning("Manual refinement requested, but no result files were generated.")

if st.session_state.get("manual_active_key"):
    render_manual_refinement_ui()
