from metalminer.main import main as run_metalminer

if __name__ == "__main__":
    run_metalminer(
        TARGET_METAL_list=["Pa","Am","Cm","Bk","Cf"],
        target_csd_ids=None,
        OXIDATION_STATES_FILTER=["all"],
        EXTRACTION_METHOD="Topological",
        R_FACTOR_LIMIT=5,
        PROCESS_LIMIT=100000,
        VISUALIZE=False,
        FILTER_POLYMERIC=False,
        FILTER_POWDER=True,
        FILTER_ALL_DISORDER=False,
        FILTER_PRIMARY_DISORDER=False,
        CORRECT_PRIMARY_DISORDER=True,
        Hydrogen_Addition_method="Geometric",
        DISORDER_RESOLVE_METHOD="Hybrid",
        num_metal_layers=1,
        GET_ABSTRACT=True,
        Extraction_Cut_off_distances={"LIMIT_NM_NM": 2.8, "LIMIT_M_NM": 3.5, "LIMIT_H_X": 1.3},
        metalloligands=["Cr", "V", "Mo", "W", "Co"],
        Edit_manual=False,
        SITE_TIMEOUT_SECONDS=240,  # optional
    )
