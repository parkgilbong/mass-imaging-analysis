# Mass Imaging Analysis Pipeline - Snakemake Workflow
# This Snakefile replaces the main.ipynb workflow

import os

# Default config file (can be overridden via CLI with --configfile)
configfile: "config/config.yaml"

# Store the config file path for passing to scripts
# This can be overridden via: snakemake --configfile config/config_custom.yaml
CONFIG_FILE = config.get("__config_file__", "config/config.yaml")

# Extract configuration
DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
LOG_DIR = os.path.join(OUTPUT_DIR, "logs")
GROUPS = [g["name"] for g in config["group_info"]]
ROIS = [r["name"] for r in config["roi_info"]]

# Helper function to get expected files based on group configuration
def get_group_files(wildcards):
    """Generate list of expected imzML files for a group"""
    files = []
    for group_dict in config["group_info"]:
        if group_dict["name"] == wildcards.group:
            n_per_group = group_dict["n_per_group"]
            num_serial = group_dict["num_serial"]
            
            # Handle both integer and list formats
            if isinstance(num_serial, int):
                individuals_info = [(n, num_serial) for n in range(1, n_per_group + 1)]
            else:  # list
                individuals_info = [(n+1, serial_count) for n, serial_count in enumerate(num_serial)]
            
            for roi_dict in config["roi_info"]:
                roi = roi_dict["name"]
                for n, num_s in individuals_info:
                    for s in range(1, num_s + 1):
                        filename = f"{wildcards.group} {n}-{s} {roi}-total ion count.imzML"
                        files.append(os.path.join(DATA_DIR, filename))
    return files

# Get all expected mean intensity files for aggregation
def get_mean_intensity_files():
    """Generate list of all mean intensity CSV files"""
    files = []
    for group_dict in config["group_info"]:
        group = group_dict["name"]
        n_per_group = group_dict["n_per_group"]
        num_serial = group_dict["num_serial"]
        
        # Handle both integer and list formats
        if isinstance(num_serial, int):
            individuals_info = [(n, num_serial) for n in range(1, n_per_group + 1)]
        else:  # list
            individuals_info = [(n+1, serial_count) for n, serial_count in enumerate(num_serial)]
        
        for roi_dict in config["roi_info"]:
            roi = roi_dict["name"]
            for n, num_s in individuals_info:
                for s in range(1, num_s + 1):
                    base_name = f"{group} {n}-{s} {roi}-total ion count"
                    files.append(os.path.join(OUTPUT_DIR, f"{base_name}_mean_intensities.csv"))
    return files

# Rule: Final target - all analysis outputs
rule all:
    input:
        # Main statistical results
        os.path.join(OUTPUT_DIR, "statistical_results_main.csv"),
        # Post-hoc results (if generated)
        os.path.join(OUTPUT_DIR, "statistical_results_posthoc.csv"),
        # Plots for each ROI
        expand(os.path.join(OUTPUT_DIR, "plot_roi_{roi}.png"), roi=ROIS),
        # Log files
        expand(os.path.join(LOG_DIR, "{step}.log"), step=["parse_data", "aggregate_data", "analyze_data"])

# Rule: Step 1 - Parse imzML files and extract m/z bin intensities
rule parse_data:
    input:
        config_file = CONFIG_FILE,
        bins_file = lambda wildcards: config["binning_settings"]["file_path"] if config["binning_settings"]["mode"] == "file" else []
    output:
        # Output all mean intensity CSV files
        mean_files = get_mean_intensity_files()
    log:
        os.path.join(LOG_DIR, "parse_data.log")
    shell:
        """
        mkdir -p {LOG_DIR}
        python src/main.py --config {input.config_file} --log-file {log} 2>&1 | tee -a {log}
        """

# Rule: Step 2 - Aggregate mean intensities across technical replicates
rule aggregate_data:
    input:
        config_file = CONFIG_FILE,
        mean_files = get_mean_intensity_files()
    output:
        aggregated = os.path.join(OUTPUT_DIR, "aggregated_mean_intensities.csv")
    log:
        os.path.join(LOG_DIR, "aggregate_data.log")
    shell:
        """
        mkdir -p {LOG_DIR}
        python src/aggregate.py --config {input.config_file} --log-file {log} 2>&1 | tee -a {log}
        """

# Rule: Step 3 - Statistical analysis and visualization
rule analyze_data:
    input:
        config_file = CONFIG_FILE,
        aggregated = os.path.join(OUTPUT_DIR, "aggregated_mean_intensities.csv")
    output:
        stats_main = os.path.join(OUTPUT_DIR, "statistical_results_main.csv"),
        stats_posthoc = os.path.join(OUTPUT_DIR, "statistical_results_posthoc.csv"),
        plots = expand(os.path.join(OUTPUT_DIR, "plot_roi_{roi}.png"), roi=ROIS)
    log:
        os.path.join(LOG_DIR, "analyze_data.log")
    shell:
        """
        mkdir -p {LOG_DIR}
        python src/analysis.py --config {input.config_file} --log-file {log} 2>&1 | tee -a {log}
        """

# Rule: Clean output directory (optional)
rule clean:
    shell:
        """
        rm -rf {OUTPUT_DIR}/*
        echo "Cleaned output directory: {OUTPUT_DIR}"
        """

# Rule: Dry run to show what would be executed
rule dry_run:
    shell:
        """
        snakemake --dry-run --printshellcmds
        """
