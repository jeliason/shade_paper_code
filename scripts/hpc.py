#!/usr/bin/env python3
"""
HPC Workflow Tool for SHADE Paper Code

This script automates syncing code to HPC, submitting SLURM jobs,
checking status, and fetching results.

Usage:
    python scripts/hpc.py sync
    python scripts/hpc.py submit <sim_dir> <step>
    python scripts/hpc.py status [sim_dir]
    python scripts/hpc.py summary [sim_dir]
    python scripts/hpc.py fetch <sim_dir>
    python scripts/hpc.py logs <sim_dir> <step> [--lines N]

Requirements:
    pip install python-dotenv jinja2 pyyaml
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

try:
    from dotenv import load_dotenv
    import yaml
    from jinja2 import Template
except ImportError:
    print("Error: Missing required packages. Please install with:")
    print("  pip install python-dotenv jinja2 pyyaml")
    sys.exit(1)

# Load environment variables from .env
PROJECT_ROOT = Path(__file__).parent.parent
ENV_FILE = PROJECT_ROOT / ".env"

if not ENV_FILE.exists():
    print(f"Error: .env file not found at {ENV_FILE}")
    print("Please copy .env.example to .env and fill in your HPC details")
    sys.exit(1)

load_dotenv(ENV_FILE, override=True)

# Read environment variables
HPC_HOST = os.getenv("HPC_HOST")
HPC_USER = os.getenv("HPC_USER")
HPC_DIR = os.getenv("HPC_DIR")
HPC_EMAIL = os.getenv("HPC_EMAIL")
HPC_COMPILER_TOOLCHAIN = os.getenv("HPC_COMPILER_TOOLCHAIN")
HPC_R_MODULE = os.getenv("HPC_R_MODULE")
HPC_DATA_PATH = os.getenv("HPC_DATA_PATH", "./data")
SHADE_SOURCE = os.getenv("SHADE_SOURCE", "")
RSYNC_EXCLUDE = os.getenv("RSYNC_EXCLUDE", ".git").split(",")

# Validate required env vars
REQUIRED_VARS = ["HPC_HOST", "HPC_USER", "HPC_DIR", "HPC_EMAIL","HPC_COMPILER_TOOLCHAIN","HPC_R_MODULE","HPC_DATA_PATH"]
missing = [var for var in REQUIRED_VARS if not os.getenv(var)]
if missing:
    print(f"Error: Missing required environment variables in .env: {', '.join(missing)}")
    sys.exit(1)


def run_command(cmd, check=True, capture=False):
    """Run a shell command"""
    print(f"Running: {cmd}")
    if capture:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if check and result.returncode != 0:
            print(f"Error: Command failed with code {result.returncode}")
            print(result.stderr)
            sys.exit(1)
        return result.stdout
    else:
        result = subprocess.run(cmd, shell=True)
        if check and result.returncode != 0:
            print(f"Error: Command failed with code {result.returncode}")
            sys.exit(1)


def sync_to_hpc():
    """Rsync code to HPC"""
    print(f"\n=== Syncing code to {HPC_HOST}:{HPC_DIR} ===\n")

    # Build exclude arguments
    exclude_args = " ".join([f"--exclude='{pattern.strip()}'" for pattern in RSYNC_EXCLUDE])

    # Rsync command
    rsync_cmd = (
        f"rsync -avz --delete "
        f"{exclude_args} "
        f"{PROJECT_ROOT}/ "
        f"{HPC_HOST}:{HPC_DIR}/"
    )

    run_command(rsync_cmd)
    print("\n✓ Sync complete\n")


def load_config(sim_dir):
    """Load SLURM config for a simulation directory"""
    config_file = PROJECT_ROOT / sim_dir / "slurm_config.yaml"

    if not config_file.exists():
        print(f"Error: Config file not found: {config_file}")
        print(f"Please create {sim_dir}/slurm_config.yaml")
        sys.exit(1)

    with open(config_file) as f:
        return yaml.safe_load(f)




def generate_slurm_script(sim_dir, step, config, force=False):
    """Generate a SLURM script from template"""

    if step not in config["jobs"]:
        print(f"Error: Step '{step}' not found in {sim_dir}/slurm_config.yaml")
        print(f"Available steps: {', '.join(config['jobs'].keys())}")
        sys.exit(1)

    job_config = config["jobs"][step]

    # Load template
    template_file = PROJECT_ROOT / "scripts" / "slurm_template.sh"
    with open(template_file) as f:
        template = Template(f.read())

    # Prepare template variables
    log_prefix = f"{sim_dir.replace('sim_', '')}_{step.replace('0', '').replace('_', '')}"

    # Check array_range and show job configuration
    array_range = job_config.get("array_range")
    if array_range:
        # Show current job configuration
        print(f"  Current job configuration in slurm_config.yaml:")
        print(f"    - job_name: {job_config.get('job_name', f'{sim_dir}_{step}')}")
        print(f"    - array_range: {array_range}")
        print(f"    - time: {job_config.get('time', '1:00:00')}")
        print(f"    - mem: {job_config.get('mem', '4G')}")
        print(f"    - cpus: {job_config.get('cpus', 1)}")
        if force:
            print(f"    - force: True (will rerun completed simulations)")
        response = input("\n  Has the job configuration been updated correctly in slurm_config.yaml? (yes/no): ").strip().lower()
        if response not in ['yes', 'y']:
            print("\n  ✗ Job submission cancelled. Please update slurm_config.yaml and try again.")
            sys.exit(0)

    variables = {
        "job_name": job_config.get("job_name", f"{sim_dir}_{step}"),
        "log_prefix": log_prefix,
        "array_range": array_range,
        "time": job_config.get("time", "1:00:00"),
        "mem": job_config.get("mem", "4G"),
        "cpus": job_config.get("cpus", 1),
        "email": HPC_EMAIL,
        "compiler_toolchain": HPC_COMPILER_TOOLCHAIN,
        "r_module": HPC_R_MODULE,
        "data_path": HPC_DATA_PATH,
        "script_path": f"{sim_dir}/{step}.R",
        "force_rerun": force
    }

    # Render template
    slurm_script = template.render(**variables)

    # Write to file
    output_file = PROJECT_ROOT / sim_dir / f"{step}.slurm"
    with open(output_file, "w") as f:
        f.write(slurm_script)

    print(f"✓ Generated {output_file}")
    return output_file


def submit_job(sim_dir, step, force=False):
    """Submit a SLURM job to HPC"""
    print(f"\n=== Submitting {sim_dir}/{step} ===\n")

    if force:
        print("  FORCE MODE: Will rerun completed simulations\n")

    # Load config
    config = load_config(sim_dir)

    # Generate SLURM script locally
    generate_slurm_script(sim_dir, step, config, force=force)

    # Sync to HPC (including the generated .slurm file)
    sync_to_hpc()

    # Clear and create logs directory at project root on HPC
    print("  Clearing logs directory on HPC...")
    ssh_cmd = f"ssh {HPC_HOST} 'cd {HPC_DIR} && rm -rf logs && mkdir -p logs'"
    run_command(ssh_cmd)

    # Create data directory on HPC
    print(f"  Creating data directory: {HPC_DATA_PATH}")
    ssh_cmd = f"ssh {HPC_HOST} 'mkdir -p {HPC_DATA_PATH}'"
    run_command(ssh_cmd)

    # Submit job
    submit_cmd = (
        f"ssh {HPC_HOST} 'cd {HPC_DIR} && sbatch {sim_dir}/{step}.slurm'"
    )

    output = run_command(submit_cmd, capture=True)
    print(output)
    print("✓ Job submitted\n")


def check_status(sim_dir=None):
    """Check SLURM job status"""
    print(f"\n=== Checking job status ===\n")

    squeue_cmd = f"ssh {HPC_HOST} 'squeue -u {HPC_USER}'"

    if sim_dir:
        # Filter by job name pattern
        squeue_cmd += f" | grep {sim_dir}"

    run_command(squeue_cmd, check=False)


def job_summary(sim_dir=None):
    """Show job summary statistics (pending, running, completed, failed)"""
    print(f"\n=== Job Summary ===\n")

    # Build sacct command to get all recent jobs
    # Format: JobID, JobName, State, ExitCode
    sacct_cmd = f"ssh {HPC_HOST} 'sacct -u {HPC_USER} --starttime=today --format=JobID,JobName%30,State,ExitCode --noheader'"

    if sim_dir:
        print(f"Filtering by simulation: {sim_dir}\n")
        sacct_cmd += f" | grep {sim_dir}"

    # Get output
    output = run_command(sacct_cmd, check=False, capture=True)

    if not output or not output.strip():
        print("No jobs found.")
        return

    # Parse job states
    lines = output.strip().split('\n')

    # Count states (ignore .batch and .extern sub-jobs)
    states = {}
    for line in lines:
        parts = line.split()
        if len(parts) >= 3:
            job_id = parts[0]
            state = parts[2]

            # Skip sub-jobs (.batch, .extern, etc.)
            if '.' in job_id:
                continue

            # Normalize state names
            if state in ['PENDING', 'PD']:
                states['pending'] = states.get('pending', 0) + 1
            elif state in ['RUNNING', 'R']:
                states['running'] = states.get('running', 0) + 1
            elif state in ['COMPLETED', 'CD']:
                states['completed'] = states.get('completed', 0) + 1
            elif state in ['FAILED', 'F', 'TIMEOUT', 'TO', 'CANCELLED', 'CA', 'NODE_FAIL', 'NF']:
                states['failed'] = states.get('failed', 0) + 1
            else:
                states['other'] = states.get('other', 0) + 1

    # Print summary
    total = sum(states.values())
    print(f"Total jobs: {total}")
    print(f"  Pending:   {states.get('pending', 0):4d}")
    print(f"  Running:   {states.get('running', 0):4d}")
    print(f"  Completed: {states.get('completed', 0):4d}")
    print(f"  Failed:    {states.get('failed', 0):4d}")
    if states.get('other', 0) > 0:
        print(f"  Other:     {states.get('other', 0):4d}")

    # Show percentage completed
    if total > 0:
        pct_complete = (states.get('completed', 0) / total) * 100
        print(f"\nProgress: {pct_complete:.1f}% complete")

    print()


def job_batch_status(job_id):
    """Show detailed status for a specific job batch (array job)"""
    print(f"\n=== Job Batch Status: {job_id} ===\n")

    # Get job array details
    sacct_cmd = f"ssh {HPC_HOST} 'sacct -j {job_id} --format=JobID,JobName%30,State,ExitCode --noheader'"

    output = run_command(sacct_cmd, check=False, capture=True)

    if not output or not output.strip():
        print(f"No job found with ID {job_id}")
        return

    lines = output.strip().split('\n')

    # Count states (ignore .batch and .extern sub-jobs)
    states = {}
    failed_jobs = []

    for line in lines:
        parts = line.split()
        if len(parts) >= 3:
            job_id_str = parts[0]
            state = parts[2]
            exit_code = parts[3] if len(parts) >= 4 else "0:0"

            # Skip sub-jobs (.batch, .extern, etc.)
            if '.batch' in job_id_str or '.extern' in job_id_str:
                continue

            # Track failed jobs
            if state in ['FAILED', 'F', 'TIMEOUT', 'TO', 'CANCELLED', 'CA', 'NODE_FAIL', 'NF']:
                failed_jobs.append((job_id_str, state, exit_code))

            # Normalize state names
            if state in ['PENDING', 'PD']:
                states['pending'] = states.get('pending', 0) + 1
            elif state in ['RUNNING', 'R']:
                states['running'] = states.get('running', 0) + 1
            elif state in ['COMPLETED', 'CD']:
                states['completed'] = states.get('completed', 0) + 1
            elif state in ['FAILED', 'F', 'TIMEOUT', 'TO', 'CANCELLED', 'CA', 'NODE_FAIL', 'NF']:
                states['failed'] = states.get('failed', 0) + 1
            else:
                states['other'] = states.get('other', 0) + 1

    # Print summary
    total = sum(states.values())
    print(f"Total array tasks: {total}")
    print(f"  Pending:   {states.get('pending', 0):4d}")
    print(f"  Running:   {states.get('running', 0):4d}")
    print(f"  Completed: {states.get('completed', 0):4d}")
    print(f"  Failed:    {states.get('failed', 0):4d}")
    if states.get('other', 0) > 0:
        print(f"  Other:     {states.get('other', 0):4d}")

    # Show percentage completed
    if total > 0:
        pct_complete = (states.get('completed', 0) / total) * 100
        print(f"\nProgress: {pct_complete:.1f}% complete")

    # Show failed jobs if any
    if failed_jobs:
        print(f"\n=== Failed Jobs ({len(failed_jobs)}) ===")
        print(f"{'Job ID':<20} {'State':<15} {'Exit Code'}")
        print("-" * 50)
        for job_id_str, state, exit_code in failed_jobs[:20]:  # Show first 20
            print(f"{job_id_str:<20} {state:<15} {exit_code}")
        if len(failed_jobs) > 20:
            print(f"\n... and {len(failed_jobs) - 20} more failed jobs")

    print()


def fetch_results(sim_dir):
    """Fetch results from HPC

    Args:
        sim_dir: Simulation directory path (e.g., 'sim_timing' or 'sim_shade_comparison/compartments')
    """
    print(f"\n=== Fetching results for {sim_dir} ===\n")

    # Fetch only the analysis summary file (not raw simulation data)
    print(f"  Fetching analysis_summary.rds from {HPC_DATA_PATH}/{sim_dir}/")

    # Parse sim_dir to determine local path structure
    # For 'sim_shade_comparison/compartments': store in sim_shade_comparison/data/compartments/
    # For 'sim_timing': store in sim_timing/data/
    parts = sim_dir.split('/')
    if len(parts) > 1:
        # Nested path: sim_shade_comparison/compartments -> sim_shade_comparison/data/compartments
        base_sim = parts[0]
        subdir = '/'.join(parts[1:])
        local_data_dir = PROJECT_ROOT / base_sim / "data" / subdir
    else:
        # Simple path: sim_timing -> sim_timing/data
        local_data_dir = PROJECT_ROOT / sim_dir / "data"

    local_data_dir.mkdir(parents=True, exist_ok=True)

    rsync_cmd = (
        f"rsync -avz "
        f"{HPC_HOST}:{HPC_DATA_PATH}/{sim_dir}/analysis_summary.rds "
        f"{local_data_dir}/"
    )

    run_command(rsync_cmd)

    # Also fetch figures if they exist (from project directory, not data directory)
    # Extract base simulation name (first part of path) for figures
    base_sim = sim_dir.split('/')[0]
    print(f"  Fetching figures from {HPC_DIR}/{base_sim}/")

    local_fig_dir = PROJECT_ROOT / base_sim / "figures"
    rsync_cmd_figs = (
        f"rsync -avz "
        f"{HPC_HOST}:{HPC_DIR}/{base_sim}/figures/ "
        f"{local_fig_dir}/ 2>/dev/null || true"
    )

    run_command(rsync_cmd_figs, check=False)

    print("\n✓ Results fetched\n")


def view_logs(sim_dir, step, lines=20):
    """View recent log files from HPC"""
    print(f"\n=== Viewing logs for {sim_dir}/{step} (last {lines} lines) ===\n")

    # Find most recent log files in project-level logs directory
    log_pattern = f"logs/*{step.replace('0', '').replace('_', '')}*.out"

    ssh_cmd = (
        f"ssh {HPC_HOST} '"
        f"cd {HPC_DIR} && "
        f"latest_log=$(ls -t {log_pattern} 2>/dev/null | head -1) && "
        f"if [ -n \"$latest_log\" ]; then "
        f"echo \"File: $latest_log\" && echo && tail -n {lines} \"$latest_log\"; "
        f"else echo \"No log files found for {step}\"; fi'"
    )

    run_command(ssh_cmd, check=False)


def install_packages(force_reinstall=False):
    """Install R packages on HPC"""
    print(f"\n=== Installing R packages on HPC ===\n")

    # Sync code first to ensure install script is present
    sync_to_hpc()

    # Set SHADE_SOURCE environment variable for the R script
    env_vars = []
    if SHADE_SOURCE:
        env_vars.append(f"export SHADE_SOURCE='{SHADE_SOURCE}'")
    if force_reinstall:
        env_vars.append("export SHADE_FORCE_REINSTALL='true'")
        print("  Force reinstall enabled for SHADE package\n")

    env_var_cmd = " && ".join(env_vars) if env_vars else ""

    # Run installation script on HPC
    install_cmd = (
        f"ssh {HPC_HOST} 'cd {HPC_DIR} && "
        f"{env_var_cmd} && "
        f"module load {HPC_COMPILER_TOOLCHAIN} && "
        f"module load {HPC_R_MODULE} && "
        f"Rscript scripts/install_r_packages.R'"
    )

    print("Running R package installation on HPC...")
    print("This may take 10-30 minutes depending on how many packages need to be installed.\n")

    run_command(install_cmd)

    print("\n✓ R package installation complete\n")


def main():
    parser = argparse.ArgumentParser(
        description="HPC workflow tool for SHADE paper code",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Install R packages on HPC (run once initially)
  python scripts/hpc.py install

  # Force reinstall SHADE package (useful after updates)
  python scripts/hpc.py install --force

  # Sync code to HPC (optional - happens automatically on submit)
  python scripts/hpc.py sync

  # Submit a job (you'll be asked to confirm array_range is correct)
  python scripts/hpc.py submit sim_timing 02_fit_models

  # Force rerun completed simulations
  python scripts/hpc.py submit sim_timing 02_fit_models --force

  # Check job status
  python scripts/hpc.py status
  python scripts/hpc.py status sim_timing

  # Show job summary (pending/running/completed/failed counts)
  python scripts/hpc.py summary
  python scripts/hpc.py summary sim_shade_comparison

  # Check status of specific job batch
  python scripts/hpc.py job-status 13659381

  # Fetch results
  python scripts/hpc.py fetch sim_timing

  # View logs
  python scripts/hpc.py logs sim_timing 02_fit_models
  python scripts/hpc.py logs sim_timing 02_fit_models --lines 50
        """
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Sync command
    subparsers.add_parser("sync", help="Sync code to HPC")

    # Submit command
    submit_parser = subparsers.add_parser("submit", help="Submit a SLURM job")
    submit_parser.add_argument("sim_dir", help="Simulation directory (e.g., sim_timing)")
    submit_parser.add_argument("step", help="Step to run (e.g., 01_generate_data)")
    submit_parser.add_argument("--force", action="store_true",
                               help="Force rerun even if output files already exist")

    # Status command
    status_parser = subparsers.add_parser("status", help="Check job status")
    status_parser.add_argument("sim_dir", nargs="?", help="Filter by simulation directory")

    # Summary command
    summary_parser = subparsers.add_parser("summary", help="Show job summary statistics")
    summary_parser.add_argument("sim_dir", nargs="?", help="Filter by simulation directory")

    # Job batch status command
    job_status_parser = subparsers.add_parser("job-status", help="Show status of a specific job batch")
    job_status_parser.add_argument("job_id", help="Job batch ID (e.g., 13659381)")

    # Fetch command
    fetch_parser = subparsers.add_parser("fetch", help="Fetch results from HPC")
    fetch_parser.add_argument("sim_dir", help="Simulation directory")

    # Logs command
    logs_parser = subparsers.add_parser("logs", help="View log files")
    logs_parser.add_argument("sim_dir", help="Simulation directory")
    logs_parser.add_argument("step", help="Step name")
    logs_parser.add_argument("--lines", type=int, default=20, help="Number of lines to show")

    # Install command
    install_parser = subparsers.add_parser("install", help="Install R packages on HPC")
    install_parser.add_argument("--force", action="store_true",
                                help="Force reinstall SHADE package even if already installed")

    args = parser.parse_args()

    if args.command == "sync":
        sync_to_hpc()
    elif args.command == "submit":
        submit_job(args.sim_dir, args.step, force=args.force)
    elif args.command == "status":
        check_status(args.sim_dir)
    elif args.command == "summary":
        job_summary(args.sim_dir)
    elif args.command == "job-status":
        job_batch_status(args.job_id)
    elif args.command == "fetch":
        fetch_results(args.sim_dir)
    elif args.command == "logs":
        view_logs(args.sim_dir, args.step, args.lines)
    elif args.command == "install":
        install_packages(force_reinstall=args.force)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
