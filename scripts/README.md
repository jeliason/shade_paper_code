# HPC Workflow Tools

Automated tools for syncing code to HPC, submitting SLURM jobs, and fetching results.

## Setup (One-Time)

### 1. Install uv (if not already installed)

[uv](https://github.com/astral-sh/uv) is a fast Python package manager.

```bash
# macOS
brew install uv

# Linux / macOS (using installer)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows (PowerShell)
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"

# Windows (using winget)
winget install --id=astral-sh.uv -e

# Windows (using scoop)
scoop install uv
```

### 2. Run the setup script

**Linux/macOS:**
```bash
./scripts/setup_env.sh
```

**Windows (PowerShell):**
```powershell
# Run these commands manually
uv venv .venv
uv pip install -r scripts/requirements.txt
Copy-Item .env.example .env
# Then edit .env with your HPC details
notepad .env
```

This will:
- Create a Python virtual environment at `.venv/`
- Install required dependencies (python-dotenv, jinja2, pyyaml)
- Create `.env` from `.env.example` (if it doesn't exist)

**Linux/macOS users** - edit `.env` with your HPC details:
```bash
nano .env
```

### 3. Configure SSH

Ensure your `~/.ssh/config` has an HPC entry:

```
Host hpc
    HostName your-hpc-cluster.edu
    User your_username
    IdentityFile ~/.ssh/id_rsa
```

Test with: `ssh hpc`

### 4. Install R packages on HPC (one-time)

Before running any simulations, install all required R packages on the HPC:

```bash
# Activate Python environment
source .venv/bin/activate  # or use uv run

# Install R packages on HPC
python scripts/hpc.py install
```

This will:
- Sync the code to HPC
- Run `scripts/install_r_packages.R` on the HPC
- Install all CRAN packages (tidyverse, spatstat, cmdstanr, etc.)
- Install CmdStan if not already present
- Install the SHADE package from the source specified in `.env`

**Note:** Make sure you set `SHADE_SOURCE` in your `.env` file:
```bash
# For GitHub installation
SHADE_SOURCE=github:yourusername/SHADE

# For local installation (if SHADE is in your HPC home directory)
SHADE_SOURCE=/home/username/path/to/SHADE
```

The installation may take 10-30 minutes depending on how many packages need to be installed.

## Usage

**Activate the virtual environment** (recommended):

**Linux/macOS:**
```bash
source .venv/bin/activate
```

**Windows (PowerShell):**
```powershell
.venv\Scripts\Activate.ps1
```

**Windows (Command Prompt):**
```cmd
.venv\Scripts\activate.bat
```

Or use `uv run` to run commands without activating (prepend `uv run` to all commands below, works on all platforms).

### Sync code to HPC

```bash
python scripts/hpc.py sync
# Or: uv run scripts/hpc.py sync
```

Syncs the entire project directory to HPC, excluding patterns specified in `.env`.

### Submit a job

```bash
python scripts/hpc.py submit <sim_dir> <step>
```

Examples:
```bash
# Submit data generation for timing sim
python scripts/hpc.py submit sim_timing 01_generate_data

# Submit model fitting
python scripts/hpc.py submit sim_timing 02_fit_models

# Submit analysis
python scripts/hpc.py submit sim_timing 03_analyze_results
```

This will:
1. Generate a `.slurm` file from the template using `slurm_config.yaml`
2. Sync code to HPC (including the generated `.slurm` file)
3. Create `logs/` directory if needed
4. Submit the job via `sbatch`

### Check job status

```bash
# Check all your jobs
python scripts/hpc.py status

# Check jobs for a specific simulation
python scripts/hpc.py status sim_timing
```

### View logs

```bash
# View last 20 lines of the most recent log
python scripts/hpc.py logs sim_timing 02_fit_models

# View more lines
python scripts/hpc.py logs sim_timing 02_fit_models --lines 50
```

### Fetch results

```bash
python scripts/hpc.py fetch sim_timing
```

Downloads:
- `sim_timing/data/*` (all result files)
- `sim_timing/*figures/*` (if present)

## Configuration

### Per-Simulation Config Files

Each simulation directory has a `slurm_config.yaml` defining resources for each step:

```yaml
# sim_timing/slurm_config.yaml
jobs:
  01_generate_data:
    job_name: timing_gen_data
    array_range: auto     # Auto-computed from R script (or use "1-140")
    time: "1:00:00"
    mem: "32G"
    cpus: 1

  02_fit_models:
    job_name: timing_fit_models
    array_range: auto     # Auto-computed from R script (or use "1-140")
    time: "4:00:00"
    mem: "64G"
    cpus: 4
```

### Template System

SLURM scripts are generated on-the-fly from:
- `scripts/slurm_template.sh` (Jinja2 template)
- `<sim_dir>/slurm_config.yaml` (resource specifications)
- `.env` (HPC connection details)

Generated `.slurm` files are gitignored.

### Auto Array Range

Instead of manually specifying array ranges, you can use `array_range: auto` in your `slurm_config.yaml`. The tool will automatically compute the array size from your R script.

To enable this, add a comment in your R script:
```r
# SLURM_ARRAY_SIZE: nrow(grid)
grid <- expand.grid(
  param1 = c(1, 2, 3),
  param2 = 1:10
)
```

The tool will:
1. Find the `# SLURM_ARRAY_SIZE:` comment
2. Extract surrounding context (grid definition)
3. Run `Rscript` to evaluate the expression
4. Generate `array_range: 1-N` where N is the computed size

This keeps your grid definitions as the single source of truth.

## Typical Workflow

```bash
# 0. Activate environment (if not already active)
source .venv/bin/activate              # Linux/macOS
# Or: .venv\Scripts\Activate.ps1      # Windows PowerShell
# Or: .venv\Scripts\activate.bat      # Windows CMD

# 1. Install R packages (first time only)
python scripts/hpc.py install

# 2. Sync latest code
python scripts/hpc.py sync

# 3. Submit jobs in sequence
python scripts/hpc.py submit sim_timing 01_generate_data
# Wait for completion or check status
python scripts/hpc.py status sim_timing

# 4. Submit next step
python scripts/hpc.py submit sim_timing 02_fit_models

# 5. Monitor progress
python scripts/hpc.py logs sim_timing 02_fit_models
python scripts/hpc.py status

# 6. Fetch results when done
python scripts/hpc.py fetch sim_timing
```

## R Package Dependencies

The `scripts/install_r_packages.R` script installs the following packages:

**CRAN packages:**
- tidyverse (data manipulation & plotting)
- Matrix (sparse matrices)
- spatstat (spatial statistics)
- cmdstanr (Stan interface)
- posterior (posterior analysis)
- ggdist (distribution visualizations)
- survival (survival analysis)
- latex2exp (LaTeX in plots)
- patchwork (combine plots)
- splines2 (spline functions)
- survminer (survival plot visualization)
- kableExtra (table formatting)
- tictoc (timing utilities)

**Additional installations:**
- CmdStan (Stan compiler, installed automatically via cmdstanr)
- SHADE package (from GitHub or local path specified in `.env`)

## Adding New Simulations

1. Create simulation directory with R scripts (01_*.R, 02_*.R, etc.)
2. Create `slurm_config.yaml` in the sim directory
3. Run `python scripts/hpc.py submit <sim_dir> <step>`

The tool will automatically generate `.slurm` files from the config.

## Troubleshooting

**"Missing required environment variables"**
- Check your `.env` file has all required fields

**"Config file not found"**
- Create `slurm_config.yaml` in the simulation directory

**SSH connection fails**
- Test: `ssh hpc` manually
- Check `~/.ssh/config` and `HPC_HOST` in `.env`

**Jobs fail immediately**
- Check logs: `python scripts/hpc.py logs <sim_dir> <step>`
- Verify R module name in `.env` matches HPC
- Check resource requests aren't too high

**R package installation fails**
- SSH into HPC and check R version: `module load R && R --version`
- Verify `SHADE_SOURCE` in `.env` points to correct location
- Try installing packages manually on HPC to see specific errors
- For SHADE: ensure GitHub repo is accessible or local path exists on HPC
