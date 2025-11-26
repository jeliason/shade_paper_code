#!/bin/bash
# Setup Python environment using uv for HPC workflow tools

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
VENV_DIR="$PROJECT_ROOT/.venv"

echo "=== Setting up Python environment for HPC workflow tools ==="
echo

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "Error: uv is not installed."
    echo
    echo "Install uv with:"
    echo "  curl -LsSf https://astral.sh/uv/install.sh | sh"
    echo
    echo "Or on macOS:"
    echo "  brew install uv"
    echo
    exit 1
fi

echo "✓ Found uv: $(uv --version)"
echo

# Create virtual environment if it doesn't exist
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment at $VENV_DIR..."
    uv venv "$VENV_DIR"
    echo "✓ Virtual environment created"
else
    echo "✓ Virtual environment already exists at $VENV_DIR"
fi
echo

# Install dependencies
echo "Installing dependencies from requirements.txt..."
uv pip install -r "$SCRIPT_DIR/requirements.txt"
echo "✓ Dependencies installed"
echo

# Create .env if it doesn't exist
if [ ! -f "$PROJECT_ROOT/.env" ]; then
    echo "Creating .env file from template..."
    cp "$PROJECT_ROOT/.env.example" "$PROJECT_ROOT/.env"
    echo "✓ Created .env (please edit with your HPC details)"
    echo
    echo "⚠️  IMPORTANT: Edit .env with your HPC configuration:"
    echo "   nano $PROJECT_ROOT/.env"
else
    echo "✓ .env file already exists"
fi
echo

echo "=== Setup complete! ==="
echo
echo "To use the HPC tools, activate the environment:"
echo "  source .venv/bin/activate"
echo
echo "Then run commands like:"
echo "  python scripts/hpc.py sync"
echo "  python scripts/hpc.py submit sim_timing 01_generate_data"
echo
echo "Or use uv run to run commands without activating:"
echo "  uv run scripts/hpc.py sync"
echo
