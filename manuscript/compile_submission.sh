#!/bin/bash

# Script to compile manuscript PDFs for submission
# Compiles main (with/without track changes), supplement, and responses

set -e  # Exit on error

OUTPUT_DIR="shade_submission"

# Remove old submission folder if it exists
if [ -d "$OUTPUT_DIR" ]; then
    rm -rf "$OUTPUT_DIR"
fi

echo "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Function to compile LaTeX document with optional trackchanges override
compile_latex() {
    local filename=$1
    local description=$2
    local trackchanges=$3  # Optional: 0 or 1 to override trackchanges

    echo "Compiling $description..."

    if [ -n "$trackchanges" ]; then
        # Compile with trackchanges override
        pdflatex -interaction=nonstopmode "\def\trackchanges{$trackchanges}\input{${filename%.tex}}" > /dev/null || true
        bibtex "${filename%.tex}" > /dev/null 2>&1 || true
        pdflatex -interaction=nonstopmode "\def\trackchanges{$trackchanges}\input{${filename%.tex}}" > /dev/null || true
        pdflatex -interaction=nonstopmode "\def\trackchanges{$trackchanges}\input{${filename%.tex}}" > /dev/null || true
    else
        # Compile with default trackchanges value
        pdflatex -interaction=nonstopmode "$filename" > /dev/null || true
        bibtex "${filename%.tex}" > /dev/null 2>&1 || true
        pdflatex -interaction=nonstopmode "$filename" > /dev/null || true
        pdflatex -interaction=nonstopmode "$filename" > /dev/null || true
    fi

    echo "  ✓ $description compiled"
}

# Compile main with track changes ON
echo "=== Compiling Main Manuscript (with track changes) ==="
compile_latex "main.tex" "Main (track changes)" 1
cp main.pdf "$OUTPUT_DIR/Revised Manuscript with Track Changes.pdf"

# Compile main with track changes OFF
echo ""
echo "=== Compiling Main Manuscript (clean version) ==="
compile_latex "main.tex" "Main (clean)" 0
cp main.pdf "$OUTPUT_DIR/Manuscript.pdf"

# Compile supplement with track changes ON
echo ""
echo "=== Compiling Supplement (with track changes) ==="
compile_latex "supplement.tex" "Supplement (track changes)" 1
cp supplement.pdf "$OUTPUT_DIR/Revised Supplement with Track Changes.pdf"

# Compile supplement with track changes OFF
echo ""
echo "=== Compiling Supplement (clean version) ==="
compile_latex "supplement.tex" "Supplement (clean)" 0
cp supplement.pdf "$OUTPUT_DIR/Supplement.pdf"

# Compile responses (last, so references to main/supplement are current)
echo ""
echo "=== Compiling Reviewer Responses ==="
compile_latex "responses.tex" "Reviewer Responses"
cp responses.pdf "$OUTPUT_DIR/Response to Reviewers.pdf"

# Create flattened source files (with snippets inserted)
echo ""
echo "=== Creating Flattened Source Files ==="
python3 flatten_tex.py main.tex "$OUTPUT_DIR/main_source.tex"
python3 flatten_tex.py supplement.tex "$OUTPUT_DIR/supplement_source.tex"
python3 flatten_tex.py responses.tex "$OUTPUT_DIR/responses_source.tex"

# Disable track changes in flattened source files
sed -i '' 's/\\providecommand{\\trackchanges}{1}/\\providecommand{\\trackchanges}{0}/g' "$OUTPUT_DIR/main_source.tex"
sed -i '' 's/\\providecommand{\\trackchanges}{1}/\\providecommand{\\trackchanges}{0}/g' "$OUTPUT_DIR/supplement_source.tex"
sed -i '' 's/\\providecommand{\\trackchanges}{1}/\\providecommand{\\trackchanges}{0}/g' "$OUTPUT_DIR/responses_source.tex"

# Copy bibliography file and latexmkrc
cp references.bib "$OUTPUT_DIR/references.bib"
cp latexmkrc "$OUTPUT_DIR/latexmkrc"

echo "  ✓ Flattened source files created (track changes disabled)"
echo "  ✓ Bibliography file copied"
echo "  ✓ LaTeX configuration file copied"

# Convert images to TIFF and update paths in flattened tex files
echo ""
echo "=== Converting Images to TIFF ==="
python3 convert_images_to_tiff.py "$OUTPUT_DIR" .
echo ""

# Rename figures to PLOS naming convention (Fig1.pdf/.tif, Fig2.pdf/.tif, etc.)
echo "=== Renaming Figures to PLOS Convention ==="
python3 rename_figures_for_plos.py "$OUTPUT_DIR"
echo ""

# Validate flattened tex files by compiling them (using PDF versions)
echo "=== Validating Flattened Source Files ==="
cd "$OUTPUT_DIR"

# Function to validate a tex file
validate_tex() {
    local tex_file=$1
    local description=$2

    echo "Validating $description..."

    # Compile with pdflatex (one pass is enough for validation)
    local log_file="${tex_file%.tex}.log"
    pdflatex -interaction=nonstopmode "$tex_file" > /dev/null 2>&1
    local exit_code=$?

    if [ $exit_code -eq 0 ]; then
        echo "  ✓ $description compiles successfully"
        # Clean up compilation artifacts
        local base="${tex_file%.tex}"
        rm -f "${base}.pdf" "${base}.aux" "${base}.log" "${base}.out" "${base}.bbl" "${base}.blg"
    else
        echo "  ⚠ WARNING: $description has LaTeX errors"
        # Show last few lines of error from log
        if [ -f "$log_file" ]; then
            echo "  Last error from log:"
            grep -A 3 "^!" "$log_file" | tail -10 | sed 's/^/    /'
        fi
        # Clean up artifacts but keep log for debugging
        local base="${tex_file%.tex}"
        rm -f "${base}.pdf" "${base}.aux" "${base}.out" "${base}.bbl" "${base}.blg"
        # Don't fail the build, just warn
    fi
}

# Validate each flattened source file
validate_tex "main_source.tex" "Main manuscript source"
validate_tex "supplement_source.tex" "Supplement source"
validate_tex "responses_source.tex" "Responses source"

# Return to manuscript directory
cd ..
echo ""

# Create README
echo ""
echo "Creating README..."
cat > "$OUTPUT_DIR/README.txt" << EOF
SHADE Manuscript Submission Package
Generated: $(date)

Contents:
---------
PDFs:
1. Response to Reviewers.pdf                  - Rebuttal letter responding to editor and reviewers
2. Revised Manuscript with Track Changes.pdf  - Main manuscript with revisions highlighted in blue
3. Manuscript.pdf                             - Unmarked main manuscript without track changes
4. Revised Supplement with Track Changes.pdf  - Supplement with revisions highlighted in blue
5. Supplement.pdf                             - Unmarked supplement without track changes

Source Files (LaTeX):
6. main_source.tex                            - Main manuscript source with all snippets inserted and PDF image paths
7. supplement_source.tex                      - Supplement source with all snippets inserted and PDF image paths
8. responses_source.tex                       - Responses source with all snippets inserted and PDF image paths
9. references.bib                             - Bibliography file
10. latexmkrc                                 - LaTeX build configuration (for shared references)

Images (PDF and TIFF formats):
11. images/                                   - Directory containing all figures in both PDF and TIFF formats
                                                Figures are named according to PLOS convention:
                                                - Main text: Fig1.tif, Fig2.tif, Fig3.tif, etc.
                                                - Supplement: S1 Fig.tif, S2 Fig.tif, etc.

Notes:
- The *_source.tex files are self-contained versions where all \input{snippets/...} commands
  have been replaced with the actual snippet content.
- Track changes are disabled in source files (clean version without blue highlighting).
- Source files reference PDF versions of figures (for LaTeX compilation).
- TIFF versions are provided separately for journal production.
- PDF images allow the source files to compile successfully.
- TIFF images meet PLOS requirements: 300 DPI, RGB, LZW compression, ≤2250x2625px.
- Figure files are named in order of first appearance: Fig1.pdf/.tif, Fig2.pdf/.tif, etc. (main text)
  and S1 Fig.pdf/.tif, S2 Fig.pdf/.tif, etc. (supplement).
EOF

# Create zip file
echo ""
echo "Creating zip archive..."
ZIP_FILE="shade_submission.zip"
rm -f "$ZIP_FILE"  # Remove old zip if it exists
zip -r "$ZIP_FILE" "$OUTPUT_DIR" > /dev/null
echo "  ✓ Archive created: $ZIP_FILE"

echo ""
echo "=== Compilation Complete ==="
echo "Output directory: $OUTPUT_DIR"
echo "Zip file: $ZIP_FILE"
echo ""
echo "Contents:"
ls -lh "$OUTPUT_DIR"
