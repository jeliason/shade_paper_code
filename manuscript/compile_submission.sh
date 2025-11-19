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
echo "  ✓ Flattened source files created"

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
6. main_source.tex                            - Main manuscript source with all snippets inserted
7. supplement_source.tex                      - Supplement source with all snippets inserted
8. responses_source.tex                       - Responses source with all snippets inserted

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
