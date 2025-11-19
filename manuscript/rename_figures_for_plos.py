#!/usr/bin/env python3
"""
Rename figures to PLOS naming convention based on order of appearance.
Main text: Fig1.tif, Fig2.tif, Fig3.tif, etc.
Supplement: S1 Fig.tif, S2 Fig.tif, S3 Fig.tif, etc.
"""

import re
import sys
from pathlib import Path
from collections import OrderedDict


def extract_figures_in_order(tex_file):
    """
    Extract figure paths in order of first appearance.
    Returns list of (image_path, is_main_figure) tuples.
    """
    with open(tex_file, 'r', encoding='utf-8') as f:
        content = f.read()

    # Find all \includegraphics commands with their paths
    pattern = r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}'
    matches = re.findall(pattern, content)

    # Return paths in order, preserving duplicates for now
    figures = []
    seen = set()
    for match in matches:
        if match not in seen:
            figures.append(match)
            seen.add(match)

    return figures


def create_figure_mapping(main_figures, supplement_figures):
    """
    Create mapping from old paths to new PLOS-compliant base names.
    Returns dict: {old_path: new_base_name} (without extension)
    """
    mapping = {}

    # Main text figures: Fig1, Fig2, etc. (will create both .pdf and .tif)
    for i, fig_path in enumerate(main_figures, start=1):
        new_base_name = f'Fig{i}'
        mapping[fig_path] = new_base_name

    # Supplement figures: S1 Fig, S2 Fig, etc.
    for i, fig_path in enumerate(supplement_figures, start=1):
        new_base_name = f'S{i} Fig'
        mapping[fig_path] = new_base_name

    return mapping


def rename_files(images_dir, mapping):
    """
    Rename both PDF and TIFF files according to mapping.
    """
    images_dir = Path(images_dir)

    for old_path, new_base_name in mapping.items():
        # Extract just the filename from the path (e.g., "images/foo.pdf" → "foo")
        # The path format is: images/filename.ext
        old_filename = Path(old_path).name  # Get just the filename with extension
        old_base = Path(old_filename).stem   # Remove extension

        # Rename both PDF and TIFF versions
        for ext in ['.pdf', '.tif']:
            old_file = images_dir / f'{old_base}{ext}'
            new_file = images_dir / f'{new_base_name}{ext}'

            if old_file.exists():
                old_file.rename(new_file)
                if ext == '.pdf':  # Only print once per figure
                    print(f'  ✓ Renamed: {old_base} → {new_base_name} (PDF + TIFF)')
            elif ext == '.tif':  # Warn only if TIFF is missing (PDF might not exist for all)
                print(f'  ⚠ Warning: TIFF not found: {old_file}')


def update_tex_references(tex_file, output_file, mapping):
    """
    Update figure references in tex file to use new PDF names (for compilation).
    """
    with open(tex_file, 'r', encoding='utf-8') as f:
        content = f.read()

    # Replace each old path with new PDF filename
    for old_path, new_base_name in mapping.items():
        # old_path is already in the format: images/filename.pdf
        # We just need to replace it with: images/Fig1.pdf
        new_pdf_path = f'images/{new_base_name}.pdf'

        # Replace the path
        content = re.sub(
            r'(\\includegraphics(?:\[[^\]]*\])?\{)' + re.escape(old_path) + r'(\})',
            r'\1' + new_pdf_path + r'\2',
            content
        )

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(content)

    print(f'  ✓ Updated: {output_file}')


def process_submission(submission_dir):
    """
    Rename figures in submission directory to PLOS convention.
    """
    submission_dir = Path(submission_dir)
    images_dir = submission_dir / 'images'

    main_tex = submission_dir / 'main_source.tex'
    supplement_tex = submission_dir / 'supplement_source.tex'
    responses_tex = submission_dir / 'responses_source.tex'

    # Extract figures in order of appearance
    print("Extracting figure order from manuscripts...")
    main_figures = extract_figures_in_order(main_tex) if main_tex.exists() else []
    supplement_figures = extract_figures_in_order(supplement_tex) if supplement_tex.exists() else []

    print(f"  Main text: {len(main_figures)} figures")
    print(f"  Supplement: {len(supplement_figures)} figures")

    # Create mapping
    print("\nCreating PLOS naming mapping...")
    mapping = create_figure_mapping(main_figures, supplement_figures)

    # Show mapping
    print("\nFigure renaming plan:")
    for old_path, new_base in mapping.items():
        old_filename = Path(old_path).name
        old_base = Path(old_filename).stem
        print(f"  {old_base} → {new_base} (PDF + TIFF)")

    # Rename files
    print("\nRenaming TIFF files...")
    rename_files(images_dir, mapping)

    # Update tex file references
    print("\nUpdating tex file references...")
    for tex_file in [main_tex, supplement_tex, responses_tex]:
        if tex_file.exists():
            # Create mapping for this specific file
            file_figures = extract_figures_in_order(tex_file)
            file_mapping = {fig: mapping[fig] for fig in file_figures if fig in mapping}

            # Update in place
            update_tex_references(tex_file, tex_file, file_mapping)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python rename_figures_for_plos.py <submission_dir>")
        print("\nExample:")
        print("  python rename_figures_for_plos.py shade_submission")
        sys.exit(1)

    submission_dir = sys.argv[1]

    if not Path(submission_dir).exists():
        print(f"Error: Directory not found: {submission_dir}")
        sys.exit(1)

    process_submission(submission_dir)
    print("\n✓ Figure renaming complete")
