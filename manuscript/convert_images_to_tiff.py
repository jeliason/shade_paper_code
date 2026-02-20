#!/usr/bin/env python3
"""
Convert images to TIFF for journal submission and update tex file references.
"""

import re
import sys
import subprocess
from pathlib import Path

def extract_image_paths(tex_file):
    r"""Extract all image paths from \includegraphics commands."""
    with open(tex_file, 'r', encoding='utf-8') as f:
        content = f.read()

    # Match \includegraphics[...]{path} and \includegraphics{path}
    pattern = r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}'
    matches = re.findall(pattern, content)

    return matches

def convert_to_tiff(image_path, output_dir, base_name):
    """
    Convert an image to TIFF format using ImageMagick.
    Also copies/converts to PDF with the same base name.
    Returns tuple: (pdf_filename, tiff_filename)
    """
    input_path = Path(image_path)

    # Handle images with or without extension
    if not input_path.suffix:
        # Try common extensions
        for ext in ['.pdf', '.png', '.jpg', '.jpeg']:
            candidate = input_path.with_suffix(ext)
            if candidate.exists():
                input_path = candidate
                break

    if not input_path.exists():
        print(f"Warning: Image not found: {input_path}")
        return None, None

    # Create output filenames using the provided base name
    pdf_name = f'{base_name}.pdf'
    tiff_name = f'{base_name}.tif'
    pdf_path = Path(output_dir) / pdf_name
    tiff_path = Path(output_dir) / tiff_name

    # First, copy/convert to PDF (for LaTeX compilation)
    try:
        import shutil
        if input_path.suffix.lower() == '.pdf':
            # Just copy the PDF
            shutil.copy2(input_path, pdf_path)
        else:
            # Convert to PDF
            pdf_cmd = ['magick', str(input_path), str(pdf_path)]
            subprocess.run(pdf_cmd, check=True, capture_output=True, text=True)
    except Exception as e:
        print(f"  ⚠ Warning: Failed to create PDF {pdf_name}: {e}")
        pdf_name = None

    # Convert to TIFF following PLOS guidelines:
    # - 300 DPI resolution
    # - RGB color mode (or grayscale)
    # - LZW compression
    # - Flattened (no alpha channels)
    # - Max dimensions: 2250x2625 pixels at 300 DPI
    # - Auto-resize if exceeds limits
    try:
        cmd = [
            'magick',
            '-density', '300',  # 300-600 DPI required
            str(input_path),
            '-colorspace', 'sRGB',  # RGB color mode
            '-alpha', 'remove',  # Remove alpha channel
            '-flatten',  # Flatten layers
            '-resize', '2250x2625>',  # Resize only if larger, maintain aspect ratio
            '-compress', 'lzw',  # LZW compression required
            '-units', 'PixelsPerInch',
            str(tiff_path)
        ]

        result = subprocess.run(cmd, check=True, capture_output=True, text=True)

        # Check file size (must be < 10 MB)
        file_size_mb = tiff_path.stat().st_size / (1024 * 1024)
        if file_size_mb > 10:
            print(f"  ⚠ Warning: {tiff_name} is {file_size_mb:.1f} MB (exceeds 10 MB limit)")

        # Check dimensions at 300 DPI (max: 2250 x 2625 pixels)
        identify_cmd = ['magick', 'identify', '-format', '%w %h', str(tiff_path)]
        dims = subprocess.run(identify_cmd, capture_output=True, text=True, check=True)
        width, height = map(int, dims.stdout.strip().split())

        print(f"  ✓ Converted: {input_path.name} -> {base_name} (PDF + TIFF, {width}x{height}px, {file_size_mb:.1f}MB)")
        return pdf_name, tiff_name

    except subprocess.CalledProcessError as e:
        print(f"  ✗ Failed to convert {input_path}: {e}")
        return None, None
    except FileNotFoundError:
        print("Error: ImageMagick not found. Please install: brew install imagemagick")
        sys.exit(1)

def update_tex_file(tex_file, image_mapping, output_file):
    """Update image paths in tex file to use PDF versions (for compilation)."""
    with open(tex_file, 'r', encoding='utf-8') as f:
        content = f.read()

    # Replace image paths
    for original_path, (pdf_name, tiff_name) in image_mapping.items():
        if pdf_name:  # Only replace if conversion succeeded
            # Replace with PDF path (in images/ subdirectory) for LaTeX compilation
            replacement = f'images/{pdf_name}'
            content = re.sub(
                r'(\\includegraphics(?:\[[^\]]*\])?\{)' + re.escape(original_path) + r'(\})',
                r'\1' + replacement + r'\2',
                content
            )

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(content)

    print(f"  ✓ Updated {output_file}")

def process_submission(manuscript_dir, output_dir):
    """Process all flattened tex files and convert images."""
    manuscript_dir = Path(manuscript_dir)
    output_dir = Path(output_dir)
    images_dir = output_dir / 'images'
    images_dir.mkdir(parents=True, exist_ok=True)

    # Find all flattened tex files (they're in output_dir, created by flatten_tex.py)
    tex_files = [
        output_dir / 'main_source.tex',
        output_dir / 'supplement_source.tex',
        output_dir / 'responses_source.tex'
    ]

    all_images = set()
    for tex_file in tex_files:
        if tex_file.exists():
            images = extract_image_paths(tex_file)
            all_images.update(images)

    print(f"Found {len(all_images)} unique images to convert")

    # Convert all images
    image_mapping = {}
    for img_path in sorted(all_images):
        # Make path relative to manuscript directory
        full_path = manuscript_dir / img_path
        # Extract base name from path (e.g., "images_foo_bar" from any path)
        relative_name = str(img_path).replace('/', '_').replace('\\', '_')
        base_name = Path(relative_name).stem
        pdf_name, tiff_name = convert_to_tiff(full_path, images_dir, base_name)
        image_mapping[img_path] = (pdf_name, tiff_name)

    # Update tex files
    print("\nUpdating tex files with TIFF paths...")
    for tex_file in tex_files:
        if tex_file.exists():
            output_file = output_dir / tex_file.name
            update_tex_file(tex_file, image_mapping, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_images_to_tiff.py <submission_dir> <manuscript_dir>")
        print("\nExample:")
        print("  python convert_images_to_tiff.py shade_submission .")
        sys.exit(1)

    submission_dir = sys.argv[1]
    manuscript_dir = sys.argv[2]

    process_submission(manuscript_dir, submission_dir)
    print("\n✓ Image conversion complete")
