#!/usr/bin/env python3
"""
Flatten LaTeX files by replacing \input{snippets/...} with actual content.
Creates self-contained .tex files for journal submission.
"""

import re
import sys
from pathlib import Path

def flatten_tex(input_file, output_file):
    """
    Read a .tex file and replace all \input{snippets/...} with actual content.
    """
    input_path = Path(input_file)
    output_path = Path(output_file)

    if not input_path.exists():
        print(f"Error: {input_file} not found")
        sys.exit(1)

    # Read the input file
    with open(input_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Pattern to match \input{snippets/filename.tex}
    # Handles both with and without .tex extension
    pattern = r'\\input\{snippets/([^}]+)\}'

    def replace_input(match):
        """Replace \input command with file contents."""
        snippet_name = match.group(1)

        # Add .tex extension if not present
        if not snippet_name.endswith('.tex'):
            snippet_name += '.tex'

        snippet_path = input_path.parent / 'snippets' / snippet_name

        if not snippet_path.exists():
            print(f"Warning: Snippet file not found: {snippet_path}")
            return match.group(0)  # Keep original if file not found

        # Read snippet content
        with open(snippet_path, 'r', encoding='utf-8') as f:
            snippet_content = f.read()

        # Add comment markers to show where snippet was inserted
        return (f"% BEGIN: {snippet_name}\n"
                f"{snippet_content}\n"
                f"% END: {snippet_name}")

    # Replace all \input{snippets/...} with actual content
    flattened_content = re.sub(pattern, replace_input, content)

    # Write flattened content to output file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(flattened_content)

    print(f"✓ Flattened {input_file} → {output_file}")

    # Count how many snippets were replaced
    num_snippets = len(re.findall(pattern, content))
    print(f"  Replaced {num_snippets} snippet(s)")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python flatten_tex.py <input.tex> <output.tex>")
        print("\nExample:")
        print("  python flatten_tex.py main.tex main_flattened.tex")
        print("  python flatten_tex.py supplement.tex supplement_flattened.tex")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    flatten_tex(input_file, output_file)
