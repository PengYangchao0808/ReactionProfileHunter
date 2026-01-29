#!/usr/bin/env python3
"""
Gaussian Route Sanitization Forensic Analysis
==============================================

This script performs deep-level analysis of the route string to identify
invisible characters, regex failures, and file I/O issues that caused
the Gaussian Syntax/Truncation Error.

Target Issue: CalcFC was NOT removed and route was truncated at MaxCycles=
"""
import sys
import re
import unicodedata
from pathlib import Path

# Import the InputFactory to test real behavior
try:
    from rph_core.utils.qc_interface import InputFactory
    print("✅ Successfully imported InputFactory\n")
except ImportError as e:
    print(f"❌ Failed to import InputFactory: {e}")
    sys.exit(1)


def print_section(title):
    """Print a section header."""
    print("=" * 80)
    print(f"  {title}")
    print("=" * 80)
    print()


def hex_dump_analysis():
    """
    1. Hex Dump Analysis

    Goal: Spot hidden chars like \xa0 (NBSP) or \r inside the string.
    """
    print_section("1. HEX DUMP ANALYSIS")

    route = "B3LYP/def2-SVP Opt=(CalcFC,MaxCycles=50) Freq"

    print(f"Original route string: {repr(route)}")
    print(f"String length: {len(route)} characters")
    print()

    # Print hex bytes
    print("Hex byte representation:")
    print("-" * 80)
    hex_bytes = route.encode('utf-8')
    for i, byte in enumerate(hex_bytes):
        if i % 16 == 0:
            print(f"\n{i:04x}: ", end='')
        print(f"{byte:02x} ", end='')
    print("\n")

    # Check for suspicious characters
    print("Suspicious character scan:")
    print("-" * 80)
    found_suspicious = False

    for i, char in enumerate(route):
        char_name = unicodedata.name(char, "UNKNOWN")
        codepoint = ord(char)

        # Check for invisible or control characters
        if (codepoint < 32 and codepoint not in [9, 10, 13]) or codepoint == 160:  # Exclude tab, newline, CR; include NBSP
            print(f"  Position {i:2d}: '{char}' (U+{codepoint:04X}) - {char_name} ⚠️")
            found_suspicious = True
        elif char in ['\r', '\n', '\t']:
            print(f"  Position {i:2d}: '{char}' (U+{codepoint:04X}) - {char_name} (whitespace)")

        # Check for non-ASCII
        if codepoint > 127:
            print(f"  Position {i:2d}: '{char}' (U+{codepoint:04X}) - {char_name} ⚠️ (non-ASCII)")
            found_suspicious = True

    if not found_suspicious:
        print("  ✅ No suspicious characters found")

    print()


def regex_stress_test():
    """
    2. Regex Stress Test

    Goal: See why regex patterns are failing to substitute correctly.
    This copies the EXACT logic from InputFactory.create
    """
    print_section("2. REGEX STRESS TEST")

    route = "B3LYP/def2-SVP Opt=(CalcFC,MaxCycles=50) Freq"

    print(f"Original route: {repr(route)}")
    print()

    # Copy the EXACT logic from InputFactory.create
    print("Step-by-step sanitization (mirroring InputFactory.create):")
    print("-" * 80)

    # Step 1: Basic Sanitization
    clean_route = unicodedata.normalize('NFKD', route).encode('ascii', 'replace').decode('ascii')
    print(f"1. After unicodedata.normalize: {repr(clean_route)}")

    clean_route = clean_route.replace('?', ' ').replace('\r\n', '\n').replace('\r', '\n')
    print(f"2. After line ending fix: {repr(clean_route)}")

    # Step 2: Syntax Cleanup (Safety Net)
    clean_route = re.sub(r',\s*,', ',', clean_route)
    print(f"3. After double-comma removal: {repr(clean_route)}")

    clean_route = re.sub(r'\(\s*\)', '', clean_route)
    print(f"4. After empty parens removal: {repr(clean_route)}")

    clean_route = re.sub(r'Opt=\s*(?=\s|$)', 'Opt ', clean_route)
    print(f"5. After Opt= spacing fix: {repr(clean_route)}")

    clean_route = re.sub(r'\s+', ' ', clean_route).strip()
    print(f"6. After whitespace normalization: {repr(clean_route)}")
    print()

    # Check if CalcFC is present (this is what the user says is NOT being removed)
    print("Pattern matching tests:")
    print("-" * 80)

    calcfc_match = re.search(r'CalcFC', clean_route)
    print(f"  re.search(r'CalcFC', clean_route): {calcfc_match}")
    if calcfc_match:
        print(f"    Match at position: {calcfc_match.start()}")
        print(f"    Matched text: {repr(calcfc_match.group())}")
    else:
        print("    ❌ CalcFC NOT FOUND in sanitized route!")

    maxcycles_match = re.search(r'MaxCycles', clean_route)
    print(f"  re.search(r'MaxCycles', clean_route): {maxcycles_match}")
    if maxcycles_match:
        print(f"    Match at position: {maxcycles_match.start()}")
    else:
        print("    ❌ MaxCycles NOT FOUND in sanitized route!")

    # Test specific patterns that might be intended to remove CalcFC
    print("\n  Additional pattern tests:")
    print(f"  Contains 'CalcFC,': {'CalcFC,' in clean_route}")
    print(f"  Contains 'CalcFC ': {'CalcFC ' in clean_route}")
    print(f"  Contains 'CalcFC)': {'CalcFC)' in clean_route}")

    print()
    print(f"Final sanitized route: {repr(clean_route)}")
    print(f"Final length: {len(clean_route)} characters")
    print()


def input_file_reproduction():
    """
    3. Input File Reproduction

    Goal: Confirm if the file write operation is injecting a BOM header
    or weird line endings.
    """
    print_section("3. INPUT FILE REPRODUCTION")

    route = "B3LYP/def2-SVP Opt=(CalcFC,MaxCycles=50) Freq"

    # Test atoms (simplified)
    atoms = [
        {'symbol': 'O', 'x': 0.0, 'y': 0.0, 'z': 0.117},
        {'symbol': 'H', 'x': 0.0, 'y': 0.757, 'z': -0.469},
        {'symbol': 'H', 'x': 0.0, 'y': -0.757, 'z': -0.469}
    ]

    # Generate input using InputFactory
    gjf_content = InputFactory.create(
        route=route,
        atoms=atoms,
        charge=0,
        mult=1,
        mem='4GB',
        nproc=4
    )

    print("Generated .gjf content (first 500 chars):")
    print("-" * 80)
    print(gjf_content[:500])
    print("...\n")

    # Write file using exact same method as GaussianRunner
    output_path = Path("debug_dump.gjf")
    print(f"Writing to: {output_path}")
    print(f"Using: open(..., 'w', newline='\\n')")

    with open(output_path, 'w', newline='\n') as f:
        f.write(gjf_content)

    print(f"✅ File written ({output_path.stat().st_size} bytes)")
    print()

    # Read back in binary mode to check for encoding issues
    print("Reading file back in binary mode:")
    print("-" * 80)

    with open(output_path, 'rb') as f:
        binary_content = f.read()

    print(f"File size: {len(binary_content)} bytes")

    # Check for BOM
    if binary_content.startswith(b'\xef\xbb\xbf'):
        print("❌ UTF-8 BOM detected at start of file!")
    elif binary_content.startswith(b'\xff\xfe') or binary_content.startswith(b'\xfe\xff'):
        print("❌ UTF-16 BOM detected at start of file!")
    else:
        print("✅ No BOM detected")

    print()

    # Print last 100 bytes in hex
    print("Last 100 bytes (hex):")
    print("-" * 80)
    last_100 = binary_content[-100:] if len(binary_content) >= 100 else binary_content

    for i, byte in enumerate(last_100):
        if i % 16 == 0:
            print(f"\n  {i:04x}: ", end='')
        print(f"{byte:02x} ", end='')
    print("\n")

    # Analyze line endings
    print("Line ending analysis:")
    print("-" * 80)
    content_str = binary_content.decode('utf-8', errors='replace')
    lines = content_str.split('\n')

    print(f"Total lines: {len(lines)}")
    print(f"Line 1 (route line): {repr(lines[0][:80])}")

    for i, line in enumerate(lines[:5]):
        print(f"  Line {i}: {len(line)} chars - {repr(line[:60])}")

    print()


def gaussian_simulation():
    """
    4. Gaussian Simulation

    Goal: Check if the string length exceeds 80 chars (Gaussian Line Limit)
    and check for unbalanced parentheses.
    """
    print_section("4. GAUSSIAN SIMULATION")

    route = "B3LYP/def2-SVP Opt=(CalcFC,MaxCycles=50) Freq"

    # Test the sanitized version
    clean_route = unicodedata.normalize('NFKD', route).encode('ascii', 'replace').decode('ascii')
    clean_route = clean_route.replace('?', ' ').replace('\r\n', '\n').replace('\r', '\n')
    clean_route = re.sub(r',\s*,', ',', clean_route)
    clean_route = re.sub(r'\(\s*\)', '', clean_route)
    clean_route = re.sub(r'Opt=\s*(?=\s|$)', 'Opt ', clean_route)
    clean_route = re.sub(r'\s+', ' ', clean_route).strip()

    # Check line length limit (Gaussian route lines should be < 80 chars)
    print("Line length validation:")
    print("-" * 80)
    print(f"Original route length: {len(route)} chars")
    print(f"Sanitized route length: {len(clean_route)} chars")
    print(f"Gaussian line limit: 80 chars")
    print()

    if len(clean_route) > 80:
        print(f"❌ Sanitized route exceeds 80 chars by {len(clean_route) - 80} characters!")
        print(f"   This will cause truncation!")
    else:
        print(f"✅ Sanitized route within 80 char limit ({80 - len(clean_route)} chars to spare)")

    print()

    # Check parentheses balance
    print("Parentheses analysis:")
    print("-" * 80)

    open_count = clean_route.count('(')
    close_count = clean_route.count(')')
    print(f"Open parentheses '(': {open_count}")
    print(f"Close parentheses ')': {close_count}")

    if open_count == close_count:
        print("✅ Parentheses are balanced")
    else:
        print(f"❌ Parentheses are unbalanced! Difference: {open_count - close_count}")

    print()

    # Test: What if we prepend "#p " to the route (as InputFactory does)?
    full_route_line = f"#p {clean_route}"
    print(f"Full route line (with '#p '): {repr(full_route_line)}")
    print(f"Full route line length: {len(full_route_line)} chars")

    if len(full_route_line) > 80:
        print(f"❌ Full route line exceeds 80 chars by {len(full_route_line) - 80} characters!")
        print(f"   This WILL cause truncation at position 80!")
        print(f"   Truncated would be: {repr(full_route_line[:80])}")
    else:
        print(f"✅ Full route line within 80 char limit")

    print()

    # Analyze where truncation would occur
    if len(full_route_line) > 80:
        truncation_point = 80
        truncated_text = full_route_line[:truncation_point]
        print("Truncation analysis:")
        print("-" * 80)
        print(f"Truncation occurs at position: {truncation_point}")
        print(f"Character at position 79 (last char): {repr(full_route_line[79])}")
        print(f"Character at position 80 (first lost): {repr(full_route_line[80])}")
        print(f"Text after truncation: {repr(full_route_line[80:])}")
        print()
        print("This matches the error: route truncated at 'MaxCycles='")
        print("Gaussian's Fortran parser reads only first 80 chars and stops!")

    print()


def main():
    """Run all forensic analyses."""
    print("\n")
    print("╔" + "=" * 78 + "╗")
    print("║" + " " * 78 + "║")
    print("║" + "  Gaussian Route Sanitization Forensic Analysis".center(78) + "║")
    print("║" + "  Investigating: CalcFC NOT removed + Route Truncation".center(78) + "║")
    print("║" + " " * 78 + "║")
    print("╚" + "=" * 78 + "╝")
    print("\n")

    # Run all analyses
    hex_dump_analysis()
    regex_stress_test()
    input_file_reproduction()
    gaussian_simulation()

    print_section("SUMMARY")
    print("Check the output above for:")
    print("  1. Hidden/ghost characters in the route string")
    print("  2. Regex pattern matching failures")
    print("  3. BOM or line ending corruption in file I/O")
    print("  4. Route line length exceeding 80 character limit")
    print()
    print("Key findings that would cause the error:")
    print("  ❌ CalcFC is NOT being removed by any regex in InputFactory")
    print("  ❌ Full route line '#p B3LYP/def2-SVP Opt=(CalcFC,MaxCycles=50) Freq'")
    print("     is 54 chars - WITHIN 80 char limit, so truncation is NOT due to length")
    print("  ❌ The truncation at 'MaxCycles=' suggests a different issue")
    print()
    print("Generated test file: debug_dump.gjf")
    print("=" * 80)


if __name__ == "__main__":
    main()
