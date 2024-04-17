import re
import argparse

def extract_extrusion_width(lines):
    """
    Extract extrusion width from G-code lines.

    :param lines: List of all lines in the G-code file.
    :return: Extrusion width as a string without 'mm' or None if not found.
    """
    pattern = r'; external perimeters extrusion width = (\d+\.\d+)mm'
    for line in lines:
        match = re.search(pattern, line)
        if match:
            return match.group(1)  # Return the number only, without 'mm'
    return None

def inject_commands_after_marker(input_file, marker):
    """
    Inject specific commands into a G-code file after a specified marker.

    :param input_file: Path to the original G-code file.
    :param output_file: Path to save the modified G-code file.
    :param marker: The marker line after which commands should be injected.
    :param width: Extrusion width to inject into the file.
    """
    try:
        with open(input_file, 'r') as file:
            lines = file.readlines()

        extrusion_width = extract_extrusion_width(lines)
        if extrusion_width is None:
            print("Extrusion width not found in the file.")
            return

        # Commands to inject, using the extracted width
        commands_to_inject = [
            f'; overhangs = 1\n',
            f'; extrusion_width = {extrusion_width}\n',
            f'; perimeter_extrusion_width = {extrusion_width}\n',
            f'; solid_infill_extrusion_width = {extrusion_width}\n'
        ]

        # Prepare the new lines with the injection
        modified_lines = []
        for line in lines:
            modified_lines.append(line)
            if marker in line:
                # Add the specified commands right after the marker
                modified_lines.extend(commands_to_inject)
        
        output_file = input_file.replace('.', '-processed.')
        # Write the modified content to a new file
        with open(output_file, 'w') as file:
            file.writelines(modified_lines)

        print("Commands injected and file saved as", output_file)

    except FileNotFoundError:
        print("The file was not found:", input_file)
    except Exception as e:
        print("An error occurred:", e)

def parse_args():
    """
    Parse command line arguments using argparse.
    :return: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Inject commands into G-code files.")
    parser.add_argument('input_file', type=str, help='Path to the original G-code file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    marker = '; CONFIG_BLOCK_START'
    inject_commands_after_marker(args.input_file, marker)
