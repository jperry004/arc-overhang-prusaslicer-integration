def ensure_width_after_perimeter(input_file, output_file):
    """
    Ensure that every 'overhang perimeter' line is followed by a 'width' setting.

    :param input_file: Path to the original G-code file.
    :param output_file: Path to save the modified G-code file.
    """
    try:
        with open(input_file, 'r') as file:
            lines = file.readlines()

        modified_lines = []
        last_width = "0.40"  # Default width value or a safe assumption if not specified
        add_width_next = False

        for i in range(len(lines)):
            line = lines[i]
            if 'M205' in line or 'asjfdsakjhfd' in line:
                continue
            # Check and update the most recent width observed
            if "WIDTH" in line and "=" in line:
                last_width = line.split('=')[-1].strip()
            
            modified_lines.append(line)

            if 'Overhang perimeter' in line:
                # Check the next line if it exists and if it does not contain a width setting
                if i + 1 < len(lines) and not ("WIDTH" in lines[i + 1]):
                    width_line = f";WIDTH:{last_width}\n"
                    modified_lines.append(width_line)  # Append the width line right after
                    h_line = f";HEIGHT:{last_width}\n"
                    modified_lines.append(h_line)

        # Write the modified content to a new file
        with open(output_file, 'w') as file:
            file.writelines(modified_lines)

        print("File has been modified and saved as", output_file)

    except FileNotFoundError:
        print("The file was not found:", input_file)
    except Exception as e:
        print("An error occurred:", e)

# Example usage of the function
input_path = '/home/josh/repos/arc-overhang-prusaslicer-integration/carl-processed.txt'
output_path = '/home/josh/repos/arc-overhang-prusaslicer-integration/carl-processed-2.gcode'
ensure_width_after_perimeter(input_path, output_path)
