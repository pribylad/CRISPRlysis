# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name

"""File contains general functions used in the CRISPRlysis project."""
import os
from pathlib import Path
import shutil

# TODO check if all functions are using detination_folder
def write_into_report(text_to_append, report_name='report.txt', destination_folder=None):
    """
    Writes text to a report file.

    Parameters:
    text_to_append (str): Text to be appended to the report.
    report_name (str): Name of the report file (default: 'report.txt').
    destination_folder (str): Path where the report should be created (default: current directory).
    """
    try:
        report_path = Path(destination_folder or os.getcwd(), report_name)
        mode = 'a' if os.path.isfile(report_path) else 'w'
        with open(report_path, mode, encoding='utf-8') as file:
            if mode == 'w':
                file.write('########## REPORT ##########\n')
                destination_folder = os.getcwd()
            file.write(text_to_append + '\n\n')
            print(text_to_append)
    except (FileNotFoundError, PermissionError) as e:
        print(f'Error: {e}')

def copy_files_with_parent_folder_to_another(source_dir, project_output_dir, prefix, suffix):
    """
    Copies files matching the criteria along with their parent folders to a new location.

    Parameters:
        source_dir (str): Source directory.
        project_output_dir (str): Destination directory.
        prefix (str): Prefix of files to copy.
        suffix (str): Suffix of files to copy.

    Outputs:
        This function does not return any value. It copies files from the source directory to the project output directory, maintaining the parent folder structure.
    """

    if not Path(source_dir).exists():
        write_into_report(f'Source path {source_dir} does not exists.', destination_folder=project_output_dir)
        return
    if not Path(project_output_dir).exists():
        write_into_report(f'Project output path {project_output_dir} does not exists and it will be created.', destination_folder=project_output_dir)
    counter=0
    try:
        for root, _, files in os.walk(source_dir):
            for file in files:
                if file.startswith(prefix) and file.endswith(suffix):
                    src_path = Path(root, file)
                    rel_path = os.path.relpath(src_path, source_dir)
                    dest_path = Path(project_output_dir, rel_path)
                    os.makedirs(dest_path.parent, exist_ok=True)
                    shutil.copyfile(src_path, dest_path)
                    counter += 1
                    write_into_report(f'Copied {src_path} to {dest_path}')
        write_into_report(f'Copied {counter} files from the {source_dir} to {project_output_dir} folder.', destination_folder=project_output_dir)
        if counter == 0:
            write_into_report(f'Error in copy_files_with_parent_folder_to_another function. Parameters: source_dir={source_dir}, project_output_dir={project_output_dir}, prefix={prefix}, suffix={suffix}', destination_folder=project_output_dir)
    except (FileNotFoundError, PermissionError, shutil.Error) as e:
        write_into_report(f'Error: {e}', destination_folder=project_output_dir)

def find_files_paths_with_criteria(directory, file_prefix, file_suffix):
    """
    Makes list of file paths matching given criteria.

    Parameters:
        project_dir (str): Path where the files should be found.
        file_prefix (str): How the file name should start.
        file_suffix (str): Sufix of the file.

    Returns:
        list: [list, of, files, matching, given, criteria]
    """
    matching_files = []
    for root, _, files in os.walk(directory):
        for filename in files:
            if filename.startswith(file_prefix) and filename.endswith(file_suffix):
                matching_files.append(os.path.join(root, filename))
    return matching_files

def check_if_path_exist_create_if_not(project_output_dir):
    """
    Checks if the path exists, and if not, creates it.

    Parameters:
        project_output_dir (str): Path to the directory.
    """
    output_path = Path(project_output_dir)
    if not output_path.exists():
        output_path.mkdir(parents=True, exist_ok=True)
        print(f'Created folder: {output_path}')
    else:
        print('The project output folder already exists.')
