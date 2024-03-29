import subprocess, os

def configureDoxyfile(input_dir, output_dir):
    with open('Doxyfile.in', 'r') as file :
        filedata = file.read()

    filedata = filedata.replace('@DOXYGEN_INPUT_DIR@', input_dir)
    filedata = filedata.replace('@DOXYGEN_OUTPUT_DIR@', output_dir)

    with open('Doxyfile', 'w') as file:
        file.write(filedata)

# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

breathe_projects = {}

if read_the_docs_build:
    input_dir = '../include'
    output_dir = '_build'
    configureDoxyfile(input_dir, output_dir)
    subprocess.call('doxygen', shell=True)
    breathe_projects['numopt'] = output_dir + '/xml'

#...
extensions = [ "breathe" ]
#...

# Breathe Configuration
breathe_default_project = "numopt"

# Theme
html_theme = "sphinx_rtd_theme"

# Project
project = "numopt"
#author
author = "Niranjan"
#copyright
copyright = "Fly, you fools"
#version
version = "v0.0"
