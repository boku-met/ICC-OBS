# coding=utf-8

import sys, traceback, warnings, os
warnings.filterwarnings("error")

try: 
    import distributed
except Exception:
    ex_type, ex, tb = sys.exc_info()
    traces = traceback.extract_tb(tb)
    fnames = [x[0] for x in traces]
    for fn in fnames:
        if "distributed/config.py" in fn:
            os.remove(fn)
            curfile = os.path.dirname(os.path.abspath(__file__))+"/lib/.config.py"
            os.system("cp "+curfile+" "+fn)
warnings.filterwarnings("ignore")

import csv

from gooey import Gooey
from gooey import GooeyParser

from lib.ICC_OBS_Main import main_function

'''
advanced: define interpolation method and parameters

falls möglich: ein button mit advanced settings, der zusätzliche felder öffnet
(oder ausgegraute felder ausfüllbar macht, wenn das einfacher ist)

zusätzliche auswahl zur interpolationsmethode (interp_method)
auswahlmöglichkeiten: inverse distance (idw) oder kriging
--> default werte (vorausgefüllt) siehe parameter der main funktion

bei wahl von idw können dann zusätzlich die parameter
    - distance (dist),
    - weight function (weights) und
    - lower limit of weights (low_lim),
bei wahl von kriging kann zusätzlich der parameter
    - variogram_model
gewählt werden
'''

defaults = {
    'lat_min': 43.2,
    'lat_max': 45,
    'lon_min': 15.7,
    'lon_max': 18.8,
    'param': 'pr',
    'start_y': 1981,
    'end_y': 2010,
    'fil_topo': '/path/to/topo_file',
    'dir_stations': '/directory/of/stationdata_files/',
    'fn_metadata': 'name_of_station_metadata',
    'fil_obsgrid': '/path/to/gridded/observations',
    'dir_save': '/directory/to/save/data/',
    'dir_save_name': 'new_data',
    'model_hist': '/hp5/Climaproof/MOD/FINAL/01DEG/ORIGINAL/pr_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_1981-2010_original.nc',
    'model_rcp': '/hp5/Climaproof/MOD/FINAL/01DEG/ORIGINAL/pr_MOHC-HadGEM2-ES_rcp85_r1i1p1_CLMcom-CCLM4-8-17_v1_2011-2100_original.nc',
    'interp_method': 'idw',
    'dist': 100,
    'neighbor': '3',
    'variogram_model': 'gaussian',
}


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)


nonbuffered_stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
sys.stdout = nonbuffered_stdout


@Gooey(
    program_name='ICC-OBS Version 0.1',
    show_stop_warning=True,
    force_stop_is_error=False,
    disable_progress_bar_animation=True,
    required_cols=1,
    default_size=(1080, 1080),
    advanced=True,
    navigation='TABBED',
    image_dir=resource_path('images'),
    menu=[{
        'name': 'File',
        'items': [{
            'type': 'AboutDialog',
            'menuTitle': 'About',
            'name': 'ICC-OBS',
            'description': 'This Tool was developed by BOKU-MET in the course of the CLIMAPROOF project, funded by the Austrian Development Agency (ADA) and co-funded by the United Nations Environmental Programme (UNEP)',
            'version': '0.1 (Beta)',
            'copyright': '(C) 2019',
            'website': 'https://github.com/chriskiehl/Gooey',
            'developer': 'Maria Wind\nInstitute of Meteorology\nUniversity of Natural Resources and Life Sciences, Vienna, Austria\n\nIgor Skoric\nSkoric IT, Graz, Austria',
            'license': 'cc-by'
        }],
    }, {
        'name': 'Help',
        'items': [{
            'type': 'Link',
            'menuTitle': 'Documentation',
            'url': 'https://data.ccca.ac.at/group/climaproof'
        }]
    }]

)
def main():
    # init arguments parser / gui generator
    parser = GooeyParser(
        description="ICC-OBS: A tool for Improving bias-corrected Climate Change scenarios with local OBServational data")
    subparsers = parser.add_subparsers(help='commands', dest='command')
    subparser_required = subparsers.add_parser('Basic', help='Uses default interpolation')
    subparser_advanced = subparsers.add_parser('Advanced', help='Possibility to modify interpolation settings')

    #    # add tab with contact information (about)
    #    subparser_about = subparsers.add_parser('About', help='Information about the tool')
    #    contact_group = subparser_about.add_argument_group(
    #			"Contact Information",
    #			gooey_options={'show_border': 1, 'columns': 1 })
    #    contact_group.add_argument('--BOKU-MET', help='University of Natural Resources and Life Sciences, Vienna')

    argument_groups_required = {}
    argument_groups_advanced = {}
    with open(resource_path('config/groups.csv')) as csv_params:
        reader = csv.DictReader(csv_params)

        for row in reader:
            argument_groups_advanced[row['id']] = subparser_advanced.add_argument_group(
                row['title'],
                row['description'],
                gooey_options={'show_border': row['show_border'] == 1, 'columns': int(row['columns'])})
            if row['id'] != '4':
                argument_groups_required[row['id']] = subparser_required.add_argument_group(
                    row['title'],
                    row['description'],
                    gooey_options={'show_border': row['show_border'] == 1, 'columns': int(row['columns'])})

    with open(resource_path('config/params.csv')) as csv_params:
        reader = csv.DictReader(csv_params)
        for row in reader:
            keyed_arguments = {
                'help': row['description'],
                'default': defaults[row['param']],
                'metavar': row['title']
            }
            if row['choices']:
                keyed_arguments['choices'] = row['choices'].split(',')
            if row['type'] == 'file':
                keyed_arguments['widget'] = 'FileChooser'
            if row['type'] == 'dir':
                keyed_arguments['widget'] = 'DirChooser'
            if row['type'] == 'float':
                keyed_arguments['type'] = float
            if row['type'] == 'int':
                keyed_arguments['type'] = int
            argument_groups_advanced[row['group']].add_argument('--' + row['param'], **keyed_arguments)

            if row['group'] != '4':
                argument_groups_required[row['group']].add_argument('--' + row['param'], **keyed_arguments)

    args = parser.parse_args()

    vars_of_args = vars(args)
    command = vars_of_args.pop('command')

    ordered_args = [vars_of_args['fil_topo'], vars_of_args['lat_min'], vars_of_args['lat_max'],
                    vars_of_args['lon_min'], vars_of_args['lon_max'], vars_of_args['param'], vars_of_args['start_y'],
                    vars_of_args['end_y'], vars_of_args['dir_stations'], vars_of_args['fn_metadata'],
                    vars_of_args['fil_obsgrid'], vars_of_args['model_hist'], vars_of_args['model_rcp'],
                    vars_of_args['dir_save'], vars_of_args['dir_save_name']]

    if command == 'Basic':
        main_function(*ordered_args)
    if command == 'Advanced':
        keyed_args = {
            'interp_method': vars_of_args['interp_method'],
            'distance': vars_of_args['dist'],
            'neighbor': vars_of_args['neighbor'],
            'variogram_model': vars_of_args['variogram_model'],
        }
        main_function(*ordered_args, **keyed_args)
        print(keyed_args)
    pass


if __name__ == "__main__":
    main()
