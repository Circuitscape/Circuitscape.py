##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import ConfigParser, os, string, copy, ast

class CSConfig:
    """Represents a Circuitscape configuration object"""
    
    DEFAULTS = {
        'Version': {
            'version': 'unknown'
        },
        'Connection scheme for raster habitat data': {
            'connect_four_neighbors_only': True, 
            'connect_using_avg_resistances': True,
        },
        'Short circuit regions (aka polygons)': {
            'use_polygons': False,
            'polygon_file': '(Browse for a short-circuit region file)'
        },                  
        'Options for advanced mode': {
            'source_file': '(Browse for a current source file)', 
            'ground_file': '(Browse for a ground point file)', 
            'ground_file_is_resistances': True, 
            'use_unit_currents': False, 
            'use_direct_grounds': False,
            'remove_src_or_gnd': 'not entered' 
        }, 
        'Mask file': {
            'mask_file': 'None', 
            'use_mask': False
        }, 
        'Calculation options': {
            'low_memory_mode': False, 
            'print_timings': False, 
            'solver': 'cg+amg'
        }, 
        'Options for one-to-all and all-to-one modes': {
            'use_variable_source_strengths': False, 
            'variable_source_file': 'None'
        }, 
        'Output options': {
            'set_null_currents_to_nodata': True, 
            'output_file': '(Choose a base name for output files)', 
            'write_cum_cur_map_only': False, 
            'log_transform_maps': False, 
            'write_max_cur_maps': False, 
            'compress_grids': False, 
            'set_null_voltages_to_nodata': True, 
            'set_focal_node_currents_to_zero': False, 
            'write_volt_maps': False, 
            'write_cur_maps': False
        }, 
        'Habitat raster or graph': {
            'habitat_map_is_resistances': True, 
            'habitat_file': '(Browse for a habitat map file)'
        }, 
        'Circuitscape mode': {
            'scenario': 'not entered', 
            'data_type': 'raster'
        }, 
        'Options for pairwise and one-to-all and all-to-one modes': {
            'use_included_pairs': False, 
            'included_pairs_file': 'None', 
            'point_file_contains_polygons': False, 
            'point_file': '(Browse for file with locations of focal points or areas)'
        }
    }
    
    CHECKS_AND_MESSAGES = {
        'scenario':                         'Please choose a scenario',
        'habitat_file':                     'Please choose a raster habitat map file',
        'habitat_map_is_resistances':       'Please choose a habitat data type', 
        'connect_four_neighbors_only':      'Please choose a cell connection scheme',
        'connect_using_avg_resistances':    'Please choose a cell connection calculation', 
        'output_file':                      'Please choose an output file name',
        'point_file':                       'Please choose a focal node file',
        'source_file':                      'Please enter a current source file',
        'ground_file':                      'Ground point file does not exist!',
        'ground_file_is_resistances':       'Please choose a ground data type'            
    }

    def __init__(self, cfgfile=None):
        o = {}
        for olist in CSConfig.DEFAULTS.values(): 
            o.update(olist)
        self.options = o
        
        if None == cfgfile:
            return
        
        if not os.path.isfile(cfgfile):
            raise RuntimeError('File %s does not exist'%(cfgfile,))
    
        config = ConfigParser.ConfigParser()
        config.read(cfgfile)
        for section in config.sections():
            for item in config.items(section):
                try:
                    self.options[item[0]] = ast.literal_eval(item[1])
                except:
                    self.options[item[0]] = item[1]
 
    def write(self, cfgfile):
        config = ConfigParser.ConfigParser()
        for section in CSConfig.DEFAULTS.keys():
            config.add_section(section)
            for option in CSConfig.DEFAULTS[section].keys():
                config.set(section, option, self.options[option])

        with open(cfgfile, 'w') as f:
            config.write(f)

    def check(self):
        """Checks to make sure sufficient options are passed to Circuitscape to complete a run."""
        
        defaults = {} 
        for olist in CSConfig.DEFAULTS.values(): 
            defaults.update(olist)
        
        # get all the checks to be done
        checks = copy.copy(CSConfig.CHECKS_AND_MESSAGES)
        
        # remove checks that are not required
        if self.options['scenario'] in ['pairwise', 'one-to-all']:
            del checks['point_file']
        elif self.options['scenario'] == 'advanced':
            for key in ['source_file', 'ground_file', 'ground_file_is_resistances']:
                del checks[key]

        if self.options['use_polygons'] == True:
            del checks['polygon_file']
        
        # check if values have been entered for the options
        for name in checks.keys(): 
            if self.options[name] == defaults[name]:
                return False,checks[name]
        return True,'None'

    def __getattr__(self, name):
        return self.options[name]
    
    def __setattr__(self, name, value):
        if name == 'options':
            self.__dict__[name] = value
        else:
            self.options[name] = value
        return value
