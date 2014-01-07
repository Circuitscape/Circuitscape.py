import ConfigParser, os, copy, ast, codecs

class CSConfig:
    """Represents a Circuitscape configuration object"""
    
    FILE_PATH_PROPS = ['polygon_file', 'source_file', 'ground_file', 'mask_file', 'output_file', 'habitat_file', 'point_file', 'reclass_file']
    
    DEFAULTS = {
        'Version': {
            'version': 'unknown'
        },
        'Connection scheme for raster habitat data': {
            'connect_four_neighbors_only': False, 
            'connect_using_avg_resistances': False,
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
            'remove_src_or_gnd': 'keepall' 
        }, 
        'Mask file': {
            'mask_file': 'None', 
            'use_mask': False
        }, 
        'Calculation options': {
            'preemptive_memory_release': False,
            'low_memory_mode': False,
            'parallelize': False,          # can parallelize if True. It may be overridden internally to False based on os support and suitability of the problem.
            'max_parallel': 0,            # passing 0 results in using all available cpus
            'print_timings': False, 
            'print_rusages': False, 
            'solver': 'cg+amg'
        }, 
        'Options for one-to-all and all-to-one modes': {
            'use_variable_source_strengths': False, 
            'variable_source_file': 'None'
        }, 
        'Output options': {
            'set_null_currents_to_nodata': False, 
            'output_file': '(Choose a base name for output files)', 
            'write_cum_cur_map_only': False, 
            'log_transform_maps': False, 
            'write_max_cur_maps': False, 
            'compress_grids': False, 
            'set_null_voltages_to_nodata': False, 
            'set_focal_node_currents_to_zero': False, 
            'write_volt_maps': False, 
            'write_cur_maps': False
        }, 
        'Habitat raster or graph': {
            'habitat_map_is_resistances': True,
            'habitat_file': '(Browse for a resistance file)'
        }, 
        'Circuitscape mode': {
            'scenario': 'not entered', 
            'data_type': 'raster'
        }, 
        'Options for pairwise and one-to-all and all-to-one modes': {
            'use_included_pairs': False, 
            'included_pairs_file': '(Browse for a file with pairs to include or exclude)', 
            'point_file': '(Browse for file with locations of focal points or regions)'
        },
        'Options for reclassification of habitat data': {
            'use_reclass_table': False,
            'reclass_file': '(Browse for file with reclassification data)'
        },
        'Logging Options': {
            'profiler_log_file': None,      # file to log timing and rusage profiling results 
            'log_file': None,               # file to log regular log messages
            'log_level': 'INFO',           # one of FATAL, ERROR, WARN, INFO, DEBUG
            'screenprint_log': False        # whether to print logs to console (stdout)
        }
    }
    
    CHECKS_AND_MESSAGES = {
        'scenario':                         'Please choose a scenario',
        'habitat_file':                     'Please choose a resistance file',
        'output_file':                      'Please choose an output file name',
        'point_file':                       'Please choose a focal node file',
        'source_file':                      'Please enter a current source file',
        'ground_file':                      'Ground point file does not exist!',
        'reclass_file':                     'Please choose a file with reclassification data',
        'polygon_file':                     'Please enter a short-circuit region file or uncheck this option in the Options menu'   
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
        try:
            config.read(cfgfile)
        except:
            # try again with utf8 bom markers
            with codecs.open(cfgfile, 'r', encoding='utf_8_sig') as fp:
                config.readfp(fp)
                
        for section in config.sections():
            for item in config.items(section):
                try:
                    self.options[item[0]] = ast.literal_eval(item[1])
                except:
                    self.options[item[0]] = item[1]

    def as_dict(self, rel_to_abs=None):
        result = {}
        for section in CSConfig.DEFAULTS.keys():
            for option in CSConfig.DEFAULTS[section].keys():
                val = self.options[option]
                if option in CSConfig.FILE_PATH_PROPS:
                    if (val == None) or (val == CSConfig.DEFAULTS[section][option]):
                        val = ''
                    elif (not os.path.isabs(val)) and (rel_to_abs != None):
                        val = os.path.join(rel_to_abs, val)
                result[option] = val
        return result
    
    def are_all_paths_relative(self):
        defaults = {}
        for olist in CSConfig.DEFAULTS.values(): 
            defaults.update(olist)
            
        for name in CSConfig.FILE_PATH_PROPS:
            if not ((name in self.options) and (self.options[name] != defaults[name]) and (self.options[name] != None)):
                continue
            
            if os.path.isabs(self.options[name]):
                return False
        return True
    
    def write(self, cfg_filename, is_filename_template=False):
        if is_filename_template:
            out_base, _out_extn = os.path.splitext(cfg_filename)
            cfg_filename = out_base + '.ini'
            out_dir = os.path.split(cfg_filename)[0]
            if (len(out_dir) > 0) and (not os.path.isdir(out_dir)):
                try:
                    os.makedirs(out_dir)
                except:
                    raise RuntimeError('Cannot create output directory: ' + out_dir + '.')    

        config = ConfigParser.ConfigParser()
        for section in CSConfig.DEFAULTS.keys():
            config.add_section(section)
            for option in CSConfig.DEFAULTS[section].keys():
                config.set(section, option, self.options[option])

        with open(cfg_filename, 'w') as f:
            config.write(f)

    def check(self):
        """Checks to make sure sufficient options are passed to Circuitscape to complete a run."""
        
        defaults = {} 
        for olist in CSConfig.DEFAULTS.values(): 
            defaults.update(olist)
        
        # get all the checks to be done
        checks = copy.copy(CSConfig.CHECKS_AND_MESSAGES)
        
        # remove checks that are not required
        if self.options['scenario'] not in ['pairwise', 'one-to-all']:
            del checks['point_file']

        if self.options['scenario'] != 'advanced':
            for key in ['source_file', 'ground_file']:#, 'ground_file_is_resistances']:
                del checks[key]

        if self.options['use_polygons'] == False:
            del checks['polygon_file']
        
        if self.options['use_reclass_table'] == False:
            del checks['reclass_file']
            
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
