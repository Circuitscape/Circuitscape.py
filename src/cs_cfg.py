##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import ConfigParser, os, string

    

# TODO: make config object
def readConfigFile(configFile):
    """Loads .ini file with run options."""    
    if os.path.isfile(configFile)==False:
        raise RuntimeError('File "'  + configFile + '" does not exist')

    config = ConfigParser.ConfigParser()
    config.read(configFile)
    options={}
    
    options['set_null_voltages_to_nodata']=True 
    options['set_null_currents_to_nodata']=True 
    options['set_focal_node_currents_to_zero']=False
    options['write_max_cur_maps']=False #fixme: need to implement for network
    options['low_memory_mode']=False        
    options['use_mask']=False
    options['mask_file']='None' 
    options['use_included_pairs']=False
    options['included_pairs_file']='None' 
    options['use_variable_source_strengths']=False
    options['variable_source_file']='None' 
    options['data_type']='raster' 
    options['version']='unknown'
        
    for section in config.sections():
        for option in config.options(section):
            try:
                options[option]=config.getboolean(section, option)
            except:
                options[option]=config.get(section, option)
    return options

    
def writeConfigFile(configFile, options):
    """Saves .ini file with run options."""    
    config = ConfigParser.ConfigParser()
 
    sections={}
    section='Version'
    sections['version']=section
    
    section='Connection scheme for raster habitat data'
    sections['connect_four_neighbors_only']=section
    sections['connect_using_avg_resistances']=section
    
    section='Short circuit regions (aka polygons)'
    sections['use_polygons']=section
    sections['polygon_file']=section
    
    section='Options for advanced mode'
    sections['source_file']=section
    sections['ground_file']=section
    sections['ground_file_is_resistances']=section
    sections['use_unit_currents']=section
    sections['use_direct_grounds']=section
    sections['remove_src_or_gnd']=section
    
    section='Calculation options'
    sections['solver']=section
    sections['print_timings']=section
    sections['low_memory_mode']=section    
    
    section='Output options'
    sections['output_file']=section
    sections['write_cur_maps']=section
    sections['write_cum_cur_map_only']=section
    sections['write_max_cur_maps']=section
    sections['log_transform_maps']=section
    sections['write_volt_maps']=section
    sections['compress_grids']=section
    sections['set_null_voltages_to_nodata']=section
    sections['set_null_currents_to_nodata']=section
    sections['set_focal_node_currents_to_zero']=section
    
    section='Mask file'
    sections['use_mask']=section
    sections['mask_file']=section  
    
    section='Options for pairwise and one-to-all and all-to-one modes'
    sections['use_included_pairs']=section
    sections['included_pairs_file']=section
    sections['point_file']=section
    sections['point_file_contains_polygons']=section
    
    section='Options for one-to-all and all-to-one modes'
    sections['use_variable_source_strengths']=section
    sections['variable_source_file']=section

    
    section='Habitat raster or graph'
    sections['habitat_file']=section
    sections['habitat_map_is_resistances']=section
    
    section="Circuitscape mode"
    sections['scenario']=section
    sections['data_type']=section

    if options['ground_file_is_resistances']=='not entered':
        options['ground_file_is_resistances'] = False
    if options['point_file_contains_polygons']=='not entered':
        options['point_file_contains_polygons'] = False
 
    for option in sections:
        try:
            config.add_section(sections[option])
        except:
            pass
    for option in sections:
        config.set(sections[option], option, options[option])

    f = open(configFile, 'w')
    config.write(f)
    f.close()
 

def setDefaultOptions():
    """Sets options to default values."""    
    options = {}
    options['data_type']='raster' 
    options['version']='unknown'
    options['low_memory_mode']=False
    options['scenario']='not entered'
    options['habitat_file']='(Browse for a habitat map file)'
    options['habitat_map_is_resistances']=True
    options['point_file']='(Browse for file with locations of focal points or areas)'
    options['point_file_contains_polygons']=False
    options['connect_four_neighbors_only']=True
    options['connect_using_avg_resistances']=True
    options['use_polygons']=False
    options['polygon_file']='(Browse for a short-circuit region file)'
    options['source_file']='(Browse for a current source file)'
    options['ground_file']='(Browse for a ground point file)'
    options['ground_file_is_resistances']=True
    options['use_unit_currents']=False
    options['use_direct_grounds']=False
    options['remove_src_or_gnd']='not entered'
    options['output_file']='(Choose a base name for output files)'
    options['write_cur_maps']=False
    options['write_cum_cur_map_only']=False
    options['log_transform_maps']=False
    options['write_volt_maps']=False
    options['solver']='cg+amg'
    options['compress_grids']=False
    options['print_timings']=False
    options['use_mask']=False
    options['mask_file']='None' 
    options['use_included_pairs']=False
    options['included_pairs_file']='None' 
    options['use_variable_source_strengths']=False
    options['variable_source_file']='None' 
    options['set_null_voltages_to_nodata']=True # Default must be false or will need to  re-do verify .ini files
    options['set_null_currents_to_nodata']=True # Default must be false or will need to re-do verify .ini files
    options['write_max_cur_maps']=False
    options['set_focal_node_currents_to_zero']=False
    
    return options

    
def checkOptions(options):
    """Checks to make sure sufficient options are passed to Circuitscape to 
    
    complete a run.
    
    """    
    if options['scenario']=='not entered':
        all_options_entered=False
        message = 'Please choose a scenario' 
        
    elif options['habitat_file']=='(Browse for a habitat map file)':
        all_options_entered=False
        message = 'Please choose a raster habitat map file'
	    
    elif options['habitat_map_is_resistances']=='not entered':
        all_options_entered=False
        message = 'Please choose a habitat data type'
            
    elif options['scenario']=='pairwise'and options['point_file']=='(Browse for a file with focal points or regions)':
        all_options_entered=False
        message = 'Please choose a focal node file'

    elif options['scenario']=='one-to-all'and options['point_file']=='(Browse for a file with focal points or regions)':
        all_options_entered=False
        message = 'Please choose a focal node file'
        
    elif options['connect_four_neighbors_only']=='not entered':
        all_options_entered=False
        message = 'Please choose a cell connection scheme'
            
    elif options['connect_using_avg_resistances']=='not entered':
        all_options_entered=False 
        message = 'Please choose a cell connection calculation'
            
    elif options['scenario']=='advanced'and options['source_file']=='(Browse for a current source file)':
        all_options_entered=False 
        message = 'Please enter a current source file'
            
    elif options['scenario']=='advanced'and options['ground_file']=='(Browse for a ground point file)':
        all_options_entered=False
        message = 'Ground point file does not exist!'
            
    elif options['scenario']=='advanced'and options['ground_file_is_resistances']=='not entered':
        all_options_entered=False 
        message = 'Please choose a ground data type'
            
    elif options['use_polygons']==True and options['polygon_file']=='(Browse for a short-circuit region file)':
        all_options_entered=False
        message = 'Please enter a short-circuit region file or uncheck this option'
            
    elif options['output_file']=='(Choose an output file name)':
        all_options_entered=False
        message = 'Please choose an output file name'
            
    else:
      all_options_entered=True
      message='None'

    return all_options_entered, message      

