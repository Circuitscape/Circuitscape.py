"""
__version__ = "$Revision: 1.7 $"
__date__ = "$Date: 2004/08/12 19:18:51 $"
"""

from PythonCard import dialog, model
from gui_options_rsrc import GUI_OPTIONS_RSRC

class OptionsWindow(model.CustomDialog):
    def __init__(self, parent):
        model.CustomDialog.__init__(self, parent, aDialogRsrc=GUI_OPTIONS_RSRC)
        self.parent = self.getParent()
        self.childOptions = self.parent.options #key- copies over options to child
        self.setChildWidgets()

    def setChildWidgets(self):
        idx = self.parent.OPTIONS_SCENARIO.index(self.parent.options.scenario)
        pairwise_enabled, advanced_enabled = self.parent.SCENARIO_PAIRWISE_ADVANCED[idx]
        is_pairwise_scenario = (self.childOptions.scenario == 'pairwise')
        networkEnabled = self.childOptions.data_type == 'network'
        rasterEnabled = not networkEnabled


        ## Calculation Options
        self.components.connectFourNBox.checked     = self.childOptions.connect_four_neighbors_only
        self.components.avgConductanceBox.checked   = not self.childOptions.connect_using_avg_resistances
        self.components.releaseMemBox.checked       = self.childOptions.preemptive_memory_release
        self.components.lowMemBox.checked           = self.childOptions.low_memory_mode
        self.components.unitSrcsBox.checked         = self.childOptions.use_unit_currents
        self.components.directGndsBox.checked       = self.childOptions.use_direct_grounds
        
        if self.childOptions.remove_src_or_gnd == 'not entered':
            self.childOptions.remove_src_or_gnd = 'keepall' # For backward compatibility- 'not entered' was old default
        idx = self.parent.OPTIONS_REMOVE_SOURCE_GROUND.index(self.childOptions.remove_src_or_gnd)
        self.components.rmvSrcGndChoice.SetSelection(idx) 

        ## Mapping Options
        self.components.compressGridsBox.checked    = self.childOptions.compress_grids
        self.components.logCurMapBox.checked        = self.childOptions.log_transform_maps
        self.components.cumMapBox.checked           = self.childOptions.write_cum_cur_map_only
        self.components.zeroCurrentsBox.checked     = self.childOptions.set_focal_node_currents_to_zero
        self.components.maxMapBox.checked           = self.childOptions.write_max_cur_maps
        
        ## Other Inputs        
        self.components.varSrcBox.checked           = self.childOptions.use_variable_source_strengths and pairwise_enabled and not is_pairwise_scenario
        self.components.maskBox.checked             = self.childOptions.use_mask and rasterEnabled
        self.components.polygonsBox.checked         = self.childOptions.use_polygons and rasterEnabled
        self.components.includePairsBox.checked     = self.childOptions.use_included_pairs and pairwise_enabled and is_pairwise_scenario

        self.components.polygonFile.text            = str(self.childOptions.polygon_file)
        if self.childOptions.included_pairs_file != None:
            self.components.includePairsFile.text   = str(self.childOptions.included_pairs_file)
        if self.childOptions.mask_file != None:
            self.components.maskFile.text           = str(self.childOptions.mask_file) 
        if self.childOptions.variable_source_file != None:
            self.components.varSrcFile.text         = str(self.childOptions.variable_source_file)

        self.enable_disable_child_widgets(pairwise_enabled, advanced_enabled)

    def enable_disable_child_widgets(self, pairwise_enabled, advancedEnabled):
        is_pairwise_scenario = (self.childOptions.scenario == 'pairwise')
        networkEnabled = self.childOptions.data_type == 'network'
        rasterEnabled = not networkEnabled
        
        self.components.connectFourNBox.enabled     = rasterEnabled
        self.components.avgConductanceBox.enabled   = rasterEnabled
        self.components.cumMapBox.enabled           = pairwise_enabled 
        self.components.zeroCurrentsBox.enabled     = pairwise_enabled and rasterEnabled
        self.components.maxMapBox.enabled           = pairwise_enabled
        self.components.lowMemBox.enabled           = pairwise_enabled and is_pairwise_scenario
        self.components.compressGridsBox.enabled    = rasterEnabled
        self.components.directGndsBox.enabled       = advancedEnabled
        self.components.unitSrcsBox.enabled         = advancedEnabled

        foreColor = self.parent.COLOR_ENABLED if advancedEnabled else self.parent.COLOR_DISABLED
        self.components.rmvSrcGndTitle.foregroundColor   = foreColor  
        self.components.rmvSrcGndChoice.enabled  = advancedEnabled  
        
        self.components.maskBox.enabled             = rasterEnabled
        self.components.maskFile.enabled            = rasterEnabled and self.components.maskBox.checked
        self.components.maskBrowse.enabled          = self.components.maskBox.checked

        self.components.polygonsBox.enabled         = rasterEnabled
        self.components.polygonFile.enabled         = rasterEnabled and self.components.polygonsBox.checked
        self.components.polygonBrowse.enabled       = self.components.polygonsBox.checked

        self.components.includePairsBox.enabled     = pairwise_enabled and is_pairwise_scenario
        self.components.includePairsFile.enabled    = pairwise_enabled and is_pairwise_scenario and self.components.includePairsBox.checked
        self.components.includePairsBrowse.enabled  = self.components.includePairsFile.enabled

        self.components.varSrcBox.enabled           = pairwise_enabled and not is_pairwise_scenario
        self.components.varSrcFile.enabled          = self.components.varSrcBox.checked
        self.components.varSrcBrowse.enabled        = self.components.varSrcBox.checked

    def on_maskBox_mouseClick(self, event):
        self.components.maskFile.enabled = self.components.maskBox.checked
        self.components.maskBrowse.enabled = self.components.maskBox.checked

    def on_polygonsBox_mouseClick(self, event):
        self.components.polygonFile.enabled = self.components.polygonsBox.checked
        self.components.polygonBrowse.enabled = self.components.polygonsBox.checked

    def on_includePairsBox_mouseClick(self, event):
        self.components.includePairsFile.enabled = self.components.includePairsBox.checked
        self.components.includePairsBrowse.enabled = self.components.includePairsBox.checked

    def on_varSrcBox_mouseClick(self, event):
        self.components.varSrcFile.enabled = self.components.varSrcBox.checked
        self.components.varSrcBrowse.enabled = self.components.varSrcBox.checked

    def on_loadPolygonBox_mouseClick(self, event):   
        print event.GetSelection()
        self.components.polygonFile.enabled = event.GetSelection()

    def on_polygonBrowse_mouseClick(self, event):
        wildcard = "ASCII Raster (*.asc)|*.asc|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Short-Circuit Region Raster', '', '', wildcard) 
        if result.accepted == True:        
            file_name = str(result.paths[0])
            self.components.polygonFile.text = file_name
            self.childOptions.polygon_file = file_name
        else:
            self.components.polygonFile.text = ''

    def on_maskBrowse_mouseClick(self, event):
        wildcard = "ASCII Raster (*.asc)|*.asc|All Files (*.*)|*.*"
        result = dialog.fileDialog(self, 'Select Raster Mask.  All cells with NODATA or non-positive-integer values will be dropped from resistance map. ', '', '', wildcard) 
        if result.accepted == True: 
            file_name = str(result.paths[0])
            self.components.maskFile.text = file_name
            self.childOptions.mask_file = file_name
        else:
            self.components.maskFile.text = ''

    def on_varSrcBrowse_mouseClick(self, event):
        wildcard = "Tab-delimited text list (*.txt)|*.txt|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select List of Source Strengths', '', '', wildcard ) 
        if result.accepted == True: 
            file_name = str(result.paths[0])
            self.components.varSrcFile.text = file_name
            self.childOptions.variable_source_file = file_name
        else:
            self.components.varSrcFile.text = ''

    def on_includePairsBrowse_mouseClick(self, event):
        wildcard = "Tab-delimited text list (*.txt)|*.txt|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Matrix of Focal Node Pairs to Include/Exclude ', '', '', wildcard ) 
        if result.accepted == True: 
            file_name = str(result.paths[0])
            self.components.includePairsFile.text = file_name
            self.childOptions.included_pairs_file = file_name
        else:
            self.components.includePairsFile.text = ''

    def get_options_in_child_boxes(self):
        self.childOptions.connect_four_neighbors_only       = self.components.connectFourNBox.checked
        self.childOptions.connect_using_avg_resistances     = not self.components.avgConductanceBox.checked
        self.childOptions.use_unit_currents                 = self.components.unitSrcsBox.checked
        self.childOptions.use_direct_grounds                = self.components.directGndsBox.checked

        self.childOptions.compress_grids                    = self.components.compressGridsBox.checked
        self.childOptions.log_transform_maps                = self.components.logCurMapBox.checked
        self.childOptions.write_cum_cur_map_only            = self.components.cumMapBox.checked
        self.childOptions.set_focal_node_currents_to_zero   = self.components.zeroCurrentsBox.checked
        self.childOptions.write_max_cur_maps                = self.components.maxMapBox.checked

        self.childOptions.preemptive_memory_release = self.components.releaseMemBox.checked
        self.childOptions.low_memory_mode = self.components.lowMemBox.checked
        self.childOptions.use_variable_source_strengths = self.components.varSrcBox.checked
        self.childOptions.use_mask = self.components.maskBox.checked 
        self.childOptions.use_polygons = self.components.polygonsBox.checked
        self.childOptions.use_included_pairs = self.components.includePairsBox.checked

        self.childOptions.polygon_file = self.components.polygonFile.text
        self.childOptions.included_pairs_file = self.components.includePairsFile.text
        self.childOptions.mask_file = self.components.maskFile.text 
        self.childOptions.variable_source_file = self.components.varSrcFile.text
        if 'browse for a' in self.childOptions.mask_file.lower():
            self.childOptions.use_mask = False
        if 'browse for a' in self.childOptions.polygon_file.lower():
            self.childOptions.use_polygons = False
        if 'browse for a' in self.childOptions.variable_source_file.lower():
            self.childOptions.use_variable_source_strengths = False
        if 'browse for a' in self.childOptions.included_pairs_file.lower():
            self.childOptions.use_included_pairs = False

    def on_rmvSrcGndChoice_select(self, event):
        self.childOptions.remove_src_or_gnd = self.parent.OPTIONS_REMOVE_SOURCE_GROUND[event.GetSelection()]    

def show_options_window(parent):
    dlg = OptionsWindow(parent)
    result = dlg.showModal()
    dlg.get_options_in_child_boxes()
    result.options = dlg.childOptions 
    dlg.destroy()
    return result



