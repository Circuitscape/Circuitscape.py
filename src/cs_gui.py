##
## Circuitscape (C) 2008, 2009, 2010, Brad McRae and Viral B. Shah. 
##
## $Id: cs_gui.py 805 2012-07-30 23:11:04Z mcrae $
##


import imp, os, sys
import ConfigParser
import traceback
import time
import traceback

import wxversion
try:
    wxversion.select('2.8')
except:
    try:
        wxversion.select('2.7')
    except:
        pass
    
import wx, PythonCard
from PythonCard import dialog, model 
from PythonCard.components import button, checkbox, choice, image, staticline, statictext, textfield

from cs_util import *
from cs_compute import *
from cs_verify import *

class cs_gui(model.Background):
    def on_initialize(self, event):
        self.state = {}        
        self.state['version']='3.5.8'
        #LOAD LAST self.options
        configFile = 'circuitscape.ini'
        self.options = self.LoadOptions(configFile) 
        self.options['version'] = self.state['version']
        ##Set all objects to reflect options
        self.setWidgets()
        self.components.calcButton.SetFocus()
        self.statusBar = self.CreateStatusBar()
        self.statusBar.SetFieldsCount(3)        
        
        if self.options['data_type']=='network':
            self.enable_disable_network_widgets(True)        
            statustext=str('V ' + self.state['version']+' BETA NETWORK MODE')            
        else:
            statustext=str('Version ' + self.state['version']+' Ready.')
            self.enable_disable_network_widgets(False)            
        self.statusBar.SetStatusText(statustext,0)


           
    ##MENU ITEMS
    def on_menuFileLoadLast_select(self, event):
        configFile='circuitscape.ini'
        self.options=self.LoadOptions(configFile)
        self.options['version'] = self.state['version']
        ##Set all objects to reflect options
        if self.options['data_type']=='network':
            self.enable_disable_network_widgets(True)
            statustext=str('V ' + self.state['version']+' BETA NETWORK MODE')
        else:
            statustext=str('Version ' + self.state['version']+' Ready.')
            self.enable_disable_network_widgets(False)            
        self.statusBar.SetStatusText(statustext,0)
        
        self.setWidgets()
        self.components.calcButton.SetFocus()


    def on_menuFileLoadPrev_select(self, event):
        wildcard = "Options Files (*.ini)|*.ini|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Circuitscape Options File', '', '', wildcard ) 
        if result.accepted==True:        
            configFile = result.paths[0]
            self.options = self.LoadOptions(configFile)
            self.options['version'] = self.state['version']
            if self.options['data_type']=='network':
                statustext=str('V ' + self.state['version']+' BETA NETWORK MODE')            
                self.enable_disable_network_widgets(True)
            else:
                statustext=str('Version ' + self.state['version']+' Ready.')
                self.enable_disable_network_widgets(False)                            
            self.statusBar.SetStatusText(statustext,0)            
            #Set all gui objects to reflect options    
            self.setWidgets()
            
            print '\n\n'
            self.components.calcButton.SetFocus()

    def on_menuFileVerifyCode_select(self, event):
        print 'Verifying code (this will take a minute or two)'
        self.statusBar.SetStatusText('Verifying code (this will take a minute or two)',0)
        self.statusBar.SetStatusText('',1)        

        try:
            testResult=cs_verifyall()
            if testResult.wasSuccessful():
                testsPassed=True
            else:
                testsPassed=False
        except:
            testsPassed=False

        if testsPassed==True:
            dial = wx.MessageDialog(None, 'All tests passed!', 'Verification complete.', wx.OK)
            dial.ShowModal()
        else:
            dial = wx.MessageDialog(None, 'Errors were found.  Please see terminal or console for details.', 'Verification failed.', wx.OK)
            dial.ShowModal()
        if self.options['data_type']=='network':
            statustext=str('V ' + self.state['version']+' BETA NETWORK MODE')
            self.enable_disable_network_widgets(True)
        else:
            statustext=str('Version ' + self.state['version']+' Ready.')
            self.enable_disable_network_widgets(False)            
        self.statusBar.SetStatusText(statustext,0)

    def on_menuFileSave_select(self, event):
        self.components.habitatFile.SetFocus() #Need to force loseFocus on text boxes to make sure they are updated.
        self.components.outFile.SetFocus()
        self.components.calcButton.SetFocus()
        #Get options set in menu bar
        self.options['use_unit_currents']=self.menuBar.getChecked('menuOptionsUnitSrcs')
        self.options['use_direct_grounds']=self.menuBar.getChecked('menuOptionsDirectGnds')
        self.options['write_cum_cur_map_only']=self.menuBar.getChecked('menuOptionsCumMap')
        self.options['write_max_cur_maps']=self.menuBar.getChecked('menuOptionsMaxMap')
        self.options['low_memory_mode']=self.menuBar.getChecked('menuOptionsLowMemory')
        self.options['log_transform_maps']=self.menuBar.getChecked('menuOptionsLogCurMap')
        self.options['compress_grids']=self.menuBar.getChecked('menuOptionsCompressGrids')
        self.options['print_timings']=self.menuBar.getChecked('menuOptionsPrintTimings')
        self.options['use_mask']=self.menuBar.getChecked('menuOptionsMask')
        self.options['use_variable_source_strengths']=self.menuBar.getChecked('menuOptionsVarSrc')
        self.options['use_included_pairs']=self.menuBar.getChecked('menuOptionsIncludePairs')



        
        rmvgnd=self.menuBar.getChecked('menuOptionsRmvGnd')
        rmvsrc=self.menuBar.getChecked('menuOptionsRmvSrc')
        if rmvgnd==True:
            if rmvsrc==True:
                self.options['remove_src_or_gnd']='rmvall'
            else:
                self.options['remove_src_or_gnd']='rmvgnd'
        elif rmvsrc==True:   
            self.options['remove_src_or_gnd']='rmvsrc'
        else:
            self.options['remove_src_or_gnd']='keepall'

        wildcard = "Options Files (*.ini)|*.ini"  #FIXME: This block of code is redundant with later code saving options.  Roll into single function.
        wildcard = '*.ini'
        result = dialog.saveFileDialog(self, 'Choose a file name', '', '', wildcard)
        if result.accepted==True:                
            optionsFileName=result.paths[0]
            outputDir, outputFile = os.path.split(optionsFileName)
            outputBase, outputExtension = os.path.splitext(outputFile)

            configFile = outputDir + '//' + outputBase + '.ini' 
            if os.path.isdir(outputDir):
                writeConfigFile(configFile, self.options)
            else:
                message=str('Output directory ' + outputDir + ' does not exist!')
                dial = wx.MessageDialog(None, message, 'Error writing configuration file', wx.OK | wx.ICON_ERROR)
                dial.ShowModal()
        return

    def on_menuFileRunBatch_select(self, event):
        wildcard = "Options Files (*.ini)|*.ini" 
        result = dialog.fileDialog(self, 'Select any number of Circuitscape Options Files within one directory', '', '', wildcard ) 
        if result.accepted==True:
            wx.BeginBusyCursor()
            print '\nRunning Circuitscape in batch mode.\n'
            startTime=time.time()
            startTimeHMS = time.strftime('%H:%M:%S')
            self.statusBar.SetStatusText('Batch start ' + str(startTimeHMS),0)
            job=0
            numjobs=len(result.paths)
            for selection in result.paths:
                job=job+1
                configDir, configFile = os.path.split(selection)
                print 'Processing',configFile,'\n'
                self.statusBar.SetStatusText('Batch start ' + str(startTimeHMS) + '. Running job ' + str(job) +'/' +str(numjobs),0)
                try:
                    cs = cs_compute(selection,self.statusBar.SetStatusText)
                except RuntimeError, error:
                    message=str(error)
                    dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)
                    dial.ShowModal()
                    return
                except MemoryError:
                    self.memory_error_feedback()
                    return
                except:
                    self.unknown_exception()
                    return

                try:
                    result,solver_failed = cs.compute()
                    print 'Finished processing',configFile,'\n'
                except RuntimeError, error:
                    message=str(error)
                    dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)
                    dial.ShowModal()
                except MemoryError:
                    self.memory_error_feedback()
                    return
                except:
                    self.unknown_exception()
            print '\nDone with batch operations.\n'
            wx.EndBusyCursor()
            self.components.calcButton.SetFocus()
            if self.options['data_type']=='network':
                self.enable_disable_network_widgets(True)            
                statustext=str('V ' + self.state['version']+' BETA NETWORK MODE')
            else:
                statustext=str('Version ' + self.state['version']+' Ready.')
                self.enable_disable_network_widgets(False)            
            self.statusBar.SetStatusText(statustext,0)

            self.statusBar.SetStatusText('',1)
            (hours,mins,secs)=elapsed_time(startTime)
            if hours>0:
                self.statusBar.SetStatusText('Batch job took ' + str(hours) +' hours ' + str(mins) + ' minutes to complete.',2)
            else:
                self.statusBar.SetStatusText('Batch job took ' + str(mins) +' minutes ' + str(secs) + ' seconds to complete.',2)


            
    def on_menuOptionsCumMap_select(self, event):
        self.menuBar.setChecked('menuOptionsCumMap', self.menuBar.getChecked('menuOptionsCumMap'))

    def on_menuOptionsMaxMap_select(self, event):
        self.menuBar.setChecked('menuOptionsMaxMap', self.menuBar.getChecked('menuOptionsMaxMap'))

    def on_menuOptionsLowMemory_select(self, event):
        self.menuBar.setChecked('menuOptionsLowMemory', self.menuBar.getChecked('menuOptionsLowMemory'))

    def on_menuOptionsUnitSrcs_select(self, event):
        self.menuBar.setChecked('menuOptionsUnitSrcs',self.menuBar.getChecked('menuOptionsUnitSrcs'))


    def on_menuOptionsDirectGnds_select(self, event):
        self.menuBar.setChecked('menuOptionsDirectGnds',self.menuBar.getChecked('menuOptionsDirectGnds'))
    

    def on_menuOptionsLogCurMap_select(self, event):
        self.menuBar.setChecked('menuOptionsLogCurMap',self.menuBar.getChecked('menuOptionsLogCurMap'))


    def on_menuOptionsCompressGrids_select(self, event):
        self.menuBar.setChecked('menuOptionsCompressGrids', self.menuBar.getChecked('menuOptionsCompressGrids'))

    def on_menuOptionsPrintTimings_select(self, event):
        self.menuBar.setChecked('menuOptionsPrintTimings', self.menuBar.getChecked('menuOptionsPrintTimings'))

    def on_menuOptionsRmvGnd_select(self, event):
        self.menuBar.setChecked('menuOptionsRmvGnd',self.menuBar.getChecked('menuOptionsRmvGnd'))

    def on_menuOptionsRmvSrc_select(self, event):
        self.menuBar.setChecked('menuOptionsRmvSrc',self.menuBar.getChecked('menuOptionsRmvSrc'))


    def on_menuOptionsMask_select(self, event):
        self.menuBar.setChecked('menuOptionsMask',self.menuBar.getChecked('menuOptionsMask'))
        if self.menuBar.getChecked('menuOptionsMask') == True:
            wildcard = "ASCII Raster (*.asc)|*.asc|All Files (*.*)|*.*"
            result = dialog.fileDialog(self, 'Select Raster Mask.  All cells with NODATA or non-positive-integer values will be dropped from habitat map. ', '', '', wildcard ) 
            if result.accepted==True: 
                file = result.paths[0]
                self.options['mask_file']=file
      

    def on_menuOptionsVarSrc_select(self, event):
        self.menuBar.setChecked('menuOptionsVarSrc',self.menuBar.getChecked('menuOptionsVarSrc'))
        if self.menuBar.getChecked('menuOptionsVarSrc') == True:
            wildcard = "Tab-delimited text list (*.txt)|*.txt|All Files (*.*)|*.*" 
            result = dialog.fileDialog(self, 'Select List of Source Strengths', '', '', wildcard ) 
            if result.accepted==True: 
                file = result.paths[0]
                self.options['variable_source_file']=file


    def on_menuOptionsIncludePairs_select(self, event):
        self.menuBar.setChecked('menuOptionsIncludePairs',self.menuBar.getChecked('menuOptionsIncludePairs'))
        if self.menuBar.getChecked('menuOptionsIncludePairs') == True:
            wildcard = "Tab-delimited text list (*.txt)|*.txt|All Files (*.*)|*.*" 
            result = dialog.fileDialog(self, 'Select Matrix of Focal Node Pairs to Include/Exclude ', '', '', wildcard ) 
            if result.accepted==True: 
                file = result.paths[0]
                self.options['included_pairs_file']=file


    def on_menuFileAbout_select(self, event):
        messagetext=str('Version ' + self.state['version']+'\n\nhttp://www.circuitscape.org/\n\nBrad McRae and Viral B. Shah\n\nCircuitscape (C) 2008-09. Licensed under LGPL.')
        dial = wx.MessageDialog(None, messagetext, 'Circuitscape', wx.OK)
        dial.ShowModal()



    ##CHOICE BOXES
    def on_scenarioChoice_select(self, event):   
        scenario=event.GetSelection() 
        if scenario == 0:
            self.options['scenario']='not entered'
            self.enable_disable_widgets(0,0)	
        elif scenario == 1:
            self.options['scenario']='pairwise'
            self.enable_disable_widgets(1,0)	
        elif scenario == 2:
            self.options['scenario']='one-to-all'
            self.enable_disable_widgets(1,0)
        elif scenario == 3:
            self.options['scenario']='all-to-one'
            self.enable_disable_widgets(1,0)            
        else:
            self.options['scenario']='advanced'
            self.enable_disable_widgets(0,1)

    def on_habResistanceChoice_select(self, event):   
        hab = event.GetSelection()
        if hab==0:
            self.options['habitat_map_is_resistances']='not entered'
        elif hab==1:
            self.options['habitat_map_is_resistances']=True
        else:
            self.options['habitat_map_is_resistances']=False            
        

    def on_connCalcChoice_select(self, event):   
        calc = event.GetSelection()
        if calc==0:
            self.options['connect_using_avg_resistances']='not entered' 
        elif calc==1:
            self.options['connect_using_avg_resistances']=True
        else: 
            self.options['connect_using_avg_resistances']=False

    def on_connSchemeChoice_select(self, event):   
        scheme = event.GetSelection()
        if scheme == 0:
            self.options['connect_four_neighbors_only']='not entered'
        elif scheme == 1:
            self.options['connect_four_neighbors_only']=True
        else:
            self.options['connect_four_neighbors_only']=False

    def on_focalNodeChoice_select(self, event):
        choice = event.GetSelection() 
        if choice == 0:
            self.options['point_file_contains_polygons']='not entered'
        elif choice==1:
            self.options['point_file_contains_polygons']=False
        else:
            self.options['point_file_contains_polygons']=True

    def on_gndResistanceChoice_select(self, event):   
        gndResistance=event.GetSelection()
        if gndResistance==0:
            self.options['ground_file_is_resistances']='not entered'
        elif gndResistance==1:
            self.options['ground_file_is_resistances']=True
        else:
            self.options['ground_file_is_resistances']=False

             
##CHECK BOXES

    def on_loadPolygonBox_mouseClick(self, event):   
        self.options['use_polygons']=event.GetSelection() 
        if self.options['use_polygons']==True:
            self.components.polygonFile.enabled = True
            self.components.polygonBrowse.enabled = True
        else:            
            self.components.polygonFile.enabled = False
            self.components.polygonBrowse.enabled = False
            
    def on_curMapBox_mouseClick(self, event):   
        self.options['write_cur_maps']=event.GetSelection() 
   
    def on_voltMapBox_mouseClick(self, event):   
        self.options['write_volt_maps']=event.GetSelection() 
               
    
##BROWSE BUTTONS
    def on_habitatBrowse_mouseClick(self, event):
        if self.options['data_type']=='network':
            wildcard = "Tab-delimited text list (*.txt)|*.txt|All Files (*.*)|*.*" 
        else:
            wildcard = "ASCII Raster (*.asc)|*.asc|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Raster Habitat Map', '', '', wildcard ) 
        if result.accepted==True: 
            file = result.paths[0]
            self.components.habitatFile.text = file
            self.options['habitat_file']=file
      
    def on_srcTargetBrowse_mouseClick(self, event):
        wildcard = "Tab-delimited text list (*.txt) or ASCII Raster (*.asc)|*.txt;*.asc|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Source/Target File', '', '', wildcard ) 
        if result.accepted==True:        
            file = result.paths[0]
            self.components.srcTargetFile.text = file                    
            self.options['point_file']=file

    def on_currentSrcBrowse_mouseClick(self, event):
        wildcard = "Tab-delimited text list (*.txt) or ASCII Raster (*.asc)|*.txt;*.asc|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Source/Target File', '', '', wildcard ) 
        if result.accepted==True:        
            file = result.paths[0]
            self.components.currentSrcFile.text = file                    
            self.options['source_file']=file

    def on_polygonBrowse_mouseClick(self, event):
        wildcard = "ASCII Raster (*.asc)|*.asc|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Short-Circuit Region Raster', '', '', wildcard ) 
        if result.accepted==True:        
            file = result.paths[0]
            self.components.polygonFile.text = file
            self.options['polygon_file']=file
        else:
            self.components.polygonFile.text = ''
            
    def on_outBrowse_mouseClick(self, event):
        wildcard = "OUT Files (*.out)|*.out|All Files (*.*)|*.*"
        result = dialog.saveFileDialog(self, 'Choose a Base Output File Name', '', '', wildcard)
        if result.accepted==True:                
            file=result.paths[0]
            self.components.outFile.text = file
            self.options['output_file']=file

    def on_gndBrowse_mouseClick(self, event):
        wildcard = "Tab-delimited text list (*.txt) or ASCII Raster (*.asc)|*.txt;*.asc|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Ground File', '', '', wildcard ) 
        if result.accepted==True:        
            file = result.paths[0]
            self.components.gndFile.text = file
            self.options['ground_file']=file
        
##TEXT BOXES
    def on_habitatFile_loseFocus(self,event):
        self.options['habitat_file']=self.components.habitatFile.text        
        
    def on_srcTargetFile_loseFocus(self,event):
        self.options['point_file']=self.components.srcTargetFile.text 
            
    def on_currentSrcFile_loseFocus(self,event):
        self.options['source_file']=self.components.currentSrcFile.text
        
    def on_polygonFile_loseFocus(self,event):
        self.options['polygon_file']=self.components.polygonFile.text
        
    def on_outFile_loseFocus(self,event):
        self.options['output_file']=self.components.outFile.text
            
    def on_gndFile_loseFocus(self,event):
        self.options['ground_file']=self.components.gndFile.text

            
##CALCULATE    
    def on_calcButton_mouseClick(self, event):
        self.components.habitatFile.SetFocus() #Need to force loseFocus on text boxes to make sure they are updated.
        self.components.outFile.SetFocus()
        self.components.calcButton.SetFocus()
        #Check to see if all inputs are chosen
        (all_options_entered, message) = checkOptions(self.options)
        if all_options_entered==True:
            #Get options set in menu bar
            self.options['use_unit_currents']=self.menuBar.getChecked('menuOptionsUnitSrcs')
            self.options['use_direct_grounds']=self.menuBar.getChecked('menuOptionsDirectGnds')
            self.options['write_cum_cur_map_only']=self.menuBar.getChecked('menuOptionsCumMap')
            self.options['write_max_cur_maps']=self.menuBar.getChecked('menuOptionsMaxMap')            
            self.options['low_memory_mode']=self.menuBar.getChecked('menuOptionsLowMemory')            
            self.options['log_transform_maps']=self.menuBar.getChecked('menuOptionsLogCurMap')
            self.options['compress_grids']=self.menuBar.getChecked('menuOptionsCompressGrids')
            self.options['print_timings']=self.menuBar.getChecked('menuOptionsPrintTimings')
            self.options['use_mask']=self.menuBar.getChecked('menuOptionsMask')
            self.options['use_variable_source_strengths']=self.menuBar.getChecked('menuOptionsVarSrc')
            self.options['use_included_pairs']=self.menuBar.getChecked('menuOptionsIncludePairs')


            rmvgnd=self.menuBar.getChecked('menuOptionsRmvGnd')
            rmvsrc=self.menuBar.getChecked('menuOptionsRmvSrc')
            if rmvgnd==True:
                if rmvsrc==True:
                    self.options['remove_src_or_gnd']='rmvall'
                else:
                    self.options['remove_src_or_gnd']='rmvgnd'
            elif rmvsrc==True:   
                self.options['remove_src_or_gnd']='rmvsrc'
            else:
                self.options['remove_src_or_gnd']='keepall'
                                 
            #save selected options in local directory
            configFile='circuitscape.ini'
            writeConfigFile(configFile, self.options)

            #save copy of options under output path
            fileName = self.options['output_file']
            outputDir, outputFile = os.path.split(fileName)
            outputBase, outputExtension = os.path.splitext(outputFile)
            configFile = outputDir + '//' + outputBase + '.ini' 
            if os.path.isdir(outputDir):
                writeConfigFile(configFile, self.options)
            else:
                message=str('Output directory ' + outputDir + ' does not exist!')
                dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)
                dial.ShowModal()
                return  
                
            
            print '\nCalling Circuitscape...\n\n'
            startTime = time.strftime('%H:%M:%S')
            self.statusBar.SetStatusText('Job started ' + str(startTime),0)
            try:
                cs = cs_compute('circuitscape.ini',self.statusBar.SetStatusText)
            except RuntimeError, error:
                message=str(error)
                dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)
                dial.ShowModal()
                return

            except:
                self.unknown_exception()
                return

            try:
                terminate=self.checkHeaders()
                if terminate==True:
                    return
            except RuntimeError, error:
                message=str(error)
                dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)
                dial.ShowModal()
                return

            if self.options['data_type']=='network':        
                print '***Running in Network (Graph) Mode***'                
                # dial = wx.MessageDialog(None, 'Network mode activated.\nThis is beta code and '
                                        # 'requires raw network data input.','Note', wx.OK | wx.ICON_INFORMATION)
                # dial.ShowModal()
                
                
            if self.options['scenario']=='pairwise':
                try:
                    wx.BeginBusyCursor()

                    resistances,solver_failed = cs.compute()
                    wx.EndBusyCursor()
                    self.components.calcButton.SetFocus()

                    if solver_failed==True:
                        print '\nPairwise resistances (-1 indicates disconnected focal node pair, -777 indicates failed solve):'
                    else:
                        print '\nPairwise resistances (-1 indicates disconnected node pair):'
                    print resistances
                    print '\nDone.\n'
                    if solver_failed==True:
                        message='At least one solve failed.  Failure is coded as -777 in output resistance matrix.'
                        dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_EXCLAMATION)
                        dial.ShowModal()
                except MemoryError:
                    self.memory_error_feedback()
                    return
                except RuntimeError, error:
                    message=str(error)
                    dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)
                    dial.ShowModal()
                    return
                except:
                    self.unknown_exception()

            elif self.options['scenario']=='advanced':
                try:
                    wx.BeginBusyCursor()
                    voltages,solver_failed=cs.compute()
                    wx.EndBusyCursor()
                    self.components.calcButton.SetFocus()
                    if solver_failed==True:
                        message='Solver failed!'
                        dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_EXCLAMATION)
                        dial.ShowModal()
                    
                    print '\nDone.\n'
                except RuntimeError, error:
                    message=str(error)
                    dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)
                    dial.ShowModal()
                    return
                except MemoryError:
                    self.memory_error_feedback()
                    return

                except:
                    self.unknown_exception()
                    return

            else:
                try:
                    wx.BeginBusyCursor()
                    resistances,solver_failed = cs.compute()
                    wx.EndBusyCursor()
                    self.components.calcButton.SetFocus()
                    if self.options['scenario']=='all-to-one':
                        if solver_failed==True:
                            print '\nResult for each focal node \n(0 indicates successful calculation, -1 indicates disconnected node, -777 indicates failed solve):\n'
                        else:
                            print '\nResult for each focal node \n(0 indicates successful calculation, -1 indicates disconnected node):\n'
                    elif solver_failed==True:
                        print '\nResistances (-1 indicates disconnected node, -777 indicates failed solve):\n'
                    else:
                        print '\nResistances (-1 indicates disconnected node):\n'
                    print resistances
                    print '\nDone.\n'
                    if solver_failed==True:
                        message='At least one solve failed.  Failure is coded as -777 in output node/resistance list.'
                        dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_EXCLAMATION)
                        dial.ShowModal()

                except RuntimeError, error:
                    message=str(error)
                    dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)
                    dial.ShowModal()
                except MemoryError:
                    self.memory_error_feedback()
                    return
                except:
                    self.unknown_exception()
            if self.options['data_type']=='network':
                statustext=str('V ' + self.state['version']+' BETA NETWORK MODE')
            else:
                statustext=str('Version ' + self.state['version']+' Ready.')
            self.statusBar.SetStatusText(statustext,0)

            self.statusBar.SetStatusText('',1)

        else:
            dial = wx.MessageDialog(None, message, 'Not all options entered', wx.OK | wx.ICON_EXCLAMATION)
            dial.ShowModal()



    def checkHeaders(self):
        """Checks to make sure headers (with cell size, number of cols, etc) match for input rasters."""  
        if self.options['data_type']=='network':
            return
        (ncols, nrows, xllcorner, yllcorner, cellsize, nodata)=read_header(self.options['habitat_file'])
        headerMismatch=False
        terminate=False
        if self.options['use_polygons']==True: 
            (ncols2, nrows2, xllcorner2, yllcorner2, cellsize2, nodata2)=read_header(self.options['polygon_file'])                        
            if (ncols2!=ncols) or (nrows2!=nrows) or (abs(xllcorner2- xllcorner) > cellsize/3) or (abs(yllcorner2- yllcorner) > cellsize/3) or (cellsize2!=cellsize):
                headerMismatch=True

        filename=self.options['point_file']
        base, extension = os.path.splitext(filename)
        if extension == ".asc":        
            (ncols3, nrows3, xllcorner3, yllcorner3, cellsize3, nodata3)=read_header(self.options['point_file'])                                    
            if (ncols3!=ncols) or (nrows3!=nrows) or (abs(xllcorner3- xllcorner) > cellsize/3) or (abs(yllcorner3- yllcorner) > cellsize/3) or (cellsize3!=cellsize):
                headerMismatch=True             
        

        if self.options['use_mask']==True: 
            (ncols2, nrows2, xllcorner2, yllcorner2, cellsize2, nodata2)=read_header(self.options['mask_file'])                        
            if (ncols2!=ncols) or (nrows2!=nrows) or (abs(xllcorner2- xllcorner) > cellsize/3) or (abs(yllcorner2- yllcorner) > cellsize/3) or (cellsize2!=cellsize):
                headerMismatch=True        

        if headerMismatch==True:
            result = wx.MessageDialog(None, "Raster map headers do not match.  Circuitscape can try to resample maps to match the habitat map (Beta code, no guarantees). \n\nNote:all maps MUST be in the same projection.  Some focal nodes or short-circuit regions may be lost. \n\nUsing the 'Export to Circuitscape' ArcGIS tool is a better bet.\n\nContinue?", "Warning", wx.YES_NO).ShowModal()
            if result == wx.ID_YES:
                terminate=False

            if result == wx.ID_NO:
                terminate=True
        return terminate


##Error handling
    def unknown_exception(self):
        try:
            dial = wx.MessageDialog(None, 'An unknown error occurred.  Please see message in terminal.', 'Error', wx.OK | wx.ICON_ERROR)
            dial.ShowModal()
            type, value, tb = sys.exc_info()
            info = traceback.extract_tb(tb)
            print 'full traceback:'
            print info
            print '***************'
            filename, lineno, function, text = info[-1] # last line only
            print "\n %s:%d: %s: %s (in %s)" %\
                  (filename, lineno, type.__name__, str(value), function)
        finally:
            type = value = tb = None # clean up
 
    def memory_error_feedback(self):
        print '\nCircuitscape ran out of memory. \nPlease see user guide for information about memory requirements.'
        message='Circuitscape ran out of memory. \nPlease see user guide for information about memory requirements.'
        dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)
        dial.ShowModal()
        try:
            type, value, tb = sys.exc_info()
            info = traceback.extract_tb(tb)
            print 'full traceback:'
            print info
            print '***************'
            filename, lineno, function, text = info[-1] # last line only
            print "\n %s:%d: %s: %s (in %s)" %\
                  (filename, lineno, type.__name__, str(value), function)
        finally:
            type = value = tb = None # clean up

        
###SUBROUTINES
    def setWidgets(self):
        if self.options['scenario'] == 'not entered':
            self.enable_disable_widgets(0,0) #May want to disable both scenario options
            self.components.scenarioChoice.SetSelection(0)
        elif self.options['scenario']=='pairwise':
            self.enable_disable_widgets(1,0)	
            self.components.scenarioChoice.SetSelection(1)
        elif self.options['scenario']=='one-to-all':
            self.enable_disable_widgets(1,0)	
            self.components.scenarioChoice.SetSelection(2)            
        elif self.options['scenario']=='all-to-one':
            self.enable_disable_widgets(1,0)	
            self.components.scenarioChoice.SetSelection(3)   
        else:
            self.enable_disable_widgets(0,1)	
            self.components.scenarioChoice.SetSelection(4)
        self.components.habitatFile.text=self.options['habitat_file']
        self.components.srcTargetFile.text=self.options['point_file']
        if self.options['point_file_contains_polygons'] == False:
            self.components.focalNodeChoice.SetSelection(1)
        elif self.options['point_file_contains_polygons'] == True:
            self.components.focalNodeChoice.SetSelection(2)
        else:    
            self.components.focalNodeChoice.SetSelection(0)
            
        self.components.polygonFile.text=self.options['polygon_file']
        if self.options['use_polygons']==True:        
            self.components.polygonBrowse.enabled = True
            self.components.polygonFile.enabled = True
        else:
            self.components.polygonBrowse.enabled = False
            self.components.polygonFile.enabled = False
        self.components.currentSrcFile.text=self.options['source_file']
        self.components.gndFile.text=self.options['ground_file']
        self.components.outFile.text=self.options['output_file']
        if self.options['habitat_map_is_resistances']=='not entered':
            self.components.habResistanceChoice.SetSelection(0)
        elif self.options['habitat_map_is_resistances']==True:
            self.components.habResistanceChoice.SetSelection(1)
        else:
            self.components.habResistanceChoice.SetSelection(2)
        if self.options['ground_file_is_resistances']==True:
            self.components.gndResistanceChoice.SetSelection(1)
        elif self.options['ground_file_is_resistances']==False:
            self.components.gndResistanceChoice.SetSelection(2)
        else:
            self.components.gndResistanceChoice.SetSelection(0)

        if self.options['connect_four_neighbors_only']=='not entered':        
            self.components.connSchemeChoice.SetSelection(0)
        elif self.options['connect_four_neighbors_only']==True: 
            self.components.connSchemeChoice.SetSelection(1)
        else:
            self.components.connSchemeChoice.SetSelection(2)
        if self.options['connect_using_avg_resistances']=='not entered':
            self.components.connCalcChoice.SetSelection(0)
        elif self.options['connect_using_avg_resistances']==True:
            self.components.connCalcChoice.SetSelection(1)
        else:
            self.components.connCalcChoice.SetSelection(2)
            
        self.components.loadPolygonBox.checked=self.options['use_polygons']
        self.menuBar.setChecked('menuOptionsUnitSrcs', self.options['use_unit_currents'])
        self.components.curMapBox.checked=self.options['write_cur_maps']
        self.components.voltMapBox.checked=self.options['write_volt_maps']

        self.menuBar.setChecked('menuOptionsCumMap', self.options['write_cum_cur_map_only'])
        self.menuBar.setChecked('menuOptionsMaxMap', self.options['write_max_cur_maps'])
        self.menuBar.setChecked('menuOptionsLowMemory', self.options['low_memory_mode'])

        self.menuBar.setChecked('menuOptionsLogCurMap',self.options['log_transform_maps'])
        self.menuBar.setChecked('menuOptionsCompressGrids',self.options['compress_grids'])
        self.menuBar.setChecked('menuOptionsPrintTimings',self.options['print_timings'])
        self.menuBar.setChecked('menuOptionsUnitSrcs', self.options['use_unit_currents'])
        self.menuBar.setChecked('menuOptionsDirectGnds',self.options['use_direct_grounds'])
        self.menuBar.setChecked('menuOptionsMask', self.options['use_mask'])
        self.menuBar.setChecked('menuOptionsVarSrc', self.options['use_variable_source_strengths'])
        self.menuBar.setChecked('menuOptionsIncludePairs',self.options['use_included_pairs'])
        
        
        if self.options['remove_src_or_gnd']=='rmvall':
            self.menuBar.setChecked('menuOptionsRmvGnd', True)   
            self.menuBar.setChecked('menuOptionsRmvSrc', True)        

        elif self.options['remove_src_or_gnd']=='rmvgnd':
            self.menuBar.setChecked('menuOptionsRmvGnd', True)   
            self.menuBar.setChecked('menuOptionsRmvSrc', False)        

        elif self.options['remove_src_or_gnd']=='rmvsrc':
            self.menuBar.setChecked('menuOptionsRmvGnd', False)   
            self.menuBar.setChecked('menuOptionsRmvSrc', True)        
        else:        
            self.menuBar.setChecked('menuOptionsRmvGnd', False)   
            self.menuBar.setChecked('menuOptionsRmvSrc', False)   
            
        
    def enable_disable_network_widgets(self, networkEnabled):
        if networkEnabled == True:
            setting = False
            self.components.focalNodeChoice.SetSelection(1)            
        else:
            setting = True
        self.components.polygonFile.enabled = setting
        self.components.polygonBrowse.enabled = setting
        self.components.connSchemeChoice.enabled = setting
        self.components.connCalcChoice.enabled = setting
        self.components.loadPolygonBox.enabled = setting                                
        
                
    def enable_disable_widgets(self, pairwiseEnabled, advancedEnabled):     
        self.components.currentSrcFile.enabled = advancedEnabled
        self.components.currentSrcBrowse.enabled = advancedEnabled
        self.components.gndFile.enabled = advancedEnabled
        self.components.gndBrowse.enabled = advancedEnabled
        self.menuBar.setEnabled('menuOptionsUnitSrcs', advancedEnabled)
        self.menuBar.setEnabled('menuOptionsDirectGnds', advancedEnabled)  
        self.menuBar.setEnabled('menuOptionsRmvGnd', advancedEnabled)   
        self.menuBar.setEnabled('menuOptionsRmvSrc', advancedEnabled)
        self.components.gndResistanceChoice.enabled = advancedEnabled
        self.components.focalNodeChoice.enabled = pairwiseEnabled
        self.components.polygonBrowse.enabled = True
        self.components.srcTargetFile.enabled = pairwiseEnabled
        self.components.srcTargetBrowse.enabled = pairwiseEnabled
        self.menuBar.setEnabled('menuOptionsCumMap',pairwiseEnabled)    
        self.menuBar.setEnabled('menuOptionsMaxMap',pairwiseEnabled)    
        if self.options['scenario']=='pairwise':
            self.menuBar.setEnabled('menuOptionsLowMemory', pairwiseEnabled)
        else:
            self.menuBar.setEnabled('menuOptionsLowMemory', False)

        self.menuBar.setEnabled('menuOptionsIncludePairs', pairwiseEnabled) 
        self.menuBar.setEnabled('menuOptionsVarSrc', pairwiseEnabled)
        if self.options['scenario']=='pairwise':
            self.menuBar.setEnabled('menuOptionsVarSrc', False)


        if pairwiseEnabled==True:
            self.components.pairwiseOptionsTitle.foregroundColor=(0, 0, 160)
            self.components.srcTargetFileText.foregroundColor=(0, 0, 160)        
        else:
            self.components.pairwiseOptionsTitle.foregroundColor=(180,180,180)
            self.components.srcTargetFileText.foregroundColor=(180,180,180)        
        if advancedEnabled==True:
            self.components.gndFileText.foregroundColor=(0, 0, 160)
            self.components.advancedOptionsTitle.foregroundColor=(0, 0, 160)
            self.components.srcFileText.foregroundColor=(0, 0, 160)
        else:                     
            self.components.gndFileText.foregroundColor=(180,180,180)
            self.components.advancedOptionsTitle.foregroundColor=(180,180,180)
            self.components.srcFileText.foregroundColor=(180,180,180)
            
    def LoadOptions(self,configFile):
        """Sets options based on configuration file from last run or sets default options if no file exists."""  
        if os.path.isfile(configFile):
            try:
                options = readConfigFile(configFile)
            except:
                options = setDefaultOptions() 

        else:
            options=setDefaultOptions()

        return options

    
if __name__ == '__main__':
    app = model.Application(cs_gui)
    app.MainLoop()

