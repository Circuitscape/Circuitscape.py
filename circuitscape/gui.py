# import sys
# if sys.modules.has_key('wx'):
    # print 'wx0'
import os, sys, traceback, logging, time, multiprocessing, tempfile, StringIO, webbrowser, string
import numpy as np

import wxversion
# print wxversion.getInstalled()

try:
    #wxversion.select(['2.9', '2.8', '2.7'])
    wxversion.select('2.8')
except:
    try:
        wxversion.select('3.0')
    except:
        try:
            wxversion.select('2.9')
        except:
            try:
                wxversion.select('2.7')
            except:
                pass

import wx
import wx.lib.newevent
from PythonCard import dialog, model
from PythonCard.components import button, checkbox, choice, image, staticline, statictext, textfield, spinner, textarea

from gui_options import show_options_window
from compute import Compute, CSIO, __version__
from compute_base import ComputeBase
from cfg import CSConfig
from gui_rsrc import GUI_RSRC

wxLogEvent, EVT_WX_LOG_EVENT = wx.lib.newevent.NewEvent()   # @UndefinedVariable
wxver = wx.VERSION[0] + 0.1*(wx.VERSION[1])                 # @UndefinedVariable

class GUILogger(logging.Handler):
    def __init__(self, dest=None):
        logging.Handler.__init__(self)
        self.dest = dest
        self.level = logging.DEBUG
        self.last_gui_yield_time = time.time()

    def flush(self):
        pass

    def emit(self, record):
        try:
            msg = self.format(record)
            msg = msg.strip('\r')
            status_msg = None
            if ('\n' not in msg and 'elapsed' not in msg): # dispay simple info messages in the status bar
                plain_msg = logging._defaultFormatter.format(record)
                if (record.levelname == 'INFO'):
                    status_msg = (plain_msg, 1)
                elif (record.levelname == 'DEBUG'):
                    status_msg = (plain_msg, 2) 
            
            evt = wxLogEvent(message=msg, levelname=record.levelname, status_msg=status_msg)
            wx.PostEvent(self.dest, evt) # @UndefinedVariable
        
            (hours,mins,secs) = ComputeBase.elapsed_time(self.last_gui_yield_time)
            if (secs > 0) or (mins > 0) or (hours > 0):
                self.last_gui_yield_time = time.time()
                wx.SafeYield(None, True)  # @UndefinedVariable
                wx.GetApp().Yield(True)  # @UndefinedVariable
            
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


class GUI(model.Background): 
    OPTIONS_SCENARIO                        = ['not entered', 'pairwise', 'advanced', 'one-to-all', 'all-to-one']
    SCENARIO_PAIRWISE_ADVANCED              = [    (0,0),        (1,0),       (0,1),        (1,0),       (1,0)  ]
    OPTIONS_DATATYPE                        = ['not entered', 'raster', 'network']
    OPTIONS_REMOVE_SOURCE_GROUND            = ['rmvsrc', 'rmvgnd', 'rmvall', 'keepall']
    OPTIONS_HABITAT_MAP_IS_RESISTANCES      = ['not entered', True, False]
    OPTIONS_CONNECT_USING_AVG_RESISTANCES   = [True, False] 
    OPTIONS_CONNECT_FOUR_NEIGHBORS_ONLY     = [True, False] 
    OPTIONS_GROUND_FILE_IS_RESISTANCES      = ['not entered', True, False]
    OPTIONS_LOG_LEVEL                       = ['DEBUG', 'INFO', 'WARN', 'ERROR']
    
    COLOR_ENABLED = (0, 0, 160)
    COLOR_DISABLED = (180,180,180)
    
    logger = None
    log_handler = None
    
    def on_initialize(self, event):
        self.state = {}
        self.state['version'] = __version__       
       
        #LOAD LAST self.options
        configFile = 'circuitscape.ini'
        self.options = self.LoadOptions(configFile) 
        self.options.version = self.state['version']
        self.options.log_level = 'INFO'
        
        ##Set all objects to reflect options
        if sys.platform.startswith('win'):
            _icon = wx.Icon('cs_logo.ico', wx.BITMAP_TYPE_ICO) # @UndefinedVariable
            self.SetIcon(_icon)        
        self.components.Image1.file = get_packaged_resource('cs_logo.jpg')
        self.setWidgets()
        self.components.calcButton.SetFocus()
        self.statusBar = self.CreateStatusBar()
        self.statusBar.SetFieldsCount(3)        
        self.reset_status_bar()
        
        GUI.log_handler = GUILogger(self)
        GUI.logger = ComputeBase._create_logger("circuitscape_gui", getattr(logging, self.options.log_level.upper()), None, False, GUI.log_handler)
        self.Bind(EVT_WX_LOG_EVENT, self.onLogEvent)    

    ##MENU ITEMS
    def on_menuFileLoadLast_select(self, event):
        configFile='circuitscape.ini'
        self.options = self.LoadOptions(configFile)
        self.options.version = self.state['version']
        self.reset_status_bar()
        ##Set all objects to reflect options
        self.setWidgets()
        self.components.calcButton.SetFocus()

    def on_menuFileLoadPrev_select(self, event):
        wildcard = "Options Files (*.ini)|*.ini|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Circuitscape Options File', '', '', wildcard ) 
        if result.accepted==True:        
            configFile = result.paths[0]
            self.options = self.LoadOptions(configFile)
            self.options.version = self.state['version']
            self.reset_status_bar()
            #Set all gui objects to reflect options    
            self.setWidgets()
            
            self.components.calcButton.SetFocus()

    def on_menuFileVerifyCode_select(self, event):
        GUI.logger.info('Verifying code (this will take a few seconds)')
        self.statusBar.SetStatusText('Verifying code (this will take a few seconds)',0)
        self.statusBar.SetStatusText('',1)

        outdir = None
        cwd = os.getcwd()
        try:
            # root_path = os.path.dirname(__file__)
            root_path = os.path.dirname(os.path.realpath(sys.argv[0]))
            outdir = tempfile.mkdtemp()
            if os.path.exists(os.path.join(root_path,'circuitscape')): # when called from csgui.exe 
                os.chdir(root_path)     
            elif os.path.exists(os.path.join(os.path.split(root_path)[0],'circuitscape')): #When called from csgui.py 
                os.chdir(os.path.split(root_path)[0])     # otherwise we are running inside a packaged folder and resources are availale at cwd
            from verify import cs_verifyall
            testResult = cs_verifyall(out_path=outdir)
            testsPassed = testResult.wasSuccessful()
        except:
            testsPassed = False
        finally:
            os.chdir(cwd)
            if None != outdir:
                for root, dirs, files in os.walk(outdir, topdown=False):
                    for name in files:
                        os.remove(os.path.join(root, name))
                    for name in dirs:
                        os.rmdir(os.path.join(root, name))
                os.rmdir(outdir)

        if testsPassed:
            dial = wx.MessageDialog(None, 'All tests passed!', 'Verification complete.', wx.OK)  # @UndefinedVariable
            dial.ShowModal()
            GUI.logger.info('All tests passed!')
        else:
            dial = wx.MessageDialog(None, 'Errors were found.  Please see terminal or console for details.', 'Verification failed.', wx.OK)  # @UndefinedVariable
            dial.ShowModal()
        self.reset_status_bar()

    def on_menuHelpAbout_select(self, event):
        messagetext = str('Version ' + self.state['version'] + '\n\nhttp://www.circuitscape.org/\n\nBrad McRae, Viral B. Shah, and Tanmay K. Mohapatra\n\nCircuitscape (C) 2008-09. Licensed under LGPL.')
        dial = wx.MessageDialog(None, messagetext, 'Circuitscape', wx.OK)  # @UndefinedVariable
        dial.ShowModal()
        self.reset_status_bar()
        
    def on_menuOptionsMoreOptions_select(self, event):
        result = show_options_window(self)
        if result.accepted == True:
            self.options = result.options
            self.report_menu_files()
        self.reset_status_bar()
        
    def on_menuHelpManual_select(self, event):
        #TODO: change to match release tag
        webbrowser.open("http://docs.circuitscape.org/circuitscape_4_0_user_guide.html?id=Desktop_v" + self.state['version'])

    def on_menuHelpFeedback_select(self, event):
        webbrowser.open("http://www.circuitscape.org/Support?id=Desktop_v" + self.state['version'])

    def save_file_dlg(self, title, wildcard):
        file_name = None
        if wxver >= 2.9:  # use wx dialogs directly as python card save dialogs are broken after wx 2.9
            dlg = wx.FileDialog(self, title, '', '', wildcard, wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) # @UndefinedVariable
            if dlg.ShowModal() != wx.ID_CANCEL: # @UndefinedVariable
                file_name = dlg.GetPath()
        else:
            result = dialog.saveFileDialog(self, 'Choose a file name', '', '', wildcard)
            if result.accepted == True:
                file_name = result.paths[0]
        self.reset_status_bar()
        return file_name
        
    def on_menuFileSave_select(self, event):
        self.components.habitatFile.SetFocus() #Need to force loseFocus on text boxes to make sure they are updated.
        self.components.outFile.SetFocus()
        self.components.calcButton.SetFocus()
        
        options_file_name = self.save_file_dlg('Choose a file name', '*.ini')
        if options_file_name != None:                
            try:
                self.options.write(options_file_name, True)
            except RuntimeError as ex:
                message = str(ex)
                dial = wx.MessageDialog(None, message, 'Error writing configuration file', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
                dial.ShowModal()
        self.reset_status_bar()

    def on_menuFileRunBatch_select(self, event):
        wildcard = "Options Files (*.ini)|*.ini" 
        result = dialog.fileDialog(self, 'Select any number of Circuitscape Options Files within one directory', '', '', wildcard ) 
        if result.accepted==True:
            wx.BeginBusyCursor()  # @UndefinedVariable
            GUI.logger.debug('Running Circuitscape in batch mode')
            startTime = time.time()
            startTimeHMS = time.strftime('%H:%M:%S')
            self.statusBar.SetStatusText('Batch start ' + str(startTimeHMS), 0)
            job = 0
            numjobs = len(result.paths)
            for selection in result.paths:
                job += 1
                _configDir, configFile = os.path.split(selection)
                GUI.logger.debug('Processing ' + configFile)
                self.statusBar.SetStatusText('Batch start ' + str(startTimeHMS) + '. Running job ' + str(job) +'/' + str(numjobs), 0)
                
                try:
                    cs = Compute(selection, self)
                except RuntimeError as error:
                    wx.EndBusyCursor()  # @UndefinedVariable
                    message = str(error)
                    dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
                    dial.ShowModal()
                    return
                except MemoryError:
                    wx.EndBusyCursor()  # @UndefinedVariable
                    self.memory_error_feedback()
                    return
                except:
                    wx.EndBusyCursor()  # @UndefinedVariable
                    self.unknown_exception()
                    return

                try:
                    self.statusBar.SetStatusText('',1)
                    self.statusBar.SetStatusText('',2)
                    result, _solver_failed = cs.compute()
                    GUI.logger.debug('Finished processing ' + configFile)
                except RuntimeError as error:
                    message = str(error)
                    dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
                    dial.ShowModal()
                except MemoryError:
                    self.memory_error_feedback()
                    return
                except:
                    self.unknown_exception()
                    
            GUI.logger.debug('Done with batch operations.')
            wx.EndBusyCursor()  # @UndefinedVariable
            
            self.components.calcButton.SetFocus()
            self.reset_status_bar()

            (hours,mins,secs) = ComputeBase.elapsed_time(startTime)
            if hours > 0:
                self.statusBar.SetStatusText('Batch job took ' + str(hours) +' hours ' + str(mins) + ' minutes to complete.',2)
            else:
                self.statusBar.SetStatusText('Batch job took ' + str(mins) +' minutes ' + str(secs) + ' seconds to complete.',2)


    ##CHOICE BOXES
    def on_dataTypeChoice_select(self, event):
        data_type0 = self.options.data_type
        data_type = event.GetSelection()
        self.options.data_type = GUI.OPTIONS_DATATYPE[data_type]        
        if self.options.data_type != data_type0: # If selection has changed
            self.options.scenario = GUI.OPTIONS_SCENARIO[0]
            self.setWidgets()
        networkEnabled = self.options.data_type == 'network'
        self.enable_disable_network_widgets(networkEnabled)      
        self.report_menu_files()
    
    def on_scenarioChoice_select(self, event):   
        scenario = event.GetSelection()
        self.options.scenario               = GUI.OPTIONS_SCENARIO[scenario]
        pairwise_enabled, advanced_enabled  = GUI.SCENARIO_PAIRWISE_ADVANCED[scenario]
        self.enable_disable_widgets(pairwise_enabled, advanced_enabled)
        self.report_menu_files()

    def on_scenarioChoiceNetwork_select(self, event):   
        scenario = event.GetSelection()
        self.options.scenario               = GUI.OPTIONS_SCENARIO[scenario]
        pairwise_enabled, advanced_enabled  = GUI.SCENARIO_PAIRWISE_ADVANCED[scenario]
        self.enable_disable_widgets(pairwise_enabled, advanced_enabled)
        self.report_menu_files()

    def on_logLevelChoice_select(self, event):
        log_lvl = event.GetSelection()
        self.options.log_level = GUI.OPTIONS_LOG_LEVEL[log_lvl]
        GUI.logger = ComputeBase._create_logger("circuitscape_gui", getattr(logging, self.options.log_level), None, False, GUI.log_handler)
        GUI.logger.setLevel(getattr(logging, self.options.log_level.upper()))
        GUI.log_handler.setLevel(getattr(logging, self.options.log_level.upper()))

##CHECK BOXES   
    def on_useConductancesBox_mouseClick(self, event):   
        self.options.habitat_map_is_resistances = not event.GetSelection() 
         
    def on_useGroundConductancesBox_mouseClick(self, event):   
        self.options.ground_file_is_resistances = not event.GetSelection() 

    def on_curMapBox_mouseClick(self, event):   
        self.options.write_cur_maps = event.GetSelection() 

    def on_voltMapBox_mouseClick(self, event):   
        self.options.write_volt_maps = event.GetSelection() 

    def on_logRusageBox_mouseClick(self, event):
        self.options.print_rusages = event.GetSelection()

    def on_printTimingsBox_mouseClick(self, event):   
        self.options.print_timings  = event.GetSelection()
    
##BROWSE BUTTONS

    def on_habitatBrowse_mouseClick(self, event):
        if self.options.data_type == 'network':
            wildcard = "Tab-delimited text list (*.txt)|*.txt;*.txt.gz|All Files (*.*)|*.*" 
        else:
            wildcard = "ASCII Raster (*.asc or *.txt)|*.asc;*.txt;*.txt.gz;*.asc.gz|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Raster Resistance Map or Network/Graph File', '', '', wildcard) 
        if result.accepted == True: 
            file_name = result.paths[0]
            self.components.habitatFile.text = file_name
            self.options.habitat_file = file_name
      
    def on_srcTargetBrowse_mouseClick(self, event):
        wildcard = "Tab-delimited text list or ASCII Raster (*.txt or *.asc)|*.txt;*.asc;*.txt.gz;*.asc.gz|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Source/Target File', '', '', wildcard) 
        if result.accepted == True:        
            file_name = result.paths[0]
            self.components.srcTargetFile.text = file_name                    
            self.options.point_file = file_name

    def on_currentSrcBrowse_mouseClick(self, event):
        if self.options.data_type == 'network':
            wildcard = "Tab-delimited text list (*.txt)|*.txt;*.txt.gz|All Files (*.*)|*.*" 
        else:
            wildcard = "Tab-delimited text list or ASCII Raster (*.txt or *.asc)|*.txt;*.asc;*.txt.gz;*.asc.gz|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Source/Target File', '', '', wildcard) 
        if result.accepted == True:        
            file_name = result.paths[0]
            self.components.currentSrcFile.text = file_name                    
            self.options.source_file = file_name
           
    def on_outBrowse_mouseClick(self, event):
        wildcard = "OUT Files (*.out)|*.out|All Files (*.*)|*.*"
        file_name = self.save_file_dlg('Choose a Base Output File Name', wildcard)
        if file_name != None:
            self.components.outFile.text = file_name
            self.options.output_file = file_name

    def on_gndBrowse_mouseClick(self, event):
        if self.options.data_type == 'network':
            wildcard = "Tab-delimited text list (*.txt)|*.txt;*.txt.gz|All Files (*.*)|*.*" 
        else:
            wildcard = "Tab-delimited text list or ASCII Raster (*.txt or *.asc)|*.txt;*.asc;*.txt.gz;*.asc.gz|All Files (*.*)|*.*" 
        result = dialog.fileDialog(self, 'Select Ground File', '', '', wildcard ) 
        if result.accepted == True:        
            file_name = result.paths[0]
            self.components.gndFile.text = file_name
            self.options.ground_file = file_name
        
##TEXT BOXES
    def on_habitatFile_loseFocus(self,event):
        self.options.habitat_file = self.components.habitatFile.text        
        
    def on_srcTargetFile_loseFocus(self,event):
        self.options.point_file = self.components.srcTargetFile.text 
            
    def on_currentSrcFile_loseFocus(self,event):
        self.options.source_file = self.components.currentSrcFile.text
        
    def on_outFile_loseFocus(self,event):
        self.options.output_file = self.components.outFile.text
            
    def on_gndFile_loseFocus(self,event):
        self.options.ground_file = self.components.gndFile.text

    def on_parallelSpin_textUpdate(self, event):
        if self.components.parallelSpin.value > 1:
            self.options.max_parallel = self.components.parallelSpin.value
            self.options.parallelize = True
        else:
            self.options.max_parallel = 1
            self.options.parallelize = False
            
##CALCULATE
    def on_calcButton_mouseClick(self, event):
        self.components.habitatFile.SetFocus() #Need to force loseFocus on text boxes to make sure they are updated.
        self.components.outFile.SetFocus()
        self.components.calcButton.SetFocus()
        
        out_base, _out_ext = os.path.splitext(self.options.output_file)

        #if (self.options.log_file == None) or (not os.path.exists(self.options.log_file)): 
        self.options.log_file = out_base + '.log' # For now, just write log file to output directory.

        if (self.options.print_rusages and not sys.platform.startswith('win')) or (self.options.print_timings): 
            self.options.profiler_log_file = out_base + '_rusages.log' #For now, always write log file to output directory.
        else:
            self.options.profiler_log_file = None

        #Check to see if all inputs are chosen
        (all_options_entered, message) = self.options.check()
        
        if not all_options_entered:
            dial = wx.MessageDialog(None, message, 'Not all options entered', wx.OK | wx.ICON_EXCLAMATION)  # @UndefinedVariable
            dial.ShowModal()
            self.reset_status_bar()
            return
                                        
        #save selected options in local directory
        configFile = 'circuitscape.ini'
        self.options.write(configFile)
        try:
            self.options.write(self.options.output_file, True)
        except RuntimeError as ex:
            dial = wx.MessageDialog(None, str(ex), 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
            dial.ShowModal()
            self.reset_status_bar()
            return  
                            
        GUI.logger.info('Calling Circuitscape...')
        startTime = time.strftime('%H:%M:%S')
        self.statusBar.SetStatusText('Job started ' + str(startTime), 0)
        try:
            cs = Compute('circuitscape.ini', GUI.log_handler)
        except RuntimeError as error:
            message = str(error)
            dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
            dial.ShowModal()
            self.reset_status_bar()
            return
        except:
            self.unknown_exception()
            self.reset_status_bar()
            return

        try:
            filetype = CSIO._guess_file_type(self.options.habitat_file)
            if self.options.data_type == 'network': 
                if filetype != CSIO.FILE_TYPE_TXTLIST:
                    raise RuntimeError('Error reading network file "' + self.options.habitat_file + '".\nPlease check file format.')
            elif filetype != CSIO.FILE_TYPE_AAGRID and filetype != CSIO.FILE_TYPE_NPY:
                    raise RuntimeError('File "' + self.options.habitat_file + '"\ndoes not appear to be a raster. Please check file format.')
            terminate = self.checkHeaders()
            if terminate == True:
                self.reset_status_bar()
                return
        except RuntimeError as error:
            message = str(error)
            dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
            dial.ShowModal()
            self.reset_status_bar()
            return
        except:
            self.unknown_exception()
            self.reset_status_bar()
            return
        if self.options.data_type == 'network':        
            GUI.logger.debug('Running in Network (Graph) Mode')
                            
        if self.options.scenario == 'pairwise':
            try:
                wx.BeginBusyCursor()  # @UndefinedVariable
                self.statusBar.SetStatusText('',1)
                self.statusBar.SetStatusText('',2)
                resistances, solver_failed = cs.compute()
                if solver_failed == True:
                    msg = '\nPairwise resistances (-1 indicates disconnected focal node pair, -777 indicates failed solve):'
                else:
                    msg = '\nPairwise resistances (-1 indicates disconnected node pair):\nNode1\tNode2\tResistance'
                
                resistances3ColumnsStr = self.formatResistanceOutput(resistances)
                GUI.logger.info(msg + "\n" + resistances3ColumnsStr)
                GUI.logger.info('Done.\n')
                if solver_failed == True:
                    message = 'At least one solve failed.  Failure is coded as -777 in output resistance matrix.'
                    dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_EXCLAMATION)  # @UndefinedVariable
                    dial.ShowModal()
                wx.EndBusyCursor()  # @UndefinedVariable
                self.components.calcButton.SetFocus()
            except MemoryError:
                wx.EndBusyCursor()  # @UndefinedVariable
                self.memory_error_feedback()
                self.reset_status_bar()
                return
            except RuntimeError as error:
                wx.EndBusyCursor()  # @UndefinedVariable
                message = str(error)
                dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
                dial.ShowModal()
                self.reset_status_bar()
                return
            except:
                wx.EndBusyCursor()  # @UndefinedVariable
                self.unknown_exception()
        elif self.options.scenario == 'advanced':
            if self.options.write_cur_maps == False and self.options.write_volt_maps == False:
                message = 'Advanced mode selected but no output maps checked.\nThere is nothing to do.'
                dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
                dial.ShowModal()
                return

            wx.BeginBusyCursor()  # @UndefinedVariable

            try:
                self.statusBar.SetStatusText('',1)
                self.statusBar.SetStatusText('',2) 
                _voltages, solver_failed = cs.compute()
                wx.EndBusyCursor()  # @UndefinedVariable
                
                self.components.calcButton.SetFocus()
                
                if solver_failed == True:
                    message = 'Solver failed!'
                    dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_EXCLAMATION)  # @UndefinedVariable
                    dial.ShowModal()
                
                GUI.logger.info('Done.\n')
            except RuntimeError as error:
                wx.EndBusyCursor()  # @UndefinedVariable
                message = str(error)
                dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
                dial.ShowModal()
                self.reset_status_bar()
                return
            except MemoryError:
                wx.EndBusyCursor()  # @UndefinedVariable
                self.memory_error_feedback()
                self.reset_status_bar()
                return
            except:
                wx.EndBusyCursor()  # @UndefinedVariable
                self.unknown_exception()
                self.reset_status_bar()
                return

        else:
            try:
                wx.BeginBusyCursor()  # @UndefinedVariable
                self.statusBar.SetStatusText('',1)
                self.statusBar.SetStatusText('',2)                                    
                resistances, solver_failed = cs.compute()
                wx.EndBusyCursor()  # @UndefinedVariable
                
                self.components.calcButton.SetFocus()
                
                if self.options.scenario == 'all-to-one':
                    if solver_failed == True:
                        msg = '\nResult for each focal node (0 indicates successful calculation, -1 indicates disconnected node, -777 indicates failed solve):'
                    else:
                        msg = '\nResult for each focal node (0 indicates successful calculation, -1 indicates disconnected node):'
                elif solver_failed == True:
                    msg = '\nResistances to ground(-1 indicates disconnected node, -777 indicates failed solve):\n      Node         Resistance'
                else:
                    msg = '\nResistances to ground (-1 indicates disconnected node):\nNode\tResistance'
                resistanceString = self.formatResistanceOutput(resistances)
                GUI.logger.info(msg + '\n' + resistanceString)
                GUI.logger.info('Done.\n')
                
                if solver_failed == True:
                    message = 'At least one solve failed.  Failure is coded as -777 in output node/resistance list.'
                    dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_EXCLAMATION)  # @UndefinedVariable
                    dial.ShowModal()
            except RuntimeError as error:
                wx.EndBusyCursor()  # @UndefinedVariable
                message = str(error)
                dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
                dial.ShowModal()
                self.reset_status_bar()
            except MemoryError:
                wx.EndBusyCursor()  # @UndefinedVariable
                self.memory_error_feedback()
                self.reset_status_bar()
                return
            except:
                wx.EndBusyCursor()  # @UndefinedVariable
                self.unknown_exception()  
                self.reset_status_bar()
                return
            # self.reset_status_bar()        

    def formatResistanceOutput(self, resistances):
        """Converts resistances from matrix to 3-column string for printing."""  
        text = ''
        if resistances.shape[1] == 2: # One-to-all or all-to-one result
            numPoints = resistances.shape[0]
            for i in range(1,numPoints):
                text = text + str(int(resistances[i,0])) + '\t' + str(int(resistances[i,1])) + '\n'
            return text

        if self.options.data_type == 'network': # Already three column
            numpairs = resistances.shape[0]
            for i in range(0,numpairs):
                text = text + str(int(resistances[i,0])) + '\t' + str(int(resistances[i,1])) + '\t' + str(resistances[i,2]) + '\n'
            return text

        numPoints = resistances.shape[0]-1
        for i in range(1,numPoints):
            for j in range(i+1,numPoints+1):
                text = (text + str(int(resistances[i,0])) +'\t' + str(int(resistances[0,j])) 
                       +'\t' + str(resistances[i,j]) + '\n') 
        return text
        

    def checkHeaders(self):
        """Checks to make sure headers (with cell size, number of cols, etc) match for input rasters."""  
        if self.options.data_type == 'network':
            return
        match_files = []
        if self.options.use_polygons == True:
            match_files.append(self.options.polygon_file)
        if 'advanced' not in self.options.scenario.lower():
            filetype = CSIO._guess_file_type(self.options.point_file)
            if filetype == CSIO.FILE_TYPE_AAGRID:
                match_files.append(self.options.point_file)
        else:
            filetype = CSIO._guess_file_type(self.options.ground_file)
            if filetype == CSIO.FILE_TYPE_AAGRID:
                match_files.append(self.options.ground_file)
            filetype = CSIO._guess_file_type(self.options.source_file)
            if filetype == CSIO.FILE_TYPE_AAGRID:
                match_files.append(self.options.source_file)
        if self.options.use_mask == True:
            match_files.append(self.options.mask_file)
        headers_match = CSIO.match_headers(self.options.habitat_file, match_files)

        if headers_match == False:
            result = wx.MessageDialog(None, "Raster map headers do not match.  Circuitscape can try to resample maps to match the resistance map (Beta code, no guarantees). \n\nNote:all maps MUST be in the same projection.  Some focal nodes or short-circuit regions may be lost. \n\nUsing the 'Export to Circuitscape' tool or the Circuitscape for ArcGIS toolbox is a better bet.\n\nContinue?", "Warning", wx.YES_NO).ShowModal()  # @UndefinedVariable
            return (result == wx.ID_NO) # @UndefinedVariable
        return False


##Error handling
    def unknown_exception(self):
        try:
            dial = wx.MessageDialog(None, 'An unknown error occurred.  Please see message in terminal.', 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
            dial.ShowModal()
            strio = StringIO.StringIO()
            traceback.print_exc(file=strio)
            GUI.logger.error(strio.getvalue())
            strio.close()
        finally:
            e_type = value = tb = None # clean up
 

    def memory_error_feedback(self):
        GUI.logger.error('Circuitscape ran out of memory. Please see user guide for information about memory requirements.')
        message='Circuitscape ran out of memory. \nPlease see user guide for information about memory requirements.'
        dial = wx.MessageDialog(None, message, 'Error', wx.OK | wx.ICON_ERROR)  # @UndefinedVariable
        dial.ShowModal()
        try:
            e_type, value, tb = sys.exc_info()
            strio = StringIO.StringIO()
            info = traceback.extract_tb(tb)
            strio.write('full traceback:\n')
            for stack_line in info:
                strio.write(str(stack_line) + '\n')
            strio.write('\n***************\n')
            filename, lineno, function, _text = info[-1] # last line only
            strio.write("\n %s:%d: %s: %s (in %s)\n" % (filename, lineno, e_type.__name__, str(value), function))
            GUI.logger.error(strio.getvalue())
            strio.close()
        finally:
            e_type = value = tb = None # clean up

        
###SUBROUTINES
    def setWidgets(self):
        idx = GUI.OPTIONS_DATATYPE.index(self.options.data_type)
        self.components.dataTypeChoice.SetSelection(idx) 
        networkEnabled = self.options.data_type == 'network'
        self.enable_disable_network_widgets(networkEnabled) 
        
        idx = GUI.OPTIONS_SCENARIO.index(self.options.scenario)
        self.components.scenarioChoice.SetSelection(idx)
        pairwise_enabled, advanced_enabled = GUI.SCENARIO_PAIRWISE_ADVANCED[idx]
        self.enable_disable_widgets(pairwise_enabled, advanced_enabled)
        if idx > 2:
            self.components.scenarioChoiceNetwork.SetSelection(0) # Only pairwise and advanced modes are available for Network data type
        else:
            self.components.scenarioChoiceNetwork.SetSelection(idx)
        
        self.components.habitatFile.text = self.options.habitat_file
        self.components.srcTargetFile.text = self.options.point_file       
        if self.options.source_file is not None:
            self.components.currentSrcFile.text = self.options.source_file
        if self.options.ground_file is not None:
            self.components.gndFile.text = self.options.ground_file
        self.components.outFile.text = self.options.output_file
       
        idx = GUI.OPTIONS_LOG_LEVEL.index(self.options.log_level)
        self.components.logLevelChoice.SetSelection(idx)
        if GUI.logger != None:
            GUI.logger.setLevel(getattr(logging, self.options.log_level.upper()))

        if GUI.log_handler != None:
            GUI.log_handler.setLevel(getattr(logging, self.options.log_level.upper()))

        self.components.parallelSpin.value = self.components.parallelSpin.max = 1
        if sys.platform.startswith('win'):
            self.options.parallelize = False
        else:
            self.components.parallelSpin.max = n_cpus = multiprocessing.cpu_count()
            
            if self.options.max_parallel == 1:
                self.options.parallelize = False
                
            if self.options.parallelize:
                if self.options.max_parallel == 0:
                    self.options.max_parallel = n_cpus
                else:
                    self.options.max_parallel = min(self.options.max_parallel, n_cpus)
                self.components.parallelSpin.value = self.options.max_parallel
            
        
        self.components.useGroundConductancesBox.checked    = not self.options.ground_file_is_resistances
        self.components.useConductancesBox.checked          = not self.options.habitat_map_is_resistances
        self.components.curMapBox.checked                   = self.options.write_cur_maps
        self.components.voltMapBox.checked                  = self.options.write_volt_maps
        self.components.logRusageBox.checked                = self.options.print_rusages and not sys.platform.startswith('win')
        self.components.printTimingsBox.checked             = self.options.print_timings
        
        self.report_menu_files()
        
    def enable_disable_network_widgets(self, networkEnabled):        
        rasterEnabled = not networkEnabled       
        self.components.scenarioChoice.enabled = rasterEnabled
        self.components.scenarioChoice.visible = rasterEnabled # The two scenario choice menus are on top of one another. Only the relevant one is visible.
        self.components.scenarioChoiceNetwork.enabled = networkEnabled
        self.components.scenarioChoiceNetwork.visible = networkEnabled # The two scenario choice menus are on top of one another. Only the relevant one is visible.      
        
    def enable_disable_widgets(self, pairwiseEnabled, advancedEnabled):
        is_pairwise_scenario = (self.options.scenario == 'pairwise')
        
        self.components.currentSrcFile.enabled      = advancedEnabled
        self.components.currentSrcBrowse.enabled    = advancedEnabled
        self.components.gndFile.enabled             = advancedEnabled
        self.components.gndBrowse.enabled           = advancedEnabled

        self.components.useGroundConductancesBox.enabled     = advancedEnabled
        self.components.srcTargetFile.enabled               = pairwiseEnabled
        self.components.srcTargetBrowse.enabled             = pairwiseEnabled
        
        self.components.parallelSpin.enabled = is_pairwise_scenario and self.options.parallelize == True
        self.components.logRusageBox.enabled = not sys.platform.startswith('win') # Resource module not availble on windows
        foreColor = GUI.COLOR_ENABLED if pairwiseEnabled else GUI.COLOR_DISABLED
        self.components.pairwiseOptionsTitle.foregroundColor    = \
            self.components.srcTargetFileText.foregroundColor       = foreColor
        
        foreColor = GUI.COLOR_ENABLED if (is_pairwise_scenario and self.options.parallelize == True) else GUI.COLOR_DISABLED
        self.components.parallelizeText.foregroundColor = foreColor
                    
        foreColor = GUI.COLOR_ENABLED if advancedEnabled else GUI.COLOR_DISABLED
        self.components.gndFileText.foregroundColor                  = \
            self.components.advancedOptionsTitle.foregroundColor     = \
            self.components.srcFileText.foregroundColor              = \
            self.components.useGroundConductancesBox.foregroundColor = foreColor


    def LoadOptions(self, config_file):
        """Sets options based on configuration file from last run or sets default options if no file exists."""
        options = CSConfig()
        try:
            options = CSConfig(config_file)
        except:
            pass
        return options   

    def on_clearLogsButton_mouseClick(self, event):
        self.components.logMessages.clear()
        self.reset_status_bar()
        self.report_menu_files()

    def onLogEvent(self, event):
        self.components.logMessages.appendText(event.message + '\n')
        if event.status_msg != None:
            self.statusBar.SetStatusText(event.status_msg[0], event.status_msg[1])
            if event.status_msg[1] == 1:
                self.statusBar.SetStatusText('', 2)
        event.Skip()

    def report_menu_files(self):
        self.components.logMessages.clear()
        msg = ''
        if self.options.data_type != 'network':
            if self.options.use_mask == True:
                msg = ' '*20 + 'Mask file set to ' + str(self.options.mask_file) + '.\n'
            if self.options.use_polygons == True:
                msg = msg + ' '*20 + 'Short-circuit region file set to ' + str(self.options.polygon_file) + '.\n'
            if (self.options.use_variable_source_strengths == True) and (self.options.scenario == 'one-to-all' or self.options.scenario == 'all-to-one'):
                msg = msg + ' '*20 + 'Variable source strength file set to ' + str(self.options.variable_source_file) + '.\n'
        if self.options.use_included_pairs == True and self.options.scenario == 'pairwise':
            msg = msg + ' '*20 + 'Include/exclude file set to ' + str(self.options.included_pairs_file) + '.\n'
        if len(msg) > 1:
            self.components.logMessages.appendText(' '*20 + '-'*150 + '\n')
            self.components.logMessages.appendText(msg)
            self.components.logMessages.appendText(' '*20 + '-'*42 + 'These inputs can be changed in the Options menu' + '-'*42 +'\n\n')

    def reset_status_bar(self):
        statustext=str('Version ' + self.state['version'] + ' Ready.')
        self.statusBar.SetStatusText(statustext,0)
        self.statusBar.SetStatusText('', 1)
        self.statusBar.SetStatusText('Please send feedback to the Circuitscape User Group', 2)
        

def get_packaged_resource(filename):
    # coding: ascii
    exec_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
    res_dir = os.path.join(os.path.split(exec_dir)[0],'circuitscape') 
    res_path = os.path.join(res_dir, filename) # we are running inside packaged app. resources packaged separately, accessible at cwd
    if not os.path.exists(res_path):
        res_dir = os.path.join(exec_dir,'circuitscape') # called from executable
        res_path = os.path.join(res_dir, filename) 
    if not os.path.exists(res_path):
        res_path = filename
    return res_path
    
def show_gui():
    app = model.Application(GUI, rsrc=GUI_RSRC)
    if hasattr(app, 'startingDirectory'):
        os.chdir(app.startingDirectory)
    app.MainLoop()
                
if __name__ == '__main__':
    show_gui()
