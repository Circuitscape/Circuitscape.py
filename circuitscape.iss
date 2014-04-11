;ENTER VERSION INFO HERE
#define Version "4.0.4"  

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{E4376ABF-DB0A-41EA-B854-54ED7A20BD38}
AppName=Circuitscape
AppVerName=Circuitscape {#Version}
AppPublisher=Circuitscape
AppPublisherURL=http://www.circuitscape.org/
AppSupportURL=http://www.circuitscape.org/
AppUpdatesURL=http://www.circuitscape.org/
DefaultDirName={pf}\Circuitscape
DefaultGroupName=Circuitscape
OutputBaseFilename=setup
Compression=lzma
SolidCompression=yes
; For 64-bit only:
ArchitecturesInstallIn64BitMode=x64 

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked
Name: "quicklaunchicon"; Description: "{cm:CreateQuickLaunchIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "dist\csgui.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "dist\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs; Permissions: users-modify
; BHM added permissions above for 64-bit install 01/17/2011
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{group}\Circuitscape"; Filename: "{app}\csgui.exe"
Name: "{commondesktop}\Circuitscape"; Filename: "{app}\csgui.exe"; IconFilename: "{app}\cs_logo.ico"; Tasks: desktopicon
Name: "{userappdata}\Microsoft\Internet Explorer\Quick Launch\Circuitscape"; Filename: "{app}\csgui.exe"; IconFilename: "{app}\cs_logo.ico"; Tasks: quicklaunchicon

[Run]
;Filename: "{app}\circuitscape_4_0.pdf"; Description: "View the user guide"; Flags: postinstall shellexec skipifsilent
Filename: http://docs.circuitscape.org/circuitscape_4_0_user_guide.html?id=Win32Install{#Version}; Description: "View the user guide"; Flags: postinstall shellexec skipifsilent
Filename: "{app}\csgui.exe"; Description: "{cm:LaunchProgram,Circuitscape}"; Flags: nowait postinstall skipifsilent

