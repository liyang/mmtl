; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

[Setup]
AppName=TNT
AppVerName=TNT Transmission Line Suite 1.2.0
AppPublisher=Mayo Clinic
AppPublisherURL=http://mmtl.sourceforge.net/
AppSupportURL=http://mmtl.sourceforge.net/
AppUpdatesURL=http://mmtl.sourceforge.net/
DefaultDirName={pf}\tnt-1.2.0
DefaultGroupName=TNT
AllowNoIcons=yes
LicenseFile=O:\apps\tnt\1.2.0\win32\COPYING.txt
Compression=lzma
SolidCompression=yes
OutputDir=.
OutputBaseFilename=setup-tnt-1.2.0


[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "O:\apps\tnt\1.2.0\win32\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{group}\TNT 1.2.0"; Filename: "{app}\bin\tclkit.exe"; Parameters: """{app}\bin\tnt.tcl"""; WorkingDir: "{app}/examples"
Name: "{group}\Users Guide (pdf)"; Filename: "{app}\doc\user-guide.pdf"
Name: "{group}\Users Guide"; Filename: "{app}\doc\user-guide\index.html"
Name: "{group}\{cm:UninstallProgram,TNT}"; Filename: "{uninstallexe}"
Name: "{userdesktop}\TNT 1.2.0"; Filename: "{app}\bin\tclkit.exe"; Parameters: """{app}\bin\tnt.tcl"""; WorkingDir: "{app}/examples"; Tasks: desktopicon

[Run]
Filename: "{app}\bin\tclkit.exe"; Parameters: """{app}\bin\tnt.tcl"""; WorkingDir: "{app}/examples"; Description: "{cm:LaunchProgram,TNT}"; Flags: nowait postinstall skipifsilent
Filename: "{app}\README.txt"; Description: "View the README file"; Flags: postinstall shellexec skipifsilent

